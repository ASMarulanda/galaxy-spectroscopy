#!/usr/bin/env python3
"""Retrieve SDSS spectra by identifier, verify them, and convert them to cases.

    python tools/download_sdss.py --plate 294 --mjd 51986 --fiber 36 --out data/raw
    gspec measure --cases data/raw/cases --out data/raw/measurements.csv

Scope is deliberately narrow: **identifiers only**. There is no catalogue query
and therefore no selection function to declare. You say which spectrum you want;
the tool fetches that spectrum.

The archive path is

    https://data.sdss.org/sas/dr18/prior-surveys/sdss2-dr8-sdss/spectro/redux/26/
            spectra/{plate:04d}/spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits

Plate and fibre are zero-padded to four digits in *both* the directory and the
file name, and MJD to five. (DATA.md quoted the directory as ``{plate}``
unpadded; that path returns 404 for plate 294, and the padded one returns 200.
The padded form is the one this tool builds.)

What the tool guarantees
------------------------
* identifiers are integers inside documented ranges before any URL exists, so a
  string can never reach the request;
* the download lands in ``<name>.part`` and is renamed only after a size check and
  a SHA-256 check, so a file with the final name is always a complete file;
* an existing file is never replaced unless ``--force`` says so;
* retries are bounded (3 attempts, exponential backoff), on 5xx and on network
  errors, never on 404;
* every object gets a provenance sidecar, and the run gets a manifest, recording
  what came from the server and what came from you.

Exit codes
----------
0  every object fetched (and converted)
1  usage error, including a bad identifier
2  verification failed (size, SHA-256) — nothing measured for that object
3  the archive could not be reached after the retries
4  the file violated the instrument contract (conversion refused)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

BASE = ("https://data.sdss.org/sas/dr18/prior-surveys/sdss2-dr8-sdss/"
        "spectro/redux/26/spectra")

#: User agent. Some archives block the default urllib string.
USER_AGENT = "gspec/0.3 (+https://github.com/ASMarulanda/galaxy-spectroscopy)"

#: Identifier ranges, from the SDSS Legacy survey's own plate and MJD span. Wide
#: enough not to reject a real object, narrow enough to catch a swapped column.
PLATE_RANGE = (1, 9999)
MJD_RANGE = (40000, 69999)
FIBER_RANGE = (1, 1000)

CHUNK = 1 << 16


class IdentifierError(ValueError):
    """An identifier is not a plausible plate/MJD/fibre triple."""


class FetchError(Exception):
    """A fetch failed. ``kind`` selects the exit code.

    kind: notfound | network | verify | fits | usage
    """

    EXIT = {"usage": 1, "verify": 2, "network": 3, "notfound": 1, "fits": 4,
            "restart": 3}

    def __init__(self, message: str, kind: str = "network"):
        self.kind = kind
        super().__init__(message)

    @property
    def exit_code(self) -> int:
        return self.EXIT.get(self.kind, 3)


class RangeIgnored(FetchError):
    """The server answered a ranged request with a whole body."""

    def __init__(self, message: str = "the archive ignored the Range request"):
        super().__init__(message, kind="restart")


# --------------------------------------------------------------- identifiers
def validate_identifier(value, name: str, bounds: tuple[int, int]) -> int:
    """Return ``value`` as an int inside ``bounds`` or raise IdentifierError.

    ``int()`` is applied first, so any string that is not an integer literal is
    rejected before it can be interpolated into a URL. Bools are rejected by name
    because ``int(True) == 1`` would otherwise pass silently.
    """
    if isinstance(value, bool):
        raise IdentifierError(f"{name} must be an integer, not a boolean")
    if isinstance(value, int):
        ival = value
    else:
        text = str(value).strip()
        if not text or not (text.lstrip("+-").isdigit()):
            raise IdentifierError(f"{name} must be a plain integer; got {value!r}")
        ival = int(text)
    lo, hi = bounds
    if not (lo <= ival <= hi):
        raise IdentifierError(f"{name}={ival} is outside the documented range "
                              f"{lo}-{hi}")
    return ival


def validate_triple(plate, mjd, fiber) -> tuple[int, int, int]:
    return (validate_identifier(plate, "plate", PLATE_RANGE),
            validate_identifier(mjd, "mjd", MJD_RANGE),
            validate_identifier(fiber, "fiber", FIBER_RANGE))


def filename_for(plate: int, mjd: int, fiber: int) -> str:
    p, m, f = validate_triple(plate, mjd, fiber)
    return f"spec-{p:04d}-{m:05d}-{f:04d}.fits"


def url_for(plate: int, mjd: int, fiber: int) -> str:
    p, m, f = validate_triple(plate, mjd, fiber)
    return f"{BASE}/{p:04d}/{filename_for(p, m, f)}"


def case_id_for(plate: int, mjd: int, fiber: int) -> str:
    p, m, f = validate_triple(plate, mjd, fiber)
    return f"SDSS-{p:04d}-{m:05d}-{f:04d}"


# -------------------------------------------------------------------- http
@dataclass
class HeadInfo:
    status: int
    content_length: int | None = None
    etag: str | None = None
    last_modified: str | None = None
    accept_ranges: str | None = None


class UrllibOpener:
    """The real network. Tests substitute an object with the same two methods."""

    def __init__(self, read_timeout: float = 60.0):
        # urllib applies one socket timeout to connect and to each read.
        self.timeout = float(read_timeout)

    def _request(self, url: str, headers: dict | None = None):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                   **(headers or {})})
        return urllib.request.urlopen(req, timeout=self.timeout)

    def head(self, url: str) -> HeadInfo:
        req = urllib.request.Request(url, method="HEAD",
                                     headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return HeadInfo(status=resp.status,
                                content_length=_int_or_none(resp.headers.get("Content-Length")),
                                etag=resp.headers.get("ETag"),
                                last_modified=resp.headers.get("Last-Modified"),
                                accept_ranges=resp.headers.get("Accept-Ranges"))
        except urllib.error.HTTPError as exc:
            return HeadInfo(status=exc.code)

    def get(self, url: str, start: int = 0):
        headers = {"Range": f"bytes={start}-"} if start else {}
        try:
            return self._request(url, headers)
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and start:          # nothing left to send
                raise FetchError(f"range {start}- not satisfiable", kind="verify")
            raise


def _int_or_none(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _with_retries(fn, attempts: int = 3, base_delay: float = 1.0,
                  sleep=time.sleep, jitter=random.uniform):
    """Run ``fn`` with bounded exponential backoff on transport failures.

    404 is not retried: the object either exists at that identifier or it does
    not, and repeating the question wastes the archive's time.
    """
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except FetchError:
            raise
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise FetchError(f"HTTP 404 at {exc.url}", kind="notfound") from exc
            last = exc
            if exc.code < 500:
                raise FetchError(f"HTTP {exc.code} at {exc.url}", kind="network") from exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
        if attempt < attempts:
            delay = base_delay * (2 ** (attempt - 1)) * (1 + 0.25 * jitter(0, 1))
            sleep(delay)
    raise FetchError(f"gave up after {attempts} attempts: {last}", kind="network")


# ------------------------------------------------------------------ fetch
@dataclass
class FetchResult:
    case_id: str
    plate: int
    mjd: int
    fiber: int
    url: str
    status: str                      # fetched | present | resumed | failed
    path: str | None = None
    bytes: int | None = None
    sha256: str | None = None
    http: dict = field(default_factory=dict)
    retrieved_utc: str | None = None
    error: str | None = None
    case: dict | None = None

    def as_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}


def sha256_of(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_one(plate, mjd, fiber, out_dir: str | Path, *,
              opener=None, expected_sha256: str | None = None,
              force: bool = False, resume: bool = True,
              attempts: int = 3, sleep=time.sleep) -> FetchResult:
    """Fetch one spectrum to ``out_dir``. Returns a FetchResult; raises FetchError."""
    p, m, f = validate_triple(plate, mjd, fiber)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    opener = opener or UrllibOpener()
    url = url_for(p, m, f)
    dest = out_dir / filename_for(p, m, f)
    part = dest.with_name(dest.name + ".part")
    result = FetchResult(case_id=case_id_for(p, m, f), plate=p, mjd=m, fiber=f, url=url,
                         status="failed")

    if dest.exists() and not force:
        if expected_sha256 is not None:
            have = sha256_of(dest)
            if have.lower() != expected_sha256.lower():
                raise FetchError(
                    f"{dest.name} exists but its SHA-256 does not match the expected "
                    f"value ({have[:16]}... vs {expected_sha256[:16]}...); "
                    f"refusing to overwrite without --force", kind="verify")
        result.status = "present"
        result.path = str(dest)
        result.bytes = dest.stat().st_size
        result.sha256 = sha256_of(dest)
        return result

    head = _with_retries(lambda: _head(opener, url), attempts=attempts, sleep=sleep)
    if head.status == 404:
        raise FetchError(f"no spectrum at plate {p}, MJD {m}, fibre {f}: {url}",
                         kind="notfound")
    if head.status >= 400:
        raise FetchError(f"HTTP {head.status} for {url}", kind="network")
    expected = head.content_length
    start = 0
    if resume and part.exists():
        start = part.stat().st_size
        if expected is not None and start > expected:
            part.unlink()
            start = 0
    if head.accept_ranges is None and start:
        start = 0                                  # the server cannot resume

    mode = "ab" if start else "wb"
    written = start
    if start:
        result.status = "resumed"
    with open(part, mode) as fh:
        offset = start
        got_any = False
        while True:
            try:
                resp, chunk = _read_chunk(opener, url, offset, attempts=attempts, sleep=sleep)
            except RangeIgnored:
                # The archive cannot resume after all. Start over, once: the part
                # file is ours, nothing is lost but the bytes we already had.
                fh.close()
                part.unlink(missing_ok=True)
                if start == 0:
                    raise
                return fetch_one(p, m, f, out_dir, opener=opener,
                                 expected_sha256=expected_sha256, force=force,
                                 resume=False, attempts=attempts, sleep=sleep)
            if not chunk:
                break
            fh.write(chunk)
            written += len(chunk)
            offset += len(chunk)
            got_any = True
            if expected is not None and written >= expected:
                break
        if not got_any and start == 0:
            raise FetchError(f"the archive returned no bytes for {url}", kind="network")
        fh.flush()
        os.fsync(fh.fileno())

    if expected is not None and written != expected:
        kept = part
        raise FetchError(
            f"incomplete download: {written} bytes of {expected} expected; "
            f"kept {kept.name} for inspection", kind="verify")

    digest = sha256_of(part)
    if expected_sha256 is not None and digest.lower() != expected_sha256.lower():
        part.unlink(missing_ok=True)
        raise FetchError(f"SHA-256 mismatch for {dest.name}: got {digest}, expected "
                         f"{expected_sha256}; the file was deleted", kind="verify")

    os.replace(part, dest)
    result.status = "resumed" if start else "fetched"
    result.path = str(dest)
    result.bytes = written
    result.sha256 = digest
    result.retrieved_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    result.http = {"status": head.status, "etag": head.etag,
                   "last_modified": head.last_modified,
                   "content_length": expected, "accept_ranges": head.accept_ranges,
                   "bytes_written": written, "resumed_from": start}
    return result


def _head(opener, url):
    info = opener.head(url)
    if info.status >= 500:
        raise urllib.error.HTTPError(url, info.status, "server error", None, None)
    return info


def _read_chunk(opener, url, offset, attempts=3, sleep=time.sleep):
    def once():
        resp = opener.get(url, offset)
        if getattr(resp, "status", 200) == 200 and offset:
            raise RangeIgnored()
        try:
            return resp, resp.read(CHUNK)
        finally:
            close = getattr(resp, "close", None)
            if callable(close):
                close()
    return _with_retries(once, attempts=attempts, sleep=sleep)


# ------------------------------------------------------------- list input
def read_list(path: str | Path) -> list[tuple[int, int, int]]:
    """Read plate/mjd/fiber triples from a CSV, matched by column name.

    Header names are matched case-insensitively; ``fibre`` and ``fiberID`` are
    accepted for the fibre column, because the SDSS catalogue calls it the latter.
    """
    path = Path(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise FetchError(f"{path} has no header row", kind="usage")
        lower = {name.strip().lower(): name for name in reader.fieldnames}
        def pick(*names):
            for n in names:
                if n in lower:
                    return lower[n]
            raise FetchError(f"{path} needs a column named one of {names}; "
                             f"it has {reader.fieldnames}", kind="usage")
        cp, cm, cf = pick("plate"), pick("mjd"), pick("fiber", "fibre", "fiberid")
        triples = []
        for i, row in enumerate(reader, start=2):
            try:
                triples.append(validate_triple(row[cp], row[cm], row[cf]))
            except IdentifierError as exc:
                raise FetchError(f"{path} line {i}: {exc}", kind="usage") from exc
    return triples


# -------------------------------------------------------------------- main
def _convert(raw: Path, cases_dir: Path, provenance_dir: Path, force: bool) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fits_to_case

    out = cases_dir / f"{raw.stem}.fits"
    record = fits_to_case.convert(raw, out, force=force)
    prov_path = provenance_dir / f"{raw.stem}.conversion.json"
    prov_path.parent.mkdir(parents=True, exist_ok=True)
    prov_path.write_text(json.dumps(record, indent=2) + "\n")
    return record


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Retrieve SDSS spectra by identifier and convert them to gspec cases.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Exit codes: 0 ok, 1 usage, 2 verification failed, 3 network, "
               "4 instrument contract violated.")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--plate", type=str, help="plate number (with --mjd and --fiber)")
    g.add_argument("--list", help="CSV with plate, mjd, fiber columns")
    ap.add_argument("--mjd", type=str, help="MJD of the observation")
    ap.add_argument("--fiber", type=str, help="fibre number")
    ap.add_argument("--out", default="data/raw", help="download directory (default data/raw)")
    ap.add_argument("--convert", dest="convert", action="store_true", default=True,
                    help="convert to the case schema (default)")
    ap.add_argument("--no-convert", dest="convert", action="store_false",
                    help="download only")
    ap.add_argument("--measure", action="store_true",
                    help="also run gspec measure on the converted cases")
    ap.add_argument("--force", action="store_true", help="replace existing files")
    ap.add_argument("--no-resume", action="store_true", help="ignore .part files")
    ap.add_argument("--attempts", type=int, default=3, help="attempts per request (default 3)")
    ap.add_argument("--timeout", type=float, default=60.0, help="read timeout, seconds")
    ap.add_argument("--sha256", default=None, help="expected SHA-256 of the archive file")
    args = ap.parse_args(argv)

    try:
        if args.list:
            triples = read_list(args.list)
        else:
            if args.mjd is None or args.fiber is None:
                ap.error("--plate needs --mjd and --fiber")
            triples = [validate_triple(args.plate, args.mjd, args.fiber)]
    except IdentifierError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except FetchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return exc.exit_code

    out = Path(args.out)
    # --out is the workspace: archive files land directly in it, with the
    # converted cases and the provenance beside them. This is the layout the
    # tutorial and the notebook both use.
    raw_dir = out
    cases_dir = out / "cases"
    prov_dir = out / "provenance"
    for d in (raw_dir, prov_dir) + ((cases_dir,) if args.convert else ()):
        d.mkdir(parents=True, exist_ok=True)

    opener = UrllibOpener(read_timeout=args.timeout)
    results: list[FetchResult] = []
    worst = 0
    for plate, mjd, fiber in triples:
        cid = case_id_for(plate, mjd, fiber)
        print(f"{cid}: {url_for(plate, mjd, fiber)}", flush=True)
        try:
            r = fetch_one(plate, mjd, fiber, raw_dir, opener=opener,
                          expected_sha256=args.sha256, force=args.force,
                          resume=not args.no_resume, attempts=args.attempts)
            digest = (r.sha256 or "")[:16]
            print(f"  {r.status}: {r.bytes} bytes  sha256 {digest}...", flush=True)
            if args.convert and r.path:
                case_path = cases_dir / f"{Path(r.path).stem}.fits"
                if case_path.exists() and not args.force:
                    # Re-running a fetch must be a no-op, not an error.
                    r.case = {"path": str(case_path), "sha256": sha256_of(case_path),
                              "status": "present"}
                    print(f"  case: present  sha256 {r.case['sha256'][:16]}...", flush=True)
                else:
                    rec = _convert(Path(r.path), cases_dir, prov_dir, args.force)
                    r.case = {"path": rec["case_path"], "sha256": rec["case"]["sha256"],
                              "status": "converted"}
                    print(f"  case: {rec['case']['sha256'][:16]}...  "
                          f"R(Halpha)~{rec['fits'].get('resolution', {}).get('R_halpha')}",
                          flush=True)
            if r.path:
                sidecar = prov_dir / f"{Path(r.path).stem}.provenance.json"
                sidecar.write_text(json.dumps(r.as_dict(), indent=2) + "\n")
        except (FetchError, ValueError, OSError) as exc:
            kind = getattr(exc, "kind", "fits")
            r = FetchResult(case_id=cid, plate=plate, mjd=mjd, fiber=fiber,
                            url=url_for(plate, mjd, fiber), status="failed",
                            error=str(exc))
            print(f"  FAILED ({kind}): {exc}", file=sys.stderr, flush=True)
            worst = max(worst, FetchError(str(exc), kind).exit_code)
        results.append(r)

    manifest = {
        "tool": "tools/download_sdss.py",
        "archive": BASE,
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "objects": [r.as_dict() for r in results],
        "ok": sum(1 for r in results if r.status != "failed"),
        "failed": sum(1 for r in results if r.status == "failed"),
    }
    (prov_dir / "fetch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    if args.measure and args.convert:
        try:
            from gspec.cli import measure_directory
            print()
            measure_directory(cases_dir, out / "measurements.csv")
        except ImportError:
            print("error: --measure needs the gspec package (pip install -e .)",
                  file=sys.stderr)
            return 1

    print(f"\n{manifest['ok']} fetched, {manifest['failed']} failed; "
          f"manifest in {prov_dir / 'fetch_manifest.json'}")
    return worst


if __name__ == "__main__":
    sys.exit(main())
