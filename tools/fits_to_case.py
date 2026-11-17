#!/usr/bin/env python3
"""Convert an SDSS spectrum file into the case schema ``gspec`` reads.

The validated pipeline is not touched by this program. It adapts the input to the
pipeline, never the other way round: a case file is a FITS with a small PRIMARY
header and one binary table of six float32/int32 columns, and this script writes
exactly that.

    python tools/fits_to_case.py spec-0294-51986-0036.fits --out data/raw/cases
    gspec measure --cases data/raw/cases --out data/raw/measurements.csv

What is copied, and what is not
-------------------------------
Copied verbatim, in their stored precision, with no interpolation, smoothing,
normalisation or re-binning:

    flux        flux density, 1e-17 erg s^-1 cm^-2 A^-1
    loglam      log10(wavelength / Angstrom), on the archive's own grid
    ivar        1/sigma^2 per pixel
    and_mask    SDSS bitmask, 32-bit
    or_mask     SDSS bitmask, 32-bit
    wdisp       instrumental sigma, **in pixels** (see ``check_wdisp_units``)

Dropped, and recorded as dropped: the ``sky`` and ``model`` columns of COADD, and
every other extension of the archive file (SPECOBJ and SPZLINE carry the survey
pipeline's own redshift and line fluxes — they are the answers to the question
this repository asks, so they are never carried downstream).

Units are not guessed. ``BUNIT`` must be present and must be one of the known
spellings of 1e-17 erg/cm^2/s/Ang, and ``VACUUM`` must be present: without them
the wavelength scale or the flux scale would be an assumption rather than an
observed fact. If a file is incompatible, the script says so and exits 4. It does
not bend the data to fit.

Exit codes
----------
0  converted
1  usage error (bad arguments, output exists and --force not given)
2  the file is not a readable FITS file
4  the file is FITS but violates the instrument contract (see the reasons printed)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

__all__ = [
    "COADD_EXTNAME",
    "REQUIRED_COLUMNS",
    "CASE_COLUMNS",
    "IncompatibleSpectrum",
    "validate_sdss",
    "convert",
    "write_case",
]

COADD_EXTNAME = "COADD"

#: Columns the case schema needs, in the order the case writes them.
CASE_COLUMNS = ("flux", "loglam", "ivar", "and_mask", "or_mask", "wdisp")

#: Same six, as the archive names them. Kept separate so the intent stays legible.
REQUIRED_COLUMNS = CASE_COLUMNS

#: FITS TFORM for each case column: float32 for the five reals, int32 for masks.
CASE_FORMATS = {"flux": "E", "loglam": "E", "ivar": "E",
                "and_mask": "J", "or_mask": "J", "wdisp": "E"}

CASE_UNITS = {"flux": "1E-17 erg/cm^2/s/Ang",
              "loglam": "log10(Angstrom)",
              "ivar": "(1e-17 erg/s/cm2/A)^-2",
              "and_mask": "",
              "or_mask": "",
              "wdisp": "pixel"}

#: Spellings of the flux-density unit observed in SDSS-era headers. Compared
#: after normalisation, so "^2" and "**2" and "s-1" and "/s" all agree.
_ACCEPTED_BUNIT = "1e-17erg/cm2/s/ang"

#: Below this fraction of finite inverse variance the file is not usable.
MIN_FINITE_IVAR = 0.90
#: Below this fraction of positive inverse variance the file is only a warning...
WARN_IVAR_POSITIVE = 0.50
#: ...and below this one it is refused outright.
MIN_IVAR_POSITIVE = 0.20


class IncompatibleSpectrum(Exception):
    """A readable FITS file that does not meet the instrument contract."""

    def __init__(self, reasons: list[str]):
        self.reasons = list(reasons)
        super().__init__("; ".join(self.reasons))


def _normalise_bunit(text: str) -> str:
    """Fold the many spellings of the SDSS flux unit onto one comparable string.

    FITS headers in the wild write this unit as any of

        1E-17 erg/cm^2/s/Ang       1e-17 erg/s/cm2/Ang
        1E-17 erg cm-2 s-1 Ang-1   1E-17erg/cm2/s/A

    The rules below are mechanical and order-sensitive: abbreviations are
    expanded before bare ``/a`` is rewritten, and ``/a`` is only rewritten as a
    final path segment, so that expanding it can never re-fire on the ``ang`` it
    just produced.
    """
    s = str(text).strip().lower().replace(" ", "").replace("_", "")
    s = s.replace("**", "^")
    # expand the long spellings first
    s = s.replace("angstroms", "ang").replace("angstrom", "ang")
    # exponents written with a caret
    s = s.replace("^2", "2").replace("^-1", "-1")
    # inverse seconds, however written
    s = s.replace("s-1", "/s")
    # negative exponents written with a sign instead of a caret (DR7-era headers)
    s = s.replace("cm-2", "cm2").replace("a-1", "/ang").replace("ergcm2/", "erg/cm2/")
    # a bare final /a or /A means Angstrom
    if s.endswith("/a"):
        s = s[:-2] + "/ang"
    # one ordering for the density part
    s = s.replace("erg/s/cm2/", "erg/cm2/s/")
    s = s.replace("//", "/")
    return s


def find_coadd(hdul: fits.HDUList) -> tuple[int, fits.BinTableHDU]:
    """Locate the COADD table by name.

    By name, not by index: the archive file has varied its extension order across
    reductions, and index 1 is a coincidence, not a contract.
    """
    for i, hdu in enumerate(hdul):
        if getattr(hdu, "name", "").strip().upper() == COADD_EXTNAME:
            if isinstance(hdu, fits.BinTableHDU):
                return i, hdu
    names = [f"{i}:{getattr(h, 'name', '?')}" for i, h in enumerate(hdul)]
    raise IncompatibleSpectrum([f"no COADD binary table found; extensions present: {names}"])


def check_wdisp_units(lam: np.ndarray, dlam: np.ndarray, wdisp: np.ndarray,
                      valid: np.ndarray) -> dict:
    """Decide, from the data, that ``wdisp`` is in pixels rather than Angstrom.

    SDSS documents COADD.wdisp as the resolution width in pixels, but the archive
    file carries no TUNIT for it. Rather than trust the memory of a document, use
    the wavelength dependence: for a spectrograph of roughly constant resolving
    power the resolution width in *Angstrom* grows with wavelength, whereas the
    width in *pixels* falls as 1/dlam. Measured on a real run2d=26 file, wdisp
    falls from 1.25 at 4000 A to 0.73 at 8800 A while wdisp*dlam stays within
    1.1-1.4 A, and 2.355*wdisp*dlam gives R = 1470 at 4000 A and 2290 at Halpha,
    which is the published shape of the Legacy spectrograph's resolution curve.
    Multiplying by dlam is therefore the correct interpretation.

    Returns a summary, and raises if the data contradict the interpretation.
    """
    hi = valid & (lam > 7000)
    lo = valid & (lam < 4200)
    if hi.sum() < 20 or lo.sum() < 20:
        return {"status": "not testable (too few pixels at the extremes)"}
    w_hi = float(np.median(wdisp[hi]))
    w_lo = float(np.median(wdisp[lo]))
    lam_ratio = float(np.median(lam[hi]) / np.median(lam[lo]))
    # In pixels, wdisp ~ const * (1/dlam) ~ const/lambda, so w_hi/w_lo < 1.
    # In Angstrom at constant R, w_hi/w_lo would track lam_ratio (> 1).
    if w_hi / w_lo > lam_ratio * 0.6:
        raise IncompatibleSpectrum([
            f"wdisp behaves like an Angstrom width, not a pixel width "
            f"(ratio {w_hi / w_lo:.3f} against a wavelength ratio of {lam_ratio:.3f}); "
            f"the pipeline interprets wdisp as pixels"])
    sigma_a = wdisp * dlam
    n = valid & (lam > 6400) & (lam < 6700)
    r_halpha = float(lam[n].mean() / (2.355 * np.median(sigma_a[n]))) if n.sum() else float("nan")
    return {"unit": "pixel", "method": "wavelength dependence of wdisp",
            "wdisp_hi_over_lo": round(w_hi / w_lo, 4),
            "wavelength_ratio": round(lam_ratio, 4),
            "sigma_lambda_median_A": round(float(np.median(sigma_a[valid])), 4),
            "R_halpha": None if np.isnan(r_halpha) else round(r_halpha, 1)}


def validate_sdss(path_or_hdul) -> dict:
    """Check one archive file against the contract. Returns a report.

    Raises :class:`IncompatibleSpectrum` listing every violated requirement at
    once, so a user sees all of them rather than one per run.
    """
    hdul = path_or_hdul if isinstance(path_or_hdul, fits.HDUList) else fits.open(path_or_hdul)
    reasons: list[str] = []
    report: dict = {}

    try:
        prim = hdul[0].header
        i, tab = find_coadd(hdul)
    except IncompatibleSpectrum as exc:
        raise exc

    cols = [c.name for c in tab.columns]
    report["coadd_hdu_index"] = i
    report["columns_present"] = cols
    missing = [c for c in REQUIRED_COLUMNS if c not in cols]
    if missing:
        raise IncompatibleSpectrum([
            f"COADD is missing required column(s) {missing}; it has {cols}"])

    data = tab.data
    n = len(data)
    report["n_pixels"] = int(n)
    if n < 100:
        reasons.append(f"only {n} pixels; the pipeline needs a real spectrum (>= 100)")

    loglam = np.asarray(data["loglam"], dtype=float)
    flux = np.asarray(data["flux"], dtype=float)
    ivar = np.asarray(data["ivar"], dtype=float)
    wdisp = np.asarray(data["wdisp"], dtype=float)
    and_mask = np.asarray(data["and_mask"], dtype=np.int32)
    or_mask = np.asarray(data["or_mask"], dtype=np.int32)

    # --- wavelength grid -------------------------------------------------
    if not np.all(np.isfinite(loglam)):
        reasons.append(f"{int((~np.isfinite(loglam)).sum())} non-finite loglam values")
    elif not np.all(np.diff(loglam) > 0):
        bad = int((np.diff(loglam) <= 0).sum())
        reasons.append(f"loglam is not strictly increasing ({bad} non-increasing steps)")
    else:
        lam = 10.0 ** loglam
        dlam = np.gradient(lam)
        report["wavelength_range_A"] = [round(float(lam[0]), 3), round(float(lam[-1]), 3)]
        report["loglam_range"] = [round(float(loglam[0]), 7), round(float(loglam[-1]), 7)]
        report["delta_lambda_A"] = [round(float(dlam[0]), 4), round(float(dlam[-1]), 4)]

        coeff0, coeff1 = prim.get("COEFF0"), prim.get("COEFF1")
        dflag = prim.get("DC-FLAG")
        report["dc_flag"] = dflag
        if dflag is not None and int(dflag) != 1:
            reasons.append(f"DC-FLAG = {dflag}: only a log-linear dispersion (DC-FLAG=1) "
                           f"is supported; this file is not on the grid gspec assumes")
        if coeff0 is not None and coeff1 is not None:
            step = float(np.diff(loglam).mean())
            dev = float(np.abs(np.diff(loglam) - float(coeff1)).max())
            report["coeff0"], report["coeff1"] = float(coeff0), float(coeff1)
            report["loglam_step_max_deviation"] = dev
            if dev > 1e-5:
                reasons.append(f"loglam steps deviate from COEFF1={coeff1} by up to {dev:.2e}")
            offset = (float(loglam[0]) - float(coeff0)) / float(coeff1)
            report["grid_offset_pixels_from_coeff0"] = round(offset, 3)
            if offset < -1e-3 or abs(offset - round(offset)) > 1e-3:
                reasons.append(f"the first pixel is {offset:.3f} pixels from COEFF0; "
                               f"not an integer offset on the documented grid")

        # --- flux density unit -------------------------------------------
        bunit = prim.get("BUNIT")
        report["bunit_raw"] = None if bunit is None else str(bunit)
        if bunit is None:
            reasons.append("BUNIT is absent: the flux density scale would be a guess")
        elif _normalise_bunit(bunit) != _ACCEPTED_BUNIT:
            reasons.append(f"BUNIT = {bunit!r} is not a recognised spelling of "
                           f"1e-17 erg/cm^2/s/Ang")

        # --- wavelength frame -------------------------------------------
        vacuum = prim.get("VACUUM")
        report["vacuum"] = None if vacuum is None else bool(vacuum)
        if vacuum is None:
            reasons.append("VACUUM is absent: cannot tell whether the wavelengths are "
                           "vacuum or air, and the pipeline's rest wavelengths are vacuum")
        elif not bool(vacuum):
            reasons.append("VACUUM = F: the pipeline's rest wavelengths are vacuum, so an "
                           "air-wavelength spectrum is not interchangeable with it")

    # --- inverse variance ------------------------------------------------
    finite = np.isfinite(ivar)
    report["ivar_finite_fraction"] = round(float(finite.mean()), 4) if n else 0.0
    if n and finite.mean() < MIN_FINITE_IVAR:
        reasons.append(f"only {finite.mean():.1%} of ivar values are finite "
                       f"(need >= {MIN_FINITE_IVAR:.0%})")
    positive = finite & (ivar > 0)
    report["ivar_positive_fraction"] = round(float(positive.mean()), 4) if n else 0.0
    report["n_ivar_zero"] = int(n - positive.sum())
    if n and positive.mean() < MIN_IVAR_POSITIVE:
        reasons.append(f"only {positive.mean():.1%} of pixels carry positive inverse "
                       f"variance (need >= {MIN_IVAR_POSITIVE:.0%}); this is not a spectrum "
                       f"the pipeline can measure")
    report["ivar_positive_warning"] = bool(n and positive.mean() < WARN_IVAR_POSITIVE)

    # --- flags and resolution -------------------------------------------
    report["n_and_mask_nonzero"] = int((and_mask != 0).sum())
    report["n_or_mask_nonzero"] = int((or_mask != 0).sum())
    if not np.all(np.isfinite(wdisp)):
        reasons.append(f"{int((~np.isfinite(wdisp)).sum())} non-finite wdisp values")
    elif positive.any():
        zero_w = int((positive & (wdisp <= 0)).sum())
        report["n_wdisp_nonpositive_on_valid_pixels"] = zero_w
        if zero_w:
            reasons.append(f"{zero_w} pixels have positive ivar but wdisp <= 0: the "
                           f"resolution is missing exactly where the data are kept")
        if not reasons:
            try:
                report["resolution"] = check_wdisp_units(lam, dlam, wdisp, positive)
            except IncompatibleSpectrum as exc:
                reasons.extend(exc.reasons)

    report["dropped_columns"] = [c for c in cols if c not in CASE_COLUMNS]
    report["dropped_extensions"] = [f"{k}:{getattr(h, 'name', '?')}"
                                    for k, h in enumerate(hdul) if k != i]

    # --- survey metadata, recorded but never relied on -------------------
    for key in ("PLATE", "MJD", "FIBERID", "RUN2D", "RUN1D", "EXPTIME", "TELESCOP",
                "INSTRUME", "SDSSVER", "SURVEY"):
        val = prim.get(key)
        if val is not None:
            report.setdefault("header", {})[key] = val

    if reasons:
        raise IncompatibleSpectrum(reasons)
    return report


def _identifiers_from_name(path: Path) -> dict:
    """plate/MJD/fibre parsed from ``spec-{plate}-{mjd}-{fibre}.fits``, if it matches."""
    stem = path.name
    if stem.endswith(".fits"):
        stem = stem[:-5]
    parts = stem.split("-")
    if len(parts) == 4 and parts[0] == "spec":
        try:
            return {"plate": int(parts[1]), "mjd": int(parts[2]), "fiber": int(parts[3]),
                    "source": "filename"}
        except ValueError:
            pass
    return {}


def build_case_hdus(data, case_id: str, proven: str, n_pixels: int,
                    extras: dict | None = None) -> fits.HDUList:
    """Assemble the HDUList a case file consists of. Pure: no I/O."""
    prim = fits.PrimaryHDU()
    h = prim.header
    h["CASEID"] = (case_id, "opaque case identifier")
    h["BUNIT"] = ("1E-17 erg/cm^2/s/Ang", "flux density unit of COADD.flux")
    h["VACUUM"] = (True, "wavelengths are vacuum")
    h["PROVEN"] = (proven, "provenance family of this case")
    h["NWAVE"] = (int(n_pixels), "number of spectral pixels")
    for key, (value, comment) in (extras or {}).items():
        h[key] = (value, comment)

    cols = []
    for name in CASE_COLUMNS:
        arr = np.asarray(data[name])
        # float32 for real columns, int32 for masks; the storage precision of the
        # archive table, which is also what the validated baseline was built from.
        if name in ("and_mask", "or_mask"):
            arr = arr.astype(np.int32, copy=False)
        else:
            arr = arr.astype(np.float32, copy=False)
        cols.append(fits.Column(name=name, format=CASE_FORMATS[name],
                                unit=CASE_UNITS[name], array=arr))
    tab = fits.BinTableHDU.from_columns(cols, name=COADD_EXTNAME)
    return fits.HDUList([prim, tab])


def write_case(hdul: fits.HDUList, out_path: str | Path) -> Path:
    """Write atomically: a temporary file, then a rename. Never a half-written case.

    No CHECKSUM/DATASUM cards: they would embed the write time and make the same
    input produce different bytes on the next run. Integrity is recorded in the
    provenance sidecar instead, which is where a reader looks for it anyway.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".part")
    try:
        hdul.writeto(tmp, overwrite=True)
        with open(tmp, "rb") as fh:
            os.fsync(fh.fileno())
        os.replace(tmp, out_path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
    return out_path


def sha256_of(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def convert(raw_path: str | Path, out_path: str | Path, case_id: str | None = None,
            proven: str = "sdss-raw", force: bool = False,
            extra_header: dict | None = None) -> dict:
    """Convert one archive file into a case file. Returns the provenance record."""
    raw_path = Path(raw_path)
    out_path = Path(out_path)
    if out_path.exists() and not force:
        raise FileExistsError(f"{out_path} exists; pass --force to replace it")

    try:
        hdul = fits.open(raw_path)
    except OSError as exc:
        raise ValueError(f"{raw_path} is not a readable FITS file: {exc}") from exc

    with hdul:
        report = validate_sdss(hdul)
        i, tab = find_coadd(hdul)
        data = {c: np.asarray(tab.data[c]) for c in CASE_COLUMNS}
        prim = hdul[0].header
        ident = _identifiers_from_name(raw_path)
        if case_id is None:
            if ident:
                case_id = "SDSS-{plate:04d}-{mjd}-{fiber:04d}".format(**ident)
            else:
                case_id = raw_path.stem

        extras = {
            "ORIGIN": ("SDSS DR18 SAS, run2d=26", "where the archive file came from"),
            "SRCNAME": (raw_path.name, "source archive file name"),
        }
        if ident:
            extras["PLATE"] = (ident["plate"], "informational; not read by gspec")
            extras["MJD"] = (ident["mjd"], "informational; not read by gspec")
            extras["FIBERID"] = (ident["fiber"], "informational; not read by gspec")
        if prim.get("RUN2D") is not None:
            extras["RUN2D"] = (prim["RUN2D"], "spectroscopic reduction version")
        for k, v in (extra_header or {}).items():
            extras[k] = (v, "user-supplied")

        hdus = build_case_hdus(data, case_id, proven, report["n_pixels"], extras)
        write_case(hdus, out_path)

    # The file is the authority on which object it is, where it says so. Compare
    # what the name claimed with what the header records; neither is "verified"
    # beyond that, and the record says which is which.
    header_ids = {k: prim.get(k) for k in ("PLATE", "MJD", "FIBERID")
                  if prim.get(k) is not None}
    id_check = {}
    for key, name in (("PLATE", "plate"), ("MJD", "mjd"), ("FIBERID", "fiber")):
        if key in header_ids and name in ident:
            id_check[key] = bool(int(header_ids[key]) == int(ident[name]))
        elif key in header_ids:
            id_check[key] = "no name to compare"
        else:
            id_check[key] = "not in header"

    record = {
        "case_id": case_id,
        "case_path": str(out_path),
        "source": {
            "file": raw_path.name,
            "sha256": sha256_of(raw_path),
            "bytes": int(raw_path.stat().st_size),
        },
        "case": {
            "sha256": sha256_of(out_path),
            "bytes": int(out_path.stat().st_size),
        },
        "conversion": {
            "tool": "tools/fits_to_case.py",
            "transformations": "none: the six columns are copied in their stored precision",
            "columns_dropped": report["dropped_columns"],
            "extensions_dropped": report["dropped_extensions"],
            "resampling": "none",
            "renormalisation": "none",
            "smoothing": "none",
        },
        "fits": report,
        "identifiers": ident,                 # observed: parsed from the file name
        "header_identifiers": header_ids,      # observed: the file's own keywords
        "identifier_check": id_check,
    }
    return record


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Convert an SDSS spectrum file into the gspec case schema.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Exit codes: 0 converted, 1 usage error, 2 not FITS, 4 violates the "
               "instrument contract.")
    ap.add_argument("raw", help="SDSS spectrum FITS file (spec-<plate>-<mjd>-<fibre>.fits)")
    ap.add_argument("--out", required=True, help="output directory for the case file")
    ap.add_argument("--case-id", default=None,
                    help="CASEID for the case (default: SDSS-<plate>-<mjd>-<fibre>)")
    ap.add_argument("--proven", default="sdss-raw", help="PROVEN value (default: sdss-raw)")
    ap.add_argument("--force", action="store_true", help="replace an existing case file")
    ap.add_argument("--report", default=None, help="write the provenance JSON here")
    args = ap.parse_args(argv)

    raw = Path(args.raw)
    if not raw.exists():
        print(f"error: {raw} does not exist", file=sys.stderr)
        return 1
    case_name = f"{args.case_id}.fits" if args.case_id else f"{raw.stem}.fits"
    out = Path(args.out) / case_name

    try:
        record = convert(raw, out, case_id=args.case_id, proven=args.proven,
                         force=args.force)
    except FileExistsError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except IncompatibleSpectrum as exc:
        print(f"error: {raw.name} is not compatible with the pipeline:", file=sys.stderr)
        for r in exc.reasons:
            print(f"  - {r}", file=sys.stderr)
        return 4

    if args.report:
        Path(args.report).parent.mkdir(parents=True, exist_ok=True)
        Path(args.report).write_text(json.dumps(record, indent=2) + "\n")

    f = record["fits"]
    print(f"wrote {record['case_path']}")
    print(f"  CASEID {record['case_id']}  {f['n_pixels']} pixels  "
          f"{f['wavelength_range_A'][0]:.1f}-{f['wavelength_range_A'][1]:.1f} A")
    print(f"  ivar > 0 on {f['ivar_positive_fraction']:.1%} of pixels  "
          f"R(Halpha) ~ {f['resolution'].get('R_halpha')}")
    print(f"  case sha256 {record['case']['sha256'][:16]}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
