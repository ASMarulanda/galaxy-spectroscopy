"""The whole workflow, with no network: bytes in, the validated numbers out.

This is the test that makes V0.3 more than plumbing. A spectrum is served by a
fake archive, downloaded through the real fetch path, converted by the real
converter, and measured by the real pipeline. The result must equal the rows the
validated baseline already holds for that object — every field, to the last
printed digit.

The fixture is rebuilt from ``R06.fits``, which DATA.md §3 records as carrying the
archive's COADD columns unchanged, so this is a statement about the actual data and
not about a toy spectrum.
"""

from __future__ import annotations

import csv
import hashlib
import json

import pytest
from astropy.io import fits

import download_sdss as D
import fits_to_case
import gspec
from conftest import CASES, ROOT

BASELINE = ROOT / "validation" / "results" / "golden_reference.csv"
FIELDS = ("redshift", "redshift_err", "flux", "flux_err", "flag")


class ServingOpener:
    """An archive that serves one real file, through UrllibOpener's interface."""

    def __init__(self, payload: bytes):
        self.payload = payload
        self.head_calls = 0
        self.get_calls = 0

    def head(self, url):
        self.head_calls += 1
        return D.HeadInfo(status=200, content_length=len(self.payload),
                          etag='"real"', last_modified="Wed, 02 May 2012 20:29:34 GMT",
                          accept_ranges="bytes")

    def get(self, url, start=0):
        self.get_calls += 1
        import io

        class R:
            status = 206 if start else 200

            def __init__(self, buf):
                self._b = io.BytesIO(buf)

            def read(self, n=-1):
                return self._b.read(n)

            def close(self):
                pass

        return R(self.payload[start:])


@pytest.fixture(scope="module")
def baseline_rows():
    with open(BASELINE, newline="") as fh:
        return {r["line"]: r for r in csv.DictReader(fh) if r["case_id"] == "R06"}


def test_end_to_end_reproduces_the_validated_baseline(make_raw, tmp_path, baseline_rows):
    raw = make_raw()
    payload = raw.read_bytes()
    opener = ServingOpener(payload)

    out = tmp_path / "run"
    result = D.fetch_one(294, 51986, 36, out / "raw", opener=opener)
    assert result.sha256 == hashlib.sha256(payload).hexdigest()

    case_path = out / "cases" / "spec-0294-51986-0036.fits"
    record = fits_to_case.convert(result.path, case_path)
    assert record["case"]["sha256"] == fits_to_case.sha256_of(case_path)

    measured = gspec.measure(gspec.read_case(case_path))
    assert measured.case_id == "SDSS-0294-51986-0036"

    # Portable: the numbers land on the baseline within tolerances that hold on any
    # machine (see tests/test_regression.py for where the numbers come from). The
    # exact-digit comparison is a characterisation test below.
    for line in measured.lines:
        row = measured.row(line)
        expected = baseline_rows[row["line"]]
        assert row["flag"] == expected["flag"], f"{line}: flag moved"
        for field, rtol in (("redshift", 1e-6), ("flux", 1e-3)):
            a, b = float(row[field]), float(expected[field])
            assert abs(a - b) <= rtol * max(abs(a), abs(b)), \
                f"{line}.{field}: {a!r} vs baseline {b!r}"
        if not measured.at_bounds:
            a, b = float(row["flux_err"]), float(expected["flux_err"])
            assert abs(a - b) <= 1e-4 * max(abs(a), abs(b)), \
                f"{line}.flux_err: {a!r} vs baseline {b!r}"

    # and the same thing through the CSV writer, which is what a reader compares
    from gspec.cli import measure_directory
    csv_path = tmp_path / "measurements.csv"
    measure_directory(out / "cases", csv_path, verbose=False)
    with open(csv_path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 4
    for r in rows:
        expected = baseline_rows[r["line"]]
        assert r["flag"] == expected["flag"]
        for field, rtol in (("redshift", 1e-6), ("flux", 1e-3)):
            a, b = float(r[field]), float(expected[field])
            assert abs(a - b) <= rtol * max(abs(a), abs(b))


def test_cli_runs_the_whole_workflow_and_writes_provenance(make_raw, tmp_path,
                                                          monkeypatch, baseline_rows):
    raw = make_raw()
    payload = raw.read_bytes()
    monkeypatch.setattr(D, "UrllibOpener", lambda **kw: ServingOpener(payload))

    code = D.main(["--plate", "294", "--mjd", "51986", "--fiber", "36",
                   "--out", str(tmp_path / "work"), "--measure"])
    assert code == 0

    # the files a reader would look for
    assert (tmp_path / "work" / "spec-0294-51986-0036.fits").exists()
    assert (tmp_path / "work" / "cases" / "spec-0294-51986-0036.fits").exists()
    assert (tmp_path / "work" / "measurements.csv").exists()

    manifest = json.loads((tmp_path / "work" / "provenance" / "fetch_manifest.json").read_text())
    assert manifest["ok"] == 1 and manifest["failed"] == 0
    obj = manifest["objects"][0]
    assert obj["status"] == "fetched"
    assert obj["sha256"] == hashlib.sha256(payload).hexdigest()
    assert obj["http"]["content_length"] == len(payload)
    assert obj["case"]["sha256"] == fits_to_case.sha256_of(
        tmp_path / "work" / "cases" / "spec-0294-51986-0036.fits")

    # the sidecars exist and agree with the files
    prov = json.loads((tmp_path / "work" / "provenance" /
                       "spec-0294-51986-0036.provenance.json").read_text())
    conv = json.loads((tmp_path / "work" / "provenance" /
                       "spec-0294-51986-0036.conversion.json").read_text())
    assert prov["sha256"] == conv["source"]["sha256"] == obj["sha256"]
    assert conv["case"]["sha256"] == obj["case"]["sha256"]

    # and the measured numbers are the baseline's, within portable tolerances
    with open(tmp_path / "work" / "measurements.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            expected = baseline_rows[r["line"]]
            a, b = float(r["flux"]), float(expected["flux"])
            assert abs(a - b) <= 1e-3 * max(abs(a), abs(b)), \
                f"{r['line']}: {a!r} vs baseline {b!r}"
            assert r["flag"] == expected["flag"]


def test_second_run_is_idempotent_offline(make_raw, tmp_path, monkeypatch):
    raw = make_raw()
    payload = raw.read_bytes()
    monkeypatch.setattr(D, "UrllibOpener", lambda **kw: ServingOpener(payload))
    args = ["--plate", "294", "--mjd", "51986", "--fiber", "36", "--out", str(tmp_path / "work")]
    assert D.main(args) == 0
    case_bytes = (tmp_path / "work" / "cases" / "spec-0294-51986-0036.fits").read_bytes()
    assert D.main(args) == 0
    assert (tmp_path / "work" / "cases" / "spec-0294-51986-0036.fits").read_bytes() == case_bytes


def test_a_corrupted_download_is_never_measured(make_raw, tmp_path, monkeypatch, capsys):
    raw = make_raw()
    payload = raw.read_bytes()
    monkeypatch.setattr(D, "UrllibOpener", lambda **kw: ServingOpener(payload))
    code = D.main(["--plate", "294", "--mjd", "51986", "--fiber", "36",
                   "--out", str(tmp_path / "work"), "--sha256", "0" * 64])
    assert code == 2                              # verification failed
    assert not (tmp_path / "work" / "cases" / "spec-0294-51986-0036.fits").exists()
    assert not (tmp_path / "work" / "spec-0294-51986-0036.fits").exists()
    assert "SHA-256 mismatch" in capsys.readouterr().err


def test_raw_and_case_hold_the_same_six_columns(make_raw, tmp_path):
    """The converter changes no number: the COADD columns survive the trip."""
    raw = make_raw()
    case = tmp_path / "case.fits"
    fits_to_case.convert(raw, case)
    a = fits.open(raw)[1].data
    b = fits.open(case)[1].data
    for col in fits_to_case.CASE_COLUMNS:
        assert (a[col] == b[col]).all()
    # the archive's sky/model columns are not carried into the case
    assert {c.name for c in fits.open(raw)[1].columns} - set(fits_to_case.CASE_COLUMNS) \
        == {"sky", "model"}
    assert [c.name for c in fits.open(case)[1].columns] == list(fits_to_case.CASE_COLUMNS)


@pytest.mark.characterisation
def test_end_to_end_reproduces_the_baseline_digit_for_digit(make_raw, tmp_path,
                                                           baseline_rows):
    """Reference environment only: every printed digit equals the baseline row.

    The portable test above allows 1e-3 on flux and 1e-6 on redshift; this one
    allows nothing. It holds on the stack the baseline was produced on and not on
    the wheel-based environment CI uses, which drifts in the 8th-10th digit and, on
    bound-pinned sigma rows, onto a different branch of the pseudo-inverse
    truncation rule (limitations L17).
    """
    raw = make_raw()
    opener = ServingOpener(raw.read_bytes())
    out = tmp_path / "run"
    result = D.fetch_one(294, 51986, 36, out / "raw", opener=opener)
    case_path = out / "cases" / "spec-0294-51986-0036.fits"
    fits_to_case.convert(result.path, case_path)

    measured = gspec.measure(gspec.read_case(case_path))
    for line in measured.lines:
        row = measured.row(line)
        expected = baseline_rows[row["line"]]
        for field in FIELDS:
            v = row[field]
            got = f"{v:.10g}" if isinstance(v, float) else str(v)
            assert got == expected[field], \
                f"{line}.{field}: measured {v!r} vs baseline {expected[field]!r}"
