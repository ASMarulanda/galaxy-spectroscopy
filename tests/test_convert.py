"""The converter: does it produce a case gspec can read, and does it refuse bad input?

Every refusal here is a requirement from the instrument contract, not a preference:
a wavelength grid that is not monotonic, a flux scale that is not stated, or an
inverse-variance column with no information would each silently corrupt a
measurement if they were allowed through.

The positive test is stronger than "it runs": converting a rebuilt archive file
must reproduce ``R06`` bit for bit, because DATA.md §3 records that the case files
keep the archive's COADD columns unchanged.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from astropy.io import fits

import fits_to_case
import gspec
from conftest import CASES


# ------------------------------------------------------------------ the good path
def test_converted_case_is_bit_identical_to_the_shipped_one(make_raw, tmp_path):
    raw = make_raw()
    out = tmp_path / "cases" / "x.fits"
    record = fits_to_case.convert(raw, out)

    new = fits.open(out)[1].data
    ref = fits.open(CASES / "R06.fits")[1].data
    for col in fits_to_case.CASE_COLUMNS:
        a, b = np.asarray(new[col]), np.asarray(ref[col])
        assert np.array_equal(a, b), f"{col} differs from the shipped case"
        assert a.dtype == b.dtype, f"{col} dtype differs: {a.dtype} vs {b.dtype}"
    assert len(new) == len(ref)
    assert record["conversion"]["transformations"].startswith("none")


def test_case_schema_is_the_one_read_case_expects(make_raw, tmp_path):
    raw = make_raw()
    out = tmp_path / "case.fits"
    fits_to_case.convert(raw, out, case_id="TEST-1")

    with fits.open(out) as hdul:
        prim = hdul[0].header
        assert prim["CASEID"] == "TEST-1"
        assert prim["BUNIT"] == "1E-17 erg/cm^2/s/Ang"
        assert prim["PROVEN"] == "sdss-raw"
        assert prim["VACUUM"] is True
        assert prim["NWAVE"] == len(hdul[1].data)
        tab = hdul[1]
        assert tab.name == "COADD"
        assert [c.name for c in tab.columns] == list(fits_to_case.CASE_COLUMNS)
        assert [c.format for c in tab.columns] == ["E", "E", "E", "J", "J", "E"]
        assert tab.header["NAXIS1"] == 24          # six 4-byte columns, as in the baseline
        # the "answers" extension is gone
        assert len(hdul) == 2

    case = gspec.read_case(out)                     # the real reader, not a copy of it
    assert case.case_id == "TEST-1"
    assert case.wave.shape == (len(case),)
    assert np.all(np.diff(case.wave) > 0)


def test_conversion_does_not_resample_or_change_length(make_raw, tmp_path):
    raw = make_raw()
    out = tmp_path / "case.fits"
    fits_to_case.convert(raw, out)
    n_in = len(fits.open(raw)[1].data)
    n_out = len(fits.open(out)[1].data)
    assert n_in == n_out


def test_conversion_is_byte_deterministic(make_raw, tmp_path):
    raw = make_raw()
    a, b = tmp_path / "a.fits", tmp_path / "b.fits"
    fits_to_case.convert(raw, a)
    fits_to_case.convert(raw, b)
    assert a.read_bytes() == b.read_bytes()


# --------------------------------------------------------------------- refusals
def test_missing_column_names_the_ones_present(make_raw, tmp_path):
    raw = make_raw(drop=("wdisp",))
    with pytest.raises(fits_to_case.IncompatibleSpectrum) as exc:
        fits_to_case.convert(raw, tmp_path / "x.fits")
    msg = str(exc.value)
    assert "wdisp" in msg and "missing" in msg
    assert "flux" in msg            # the actual column list is reported, not a generic error


def test_non_monotonic_loglam_is_refused(make_raw, tmp_path):
    lam = np.full(200, 3.8, dtype=np.float32)
    lam[100:] = 3.79                 # a step backwards
    raw = make_raw(loglam=lam, n_rows=200)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="not strictly increasing"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_absent_bunit_is_refused(make_raw, tmp_path):
    raw = make_raw(bunit="")
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="BUNIT"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


@pytest.mark.parametrize("bunit", ["counts", "1e-16 erg/cm^2/s/Ang", "ADU",
                                   "erg/cm^2/s/Ang"])
def test_foreign_flux_units_are_refused(make_raw, tmp_path, bunit):
    raw = make_raw(bunit=bunit)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="BUNIT"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


@pytest.mark.parametrize("bunit", ["1E-17 erg/cm^2/s/Ang", "1e-17 erg/s/cm2/Ang",
                                   "1E-17 erg/cm2/s/angstrom", "1e-17erg/cm2/s/A"])
def test_every_spelling_of_the_expected_unit_is_accepted(make_raw, tmp_path, bunit):
    raw = make_raw(bunit=bunit)
    rec = fits_to_case.convert(raw, tmp_path / "x.fits")
    assert rec["fits"]["bunit_raw"] == bunit


def test_absent_vacuum_keyword_is_refused(make_raw, tmp_path):
    raw = make_raw(vacuum=None)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="VACUUM is absent"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_air_wavelengths_are_refused(make_raw, tmp_path):
    raw = make_raw(vacuum=False)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="VACUUM = F"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_all_zero_ivar_is_refused(make_raw, tmp_path):
    n = 300
    raw = make_raw(ivar=np.zeros(n), n_rows=n)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="positive inverse"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_linear_dispersion_is_refused(make_raw, tmp_path):
    raw = make_raw(dflag=0)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="DC-FLAG"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_wrong_grid_step_is_refused(make_raw, tmp_path):
    raw = make_raw(coeff1=2e-4)          # header claims a step the pixels do not have
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="COEFF1"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_non_integer_grid_offset_is_refused(make_raw, tmp_path):
    raw = make_raw(coeff0=3.99999)       # 3828 pixels away, and not a whole pixel
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="pixels from COEFF0"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_wdisp_measured_in_angstrom_is_refused(make_raw, tmp_path):
    """If wdisp grew with wavelength it would be an Angstrom width, not a pixel one."""
    step = 1e-3
    n = 353                                   # 4000 A to 9000 A on a log grid
    loglam = (np.log10(4000.0) + np.arange(n) * step).astype(np.float32)
    lam = 10.0 ** loglam
    # an Angstrom-width resolution: grows with wavelength, unlike a pixel width
    raw = make_raw(loglam=loglam, n_rows=n, coeff0=float(np.log10(4000.0)),
                   coeff1=step, wdisp=(lam / 3000.0).astype(np.float32))
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="Angstrom width"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_zero_wdisp_on_valid_pixels_is_refused(make_raw, tmp_path):
    n = 300
    w = np.ones(n, dtype=np.float32)
    w[:50] = 0.0
    raw = make_raw(wdisp=w, n_rows=n)
    with pytest.raises(fits_to_case.IncompatibleSpectrum, match="wdisp <= 0"):
        fits_to_case.convert(raw, tmp_path / "x.fits")


def test_refuses_to_overwrite_without_force(make_raw, tmp_path):
    raw = make_raw()
    out = tmp_path / "case.fits"
    fits_to_case.convert(raw, out)
    before = out.read_bytes()
    with pytest.raises(FileExistsError):
        fits_to_case.convert(raw, out)
    assert out.read_bytes() == before
    fits_to_case.convert(raw, out, force=True)          # explicit replacement is fine


def test_a_non_fits_file_is_an_error_not_a_crash(tmp_path):
    bogus = tmp_path / "not.fits"
    bogus.write_text("<html>404</html>")
    with pytest.raises(ValueError, match="not a readable FITS"):
        fits_to_case.convert(bogus, tmp_path / "x.fits")


# ------------------------------------------------------------------- provenance
def test_provenance_record_separates_observed_from_supplied(make_raw, tmp_path):
    raw = make_raw()
    rec = fits_to_case.convert(raw, tmp_path / "case.fits", case_id="USER-CHOSEN")

    assert rec["case_id"] == "USER-CHOSEN"
    assert rec["source"]["sha256"] == fits_to_case.sha256_of(raw)
    assert rec["case"]["sha256"] == fits_to_case.sha256_of(tmp_path / "case.fits")
    assert rec["conversion"]["resampling"] == "none"
    assert rec["conversion"]["smoothing"] == "none"
    assert rec["conversion"]["renormalisation"] == "none"

    fits_report = rec["fits"]
    assert fits_report["n_pixels"] > 100
    assert 0.0 < fits_report["ivar_positive_fraction"] <= 1.0
    assert fits_report["resolution"]["unit"] == "pixel"
    assert 1500 < fits_report["resolution"]["R_halpha"] < 3500
    # the identifier came from the file name, and says so
    assert rec["identifiers"]["source"] == "filename"
    assert rec["identifiers"]["plate"] == 294
    # SPECOBJ was present and is recorded as dropped
    assert any("SPECOBJ" in e for e in rec["conversion"]["extensions_dropped"])


def test_provenance_json_round_trips(make_raw, tmp_path):
    raw = make_raw()
    rec = fits_to_case.convert(raw, tmp_path / "case.fits")
    text = json.dumps(rec)
    assert json.loads(text)["case"]["sha256"] == rec["case"]["sha256"]


# ------------------------------------------------------------------------- CLI
def test_cli_reports_incompatibility_with_exit_4(make_raw, tmp_path, capsys):
    raw = make_raw(bunit="counts")
    code = fits_to_case.main([str(raw), "--out", str(tmp_path / "out")])
    captured = capsys.readouterr()
    assert code == 4
    assert "not compatible" in captured.err
    assert "BUNIT" in captured.err


def test_cli_success_writes_case_and_report(make_raw, tmp_path, capsys):
    raw = make_raw()
    report = tmp_path / "prov" / "x.json"
    code = fits_to_case.main([str(raw), "--out", str(tmp_path / "cases"),
                              "--report", str(report)])
    assert code == 0
    assert report.exists()
    assert json.loads(report.read_text())["fits"]["n_pixels"] > 100
    assert "R(Halpha)" in capsys.readouterr().out


def test_cli_missing_input_is_exit_1(tmp_path, capsys):
    code = fits_to_case.main([str(tmp_path / "nope.fits"), "--out", str(tmp_path)])
    assert code == 1
    assert "does not exist" in capsys.readouterr().err


def test_identifier_cross_check_is_recorded(make_raw, tmp_path):
    """The name and the header are both observed; the record says which agreed."""
    raw = make_raw()                       # header carries MJD 51986 and FIBERID 36
    rec = fits_to_case.convert(raw, tmp_path / "case.fits")
    assert rec["identifiers"]["mjd"] == 51986 and rec["identifiers"]["source"] == "filename"
    assert rec["header_identifiers"] == {"MJD": 51986, "FIBERID": 36}
    assert rec["identifier_check"]["MJD"] is True
    assert rec["identifier_check"]["FIBERID"] is True
    assert rec["identifier_check"]["PLATE"] == "not in header"


def test_identifier_mismatch_is_reported_not_hidden(make_raw, tmp_path):
    """A file whose header disagrees with its name is flagged, not silently accepted."""
    raw = make_raw(name="spec-0294-51986-0042.fits")   # name says fibre 42
    rec = fits_to_case.convert(raw, tmp_path / "case.fits")
    assert rec["identifiers"]["fiber"] == 42
    assert rec["identifier_check"]["FIBERID"] is False
