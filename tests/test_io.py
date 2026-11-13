"""The FITS contract and the mask policy.

These tests are about the *interface with the data*, not about the science. They
exist because every downstream number depends on them and because the two most
common silent failures in this kind of code live here: a reader that quietly
re-casts or re-bins, and a mask policy that treats "no information" as "zero
variance".
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

import gspec
from gspec import constants as K


def test_s01_reads_with_the_documented_geometry(cases):
    case = cases["S01"]
    assert case.case_id == "S01"
    assert len(case) == 3849
    assert case.wave[0] == pytest.approx(3800.0, abs=1e-3)
    assert case.wave[-1] == pytest.approx(9200.0, abs=1e-2)
    # constant in log-lambda, which is what makes the local dlambda vary
    dlog = np.median(np.diff(case.loglam))
    assert dlog == pytest.approx(9.99e-5, rel=1e-3)


def test_units_are_only_in_the_primary_header(data_dir):
    """A reader that trusts the table's own header gets numbers with no unit."""
    with fits.open(data_dir / "cases" / "S01.fits") as hdul:
        assert hdul[0].header["BUNIT"] == "1E-17 erg/cm^2/s/Ang"
        assert hdul[0].header["VACUUM"] is True
        assert "BUNIT" not in hdul[1].header


def test_columns_are_stored_as_float32(data_dir):
    """The precision limitation is a measured fact, not a footnote.

    Everything is computed in float64, but the *input* carries ~1e-7 relative
    quantisation, so no reproducibility claim tighter than that can be tested
    against these files (DATA.md §5).
    """
    with fits.open(data_dir / "cases" / "S01.fits") as hdul:
        formats = {c.name: c.format for c in hdul[1].columns}
    assert formats == {
        "flux": "E", "loglam": "E", "ivar": "E",
        "and_mask": "J", "or_mask": "J", "wdisp": "E",
    }
    # and the float32 nature is visible in the data, not just the header
    with fits.open(data_dir / "cases" / "S01.fits") as hdul:
        assert hdul[1].data["flux"].dtype.itemsize == 4


def test_wave_is_exactly_ten_to_the_loglam(cases):
    case = cases["S01"]
    np.testing.assert_array_equal(case.wave, 10.0 ** case.loglam)


def test_dlam_is_local_not_global(cases):
    """A single global delta-lambda would bias every width by ~10 %."""
    case = cases["S02"]
    ratio = case.dlam[-1] / case.dlam[0]
    assert ratio == pytest.approx(9200.0 / 3800.0, rel=2e-2)
    assert ratio > 2.0


def test_fatal_mask_bits_are_the_documented_three():
    assert K.FATAL_MASK_BITS == K.MASK_FULLREJECT | K.MASK_NOSKY | K.MASK_NODATA
    # the caution bits must NOT be fatal, or usable pixels are discarded
    for bit in (K.MASK_NEARBADPIXEL, K.MASK_LOWFLAT, K.MASK_BRIGHTSKY):
        assert not (K.FATAL_MASK_BITS & bit)


def test_clean_rejects_zero_ivar_and_fatal_bits(data_dir):
    # Deliberately a fresh read: the session `cases` fixture is shared, and a test
    # that mutates it would corrupt every other test's data.
    case = gspec.read_case(data_dir / "cases" / "S03.fits")

    case.ivar[10] = 0.0
    case.and_mask[20] |= K.MASK_FULLREJECT
    case.and_mask[30] |= K.MASK_BRIGHTSKY      # informational, must stay usable
    mask = gspec.clean(case)
    assert not mask[10], "ivar == 0 must not constrain the fit"
    assert not mask[20], "FULLREJECT is fatal"
    assert mask[30], "BRIGHTSKY is informational and must remain usable"


def test_masked_pixels_have_zero_inverse_variance(data_dir):
    """The synthetic generator sets ivar = 0 on its bad pixels; verify the coupling."""
    with fits.open(data_dir / "cases" / "S08.fits") as hdul:
        ivar = np.asarray(hdul[1].data["ivar"], float)
        and_mask = np.asarray(hdul[1].data["and_mask"], np.int32)
    bad = (and_mask & K.MASK_FULLREJECT) != 0
    assert bad.sum() > 0, "S08 is the case designed with masked pixels"
    assert np.all(ivar[bad] == 0.0)
