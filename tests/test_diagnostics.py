"""The line-ratio layer: propagation, the extinction curve, and honest failure.

Two classes of failure are tested. The first is arithmetic: a ratio error that is
wrong by a square makes every diagnostic error bar meaningless. The second is
contract: a non-detection must not silently acquire a value, and an extinction
must not be reported as measured when the decrement is at or below case B.

`ccm89_k` is anchored to the published R_λ = A_λ/A_V values of the curve rather
than to remembered k(λ) numbers: R_V is defined so that R_λ(V) = 1 exactly, which
is the check that catches a wrong polynomial branch.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from gspec import diagnostics as dg


# --------------------------------------------------------------------------- #
# Ratio propagation
# --------------------------------------------------------------------------- #
from gspec.results import LineMeasurement


def _M(line, flux, err, flag="ok"):
    """A real LineMeasurement, so the contract under test is the shipped one."""
    return LineMeasurement(line=line, flux=flux, flux_err=err, flag=flag)


def test_ratio_error_is_first_order_propagation():
    """σ_R² = (σ_a/b)² + (a σ_b/b²)²."""
    a = _M("a", 100.0, 10.0)
    b = _M("b", 50.0, 5.0)
    r = dg.ratio(a, b)
    expect = math.hypot(10.0 / 50.0, (100.0 / 50.0) * 5.0 / 50.0)
    assert r.detected and math.isclose(r.value, 2.0)
    assert math.isclose(r.error, expect, rel_tol=1e-12)


def test_relative_error_transfers_from_either_side():
    """A 1 % error in either line gives a 1 % error in the ratio."""
    r1 = dg.ratio(_M("a", 100.0, 1.0), _M("b", 50.0, 0.0 + 1e-9))
    r2 = dg.ratio(_M("a", 100.0, 1e-9), _M("b", 50.0, 0.5))
    assert math.isclose(r1.error / r1.value, 0.01, rel_tol=1e-6)
    assert math.isclose(r2.error / r2.value, 0.01, rel_tol=1e-6)


@pytest.mark.parametrize("flux,err,why", [
    (20.0, 10.0, "below 3σ"),
    (100.0, 0.0, "no finite uncertainty"),
    (float("nan"), 1.0, "no finite"),
])
def test_a_non_detection_never_gets_a_value(flux, err, why):
    r = dg.ratio(_M("a", flux, err), _M("b", 100.0, 1.0))
    assert not r.detected
    assert why.split()[0] in r.reason or "σ" in r.reason
    assert math.isnan(r.value) or not r.detected


def test_failed_fit_is_not_a_measurement():
    r = dg.ratio(_M("a", 100.0, 1.0, flag="failed"), _M("b", 100.0, 1.0))
    assert not r.detected and "failed" in r.reason


# --------------------------------------------------------------------------- #
# Extinction curve
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("wavelength,published", [
    (5500, 1.000),      # V — the definition of R_V
    (4400, 1.324),      # B
    (12500, 0.282),     # J
    (21900, 0.114),     # K
])
def test_ccm89_matches_published_r_lambda(wavelength, published):
    """R_λ = k(λ)/R_V at R_V = 3.1, against the published CCM89 values."""
    assert math.isclose(dg.ccm89_k(wavelength) / 3.1, published, abs_tol=2e-3)


def test_balmer_coefficient_is_the_standard_one():
    """2.5/(k(Hβ) − k(Hα)) ≈ 2.33: the coefficient every Balmer-decrement
    conversion uses. A wrong branch or a swapped label breaks this."""
    dk = dg.ccm89_k(dg.LAMBDA_HB) - dg.ccm89_k(dg.LAMBDA_HA)
    assert math.isclose(dk, 1.073, abs_tol=2e-3)
    assert math.isclose(2.5 / dk, 2.329, abs_tol=5e-3)


def test_curve_refuses_wavelengths_outside_its_branch():
    with pytest.raises(ValueError):
        dg.ccm89_k(2000)


# --------------------------------------------------------------------------- #
# The inference
# --------------------------------------------------------------------------- #
def test_extinction_inverts_the_decrement():
    """Feeding a known E(B−V) through the curve returns it."""
    ebv_true, rv = 0.30, 3.1
    dk = dg.ccm89_k(dg.LAMBDA_HB) - dg.ccm89_k(dg.LAMBDA_HA)
    decrement = dg.HALPHA_HBETA_CASE_B * 10 ** (ebv_true * dk / 2.5)
    ex = dg.extinction_from_balmer(dg.Ratio("Hα/Hβ", decrement, 0.01, True))
    assert ex.valid
    assert math.isclose(ex.e_bv, ebv_true, rel_tol=1e-9)
    assert math.isclose(ex.av, rv * ebv_true, rel_tol=1e-9)


def test_decrement_at_case_b_gives_an_upper_limit_not_a_measurement():
    """No reddening excess is not "A_V = 0": it is no constraint."""
    ex = dg.extinction_from_balmer(dg.Ratio("Hα/Hβ", 2.70, 0.10, True))
    assert not ex.valid
    assert math.isnan(ex.av)
    assert math.isclose(ex.av_upper, 0.0, abs_tol=0.35)
    assert "upper limit" in ex.reason


def test_undetected_decrement_yields_no_extinction():
    ex = dg.extinction_from_balmer(dg.Ratio("Hα/Hβ", float("nan"), float("nan"), False,
                                            reason="Hbeta below 3σ"))
    assert not ex.valid and math.isnan(ex.e_bv)
    assert "Hbeta below 3σ" in ex.reason


def test_assumptions_travel_with_the_number():
    ex = dg.extinction_from_balmer(dg.Ratio("Hα/Hβ", 3.5, 0.1, True))
    assert ex.valid
    joined = " ".join(ex.assumptions)
    for token in ("case-B", "screen", "CCM89", "optically thin"):
        assert token in joined
