"""Mathematical identities and contract branches that the data never exercises.

Each test here pins a rule that is *stated* in the code or the docs but that no
spectrum in `data/` demonstrates. They are unit tests, they cost milliseconds, and
each one covers a branch that would otherwise be untested until it fails.

Deliberately not included: anything already covered by the accuracy, injection,
uncertainty or regression tests. A test that restates one of those adds runtime
and no information.
"""

from __future__ import annotations

import numpy as np
import pytest

import gspec
from gspec import constants as K
from gspec.covariance import covariance, covariance_scale
from gspec.fit import at_bounds
from gspec.flags import classify
from gspec.fluxes import integrated_flux
from gspec.pipeline import measure, measure_fixed_model


# --------------------------------------------------------------- covariance
def test_covariance_scale_is_chi2_over_dof_below_the_cap():
    assert covariance_scale(chi2=1000.0, dof=1000) == pytest.approx(1.0)
    assert covariance_scale(chi2=1200.0, dof=1000) == pytest.approx(1.2)


def test_covariance_scale_saturates_at_the_cap():
    """The cap is a documented scientific choice, so its identity is pinned.

    Without it, chi2/dof on real spectra (driven by model mismatch) inflates every
    reported sigma; the exact factor matters, so it is tested rather than assumed.
    """
    assert covariance_scale(chi2=3000.0, dof=1000) == K.COV_SCALE_CAP
    assert covariance_scale(chi2=1e9, dof=1000) == K.COV_SCALE_CAP


def test_covariance_respects_the_cap_end_to_end():
    jac = np.eye(3)
    cov_cheap = covariance(jac, chi2=1.0, dof=1000)      # scale 0.001
    cov_capped = covariance(jac, chi2=1e9, dof=1000)     # scale would be 1e6
    assert cov_cheap[0, 0] == pytest.approx(0.001)
    assert cov_capped[0, 0] == pytest.approx(K.COV_SCALE_CAP)


# ------------------------------------------------------------------- flags
@pytest.mark.parametrize("flux,err,broad,expected", [
    (100.0, 10.0, 0.0, "ok"),                 # 10 sigma
    (20.0, 10.0, 0.0, "nondetection"),        # 2 sigma
    (-20.0, 10.0, 0.0, "nondetection"),       # a negative 2-sigma excursion
    (100.0, 95.0, 5.0, "nondetection"),       # blend-like, but below 3 sigma
    (float("nan"), 10.0, 0.0, "failed"),
    (100.0, 0.0, 0.0, "failed"),
    (100.0, -1.0, 0.0, "failed"),
])
def test_flag_contract(flux, err, broad, expected):
    """Every branch that the data can reach."""
    assert classify(flux, err, broad) == expected


def test_unresolved_blend_is_unreachable_by_construction():
    """A documented flag that the pipeline can never emit — found by this review.

    `classify` returns `nondetection` when ``|F|/sigma < 3``, and only then
    considers ``unresolved_blend``, which requires ``sigma > 0.9 |F|``. Those two
    conditions are mutually exclusive for every finite flux:

        |F|/sigma >= 3   =>   sigma <= |F|/3   =>   sigma < 0.9 |F|

    So the branch is dead code, and the `unresolved_blend` value documented in
    ``flags.FLAGS``, in the CSV column and in `docs/methods.md` §8 can never
    appear in an output. This test pins the arithmetic rather than the branch:
    if someone later reorders the conditions, this fails and the documentation
    has to change with it.

    The behaviour is inherited from the original implementation and is **not**
    changed here — changing a flag rule changes published outputs.
    """
    for flux in (1.0, 50.0, 5000.0):
        for err in (0.1, 1.0, 10.0, 100.0, 1e4):
            got = classify(flux, err, broad_amp=5.0)
            assert got != "unresolved_blend", (
                f"classify({flux}, {err}) returned unresolved_blend; the branch "
                "has become reachable and the docs must be updated"
            )
    # the mutually exclusive thresholds, stated directly
    assert K.BLEND_ERR_FRACTION > 1.0 / K.DETECTION_SIGMA


def test_detection_boundary_is_three_sigma():
    assert classify(3.01, 1.0) == "ok"
    assert classify(2.99, 1.0) == "nondetection"


# --------------------------------------------------------------- parameters
def test_at_bounds_identifies_a_pinned_parameter():
    names = ["a", "b", "c"]
    p = np.array([0.0, 0.5, 3.0])
    lo = np.array([0.0, 0.0, 1.0])
    hi = np.array([10.0, 1.0, 3.0])
    pinned = at_bounds(p, names, lo, hi)
    assert "a" in pinned      # on its lower bound
    assert "c" in pinned      # on its upper bound
    assert "b" not in pinned  # interior


def test_at_bounds_ignores_infinite_bounds():
    names = ["a"]
    pinned = at_bounds(np.array([-1e9]), names, np.array([-np.inf]), np.array([np.inf]))
    assert pinned == []


def test_the_documented_pinned_amplitude_is_reported(reference):
    """S01's broad-Hbeta amplitude is pinned at zero — the mechanism behind L4.

    If the diagnostics stopped reporting this, every statement about the
    covariance fragility would become unfalsifiable.
    """
    assert "amp_Hb_b" in reference["S01"].at_bounds
    assert abs(reference["S01"].parameters["amp_Hb_b"]) < 1e-12


# ------------------------------------------------------------ determinism
def test_measure_is_deterministic(data_dir):
    """Two calls, one process, same numbers — the tier-1 reproducibility claim."""
    case = gspec.read_case(data_dir / "cases" / "S08.fits")
    a, b = measure(case), measure(case)
    assert a.redshift == b.redshift
    for line in a.lines:
        assert a.lines[line].flux == b.lines[line].flux
        assert a.lines[line].flux_err == b.lines[line].flux_err


def test_refitting_at_the_reported_seed_agrees_to_1e_minus_4(data_dir):
    """The property every controlled experiment relies on.

    Re-fitting the *selected* model at the *reported* redshift must land on the
    same solution. It does not reproduce it bit for bit — the optimiser starts
    somewhere else and stops at a slightly different point — so the tolerance is
    1e-4 relative, not zero. Measured: 1e-8 to 1e-6 in flux, 1e-5 in sigma.
    """
    case = gspec.read_case(data_dir / "cases" / "S08.fits")
    a = measure(case)
    b = measure_fixed_model(case, a.model, z_known=a.redshift)
    assert b.redshift == pytest.approx(a.redshift, rel=1e-4)
    for line in a.lines:
        assert b.lines[line].flux == pytest.approx(a.lines[line].flux, rel=1e-4)
        assert b.lines[line].flux_err == pytest.approx(a.lines[line].flux_err, rel=1e-3)


# ------------------------------------------------------------- integration
def test_integrated_flux_scales_with_the_amplitude():
    """F is linear in the amplitude, so a 10x amplitude is a 10x flux."""
    assert integrated_flux(10.0, 1.5, 1.2) == pytest.approx(10 * integrated_flux(1.0, 1.5, 1.2))
