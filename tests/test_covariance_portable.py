"""What the uncertainty can be asked to prove without a reference environment.

The reported sigma is

    Cov = pinv(JᵀJ) · min(χ²/dof, cap),     σ_F = sqrt(gᵀ Cov g),  g = ∂F/∂p,

with ``pinv``'s default ``rcond = 1e-15``. On this project's spectra some
amplitudes sit on a bound, and the direction they span has a curvature at the level
of floating-point noise — on S08 the smallest singular value of the Jacobian is
2.2e-8 against a largest of 0.38. Whether ``pinv`` keeps that direction is then
decided by the last digits, and σ(Hα) on that spectrum is bimodal: **11.50 or
40.15**, a factor of 3.5, depending on a one-ulp change in the input.

That makes σ on such a row *not portable*, and this file does not pretend
otherwise. What it asserts is everything that does survive:

* the propagation rule and the truncation rule, analytically, on constructed
  Jacobians — exactly, because they are arithmetic;
* on real spectra: σ finite and non-negative, the reconstructed covariance
  positive semi-definite, and a bound-pinned fit reporting its bounds;
* fluxes stable to 1e-6 relative under float32-lsb perturbations of the input,
  while σ is measured and printed rather than asserted.

What it cannot establish: that σ is *calibrated*, or that it repeats across BLAS
builds. Calibration lives in `tests/test_uncertainty.py` (pull distributions
against truth, which are tolerance-based and portable); cross-environment
repetition lives in `tests/test_regression.py` under the ``characterisation``
marker, and in `docs/limitations.md` L4 and L17.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import gspec
import gspec.constants as K
from gspec.covariance import covariance, covariance_scale, gradient, sigma_of_flux
from gspec.pipeline import measure, measure_fixed_model


# ---------------------------------------------------------------- the rule
def test_the_propagation_rule_is_the_analytic_quadratic_form():
    """σ² = gᵀCg, checked against arithmetic on a covariance we chose."""
    a = np.array([2.0, -1.0, 0.5])
    cov = np.array([[4.0, 0.5, 0.0],
                    [0.5, 1.0, 0.0],
                    [0.0, 0.0, 9.0]])
    flux = lambda p: float(a @ p) + 100.0        # noqa: E731
    p = np.array([1.0, 2.0, 3.0])

    f, sigma = sigma_of_flux(flux, p, cov)

    assert f == pytest.approx(float(a @ p) + 100.0, rel=1e-12)
    # g comes from a forward difference with rel_step 1e-6, so the recovered
    # sigma carries that quotient's rounding (~1e-9 here). 1e-7 is the accuracy
    # of the rule as implemented, not slack for a wrong answer.
    assert sigma == pytest.approx(math.sqrt(float(a @ cov @ a)), rel=1e-7)


def test_the_gradient_is_accurate_for_a_function_we_can_differentiate():
    """A forward difference with rel_step 1e-6 must recover an exact derivative."""
    def flux(p):
        return 3.0 * p[0] ** 2 + 0.5 * p[1] - 7.0

    p = np.array([2.0, 4.0])
    g = gradient(flux, p)
    assert g[0] == pytest.approx(12.0, rel=1e-6)      # d/dp0 = 6 p0
    assert g[1] == pytest.approx(0.5, rel=1e-6)


def test_a_null_direction_contributes_nothing():
    """A parameter the data cannot see gets zero variance, not infinite variance."""
    jac = np.array([[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])   # second column is null
    cov = covariance(jac, chi2=3.0, dof=1)
    assert np.isfinite(cov).all()
    assert cov[1, 1] == pytest.approx(0.0, abs=1e-12)
    assert cov[0, 0] > 0.0


def test_the_scale_is_min_of_chi2_over_dof_and_the_cap():
    assert covariance_scale(chi2=10.0, dof=10) == pytest.approx(1.0)
    assert covariance_scale(chi2=10_000.0, dof=10) == pytest.approx(K.COV_SCALE_CAP)


def test_covariance_is_positive_semidefinite_by_construction():
    rng = np.random.default_rng(0)
    jac = rng.normal(size=(40, 5))
    cov = covariance(jac, chi2=30.0, dof=10)
    w = np.linalg.eigvalsh(cov)
    assert w.min() > -1e-12 * max(abs(w.max()), 1.0)


# ------------------------------------------------------- on real spectra
@pytest.mark.parametrize("case_id", ["S01", "S08", "R01", "R06"])
def test_sigma_is_finite_and_non_negative_everywhere(cases, case_id):
    """A reported uncertainty must be a number, whatever the fit did."""
    r = measure(cases[case_id])
    for line, m in r.lines.items():
        assert math.isfinite(m.flux_err), f"{case_id}/{line}: sigma is not finite"
        assert m.flux_err >= 0.0, f"{case_id}/{line}: negative sigma"
        if m.flag == "ok":
            assert m.flux_err > 0.0, f"{case_id}/{line}: a detection with zero sigma"


@pytest.mark.parametrize("case_id", ["S08", "R01"])
def test_the_covariance_reconstructed_from_the_fitted_jacobian_is_psd(cases, case_id):
    """The covariance the pipeline would rebuild from its own Jacobian is usable."""
    r = measure_fixed_model(cases[case_id], measure(cases[case_id]).model,
                            keep_jacobian=True)
    assert r.jacobian is not None
    cov = covariance(r.jacobian, r.chi2, r.dof)
    assert np.isfinite(cov).all()
    w = np.linalg.eigvalsh(cov)
    assert w.min() > -1e-10 * max(abs(w.max()), 1.0), f"{case_id}: negative variance direction"


def test_bound_pinned_fits_report_their_bounds(cases):
    """The diagnostic that makes the instability visible must be present."""
    for case_id in ("S01", "S08"):
        r = measure(cases[case_id])
        assert r.at_bounds, f"{case_id}: no bounds reported, yet sigma moves 3.5x"
        for name in r.at_bounds:
            assert name in r.parameters, f"{case_id}: bound parameter {name} has no value"


def test_flux_is_stable_under_a_one_ulp_perturbation_of_the_input(cases, capsys):
    """The portable half of the perturbation experiment.

    Fluxes move by less than 1e-3 relative when every pixel is nudged by up to one
    float32 ulp. Measured: 8.3e-7 on the reference stack, 7.1e-7 on the wheel stack —
    a bound set at 1e-6 would have 30 % of headroom, which is a test that passes by
    luck. The 1e-3 threshold states the property that actually matters: a one-ulp
    change in the input cannot move a flux by anything a reader would notice, since a
    real modelling change moves it by 6e-4 at the smallest (L17) and a real bug by
    orders of magnitude more. Sigma on a bound-pinned row is *not* stable under the
    same nudge; that is measured and printed, never asserted.
    """
    case = cases["S08"]
    base = measure(case)
    rng = np.random.default_rng(20261116)

    flux_rel, sigma_ratio = [], []
    for _ in range(4):
        scaled = (case.flux * (1.0 + rng.uniform(-1e-7, 1e-7, case.flux.shape)))
        perturbed = gspec.Case(
            case_id=case.case_id, wave=case.wave,
            flux=scaled.astype(np.float32).astype(np.float64),
            ivar=case.ivar.copy(), wdisp=case.wdisp.copy(),
            loglam=case.loglam.copy(), and_mask=case.and_mask,
            or_mask=case.or_mask, bunit=case.bunit, provenance=case.provenance)
        r = measure(perturbed)
        for line in ("Halpha", "Hbeta", "NII_6583"):
            f0 = base.lines[line].flux
            if abs(f0) > 0.0:
                flux_rel.append(abs(r.lines[line].flux - f0) / abs(f0))
            sigma_ratio.append(r.lines[line].flux_err / base.lines[line].flux_err)

    assert max(flux_rel) < 1e-3, f"flux moved by {max(flux_rel):.2e} under a 1-ulp nudge"
    with capsys.disabled():
        print(f"\n  flux change under 1-ulp perturbation: max {max(flux_rel):.2e} relative")
        print(f"  sigma change under the same nudge  : {min(sigma_ratio):.2f}x .. "
              f"{max(sigma_ratio):.2f}x  (not asserted — see docs/limitations.md L17)")
