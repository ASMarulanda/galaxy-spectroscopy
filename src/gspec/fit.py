"""The bounded non-linear fit: parameters, bounds, solve, and results.

This is the only place the optimiser is called. It builds the parameter vector
in a fixed order, bounds it, hands the residual function to
``scipy.optimize.least_squares`` and turns the solution into a
:class:`~gspec.results.SpectrumResult`.

Parameter order (fixed, and part of the reproducibility contract)::

    z  narrowA  narrowB  narrowO  broad
    [abs_scale  abs_ratio  gband]          only when absorption is fitted
    amp_<key> for every untied emission component, in table order
    c0 .. c_{n_cont-1}                     (ndeg+1) per complex

Why the fit is joint
--------------------
Halpha and Hbeta are fitted **simultaneously with one common redshift**. A
single complex cannot separate redshift from line width; two complexes 1700 A
apart can, because the redshift shifts both by the same fraction while the
widths are independent. The continuum is allowed to differ between the two
windows (independent polynomials), which absorbs the real shape difference
without forcing a global model.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import least_squares

from . import constants as K
from .covariance import covariance, covariance_scale, sigma_of_flux
from .flags import classify
from .fluxes import SIGMA_PARAM_FOR_LINE, broad_amplitude, components, line_total
from .io import Case, clean
from .lines import Complex, alpha_complex, beta_complex
from .model import residuals
from .redshift import redshift_peaks
from .results import ComponentMeasurement, LineMeasurement, ModelConfig, SpectrumResult

#: Solver settings, preserved from the validated implementation. `x_scale="jac"`
#: makes the step size scale with each parameter's own derivative, which matters
#: because the vector mixes a redshift (~0.05), log-widths (~4.5), amplitudes
#: (~10-1000) and continuum coefficients (~1-100).
X_SCALE = "jac"
MAX_NFEV = 3000
FTOL = 1e-10
XTOL = 1e-10

#: Minimum usable pixels in a window, below which the fit is not identifiable.
MIN_WINDOW_PIXELS = 40


def build_complexes(config: ModelConfig) -> list[Complex]:
    """Instantiate the window/line table for one model configuration."""
    return [
        alpha_complex(pad=config.pad),
        beta_complex(window=config.hb_window,
                     gamma_weight=config.gamma_weight,
                     absorbers=config.absorbers_hb),
    ]


def _select_windows(case: Case, complexes: list[Complex], z0: float) -> list[np.ndarray] | None:
    good = clean(case)
    windows = []
    for cx in complexes:
        m = good & (case.wave > cx.lo_rest * (1 + z0)) & (case.wave < cx.hi_rest * (1 + z0))
        if m.sum() < MIN_WINDOW_PIXELS:
            return None
        windows.append(m)
    return windows


def build_parameters(
    case: Case,
    complexes: list[Complex],
    windows: list[np.ndarray],
    config: ModelConfig,
    z0: float,
) -> tuple[list[str], dict[str, int], np.ndarray, np.ndarray, np.ndarray]:
    """Assemble names, index map, starting point and bounds, in the fixed order."""
    ndeg = config.ndeg
    use_abs = config.use_absorption

    names = [K.NAME_Z, K.NAME_NARROW_A, K.NAME_NARROW_B, K.NAME_NARROW_O, K.NAME_BROAD]
    if use_abs:
        names += [K.NAME_ABS_SCALE, K.NAME_ABS_RATIO, K.NAME_GBAND]

    free_keys: list[str] = []
    for cx in complexes:
        for e in cx.emis:
            if e.tie is None and e.key not in free_keys:
                free_keys.append(e.key)
    names += [f"amp_{k}" for k in free_keys]

    ncont = (ndeg + 1) * len(complexes)
    names += [f"c{i}" for i in range(ncont)]
    idx = {n: i for i, n in enumerate(names)}

    p0 = np.zeros(len(names))
    lo = np.full(len(names), -np.inf)
    hi = np.full(len(names), np.inf)

    p0[idx[K.NAME_Z]] = z0
    lo[idx[K.NAME_Z]], hi[idx[K.NAME_Z]] = z0 - 0.0035, z0 + 0.0035

    for name, v0 in ((K.NAME_NARROW_A, 90.0), (K.NAME_NARROW_B, 90.0), (K.NAME_NARROW_O, 90.0)):
        p0[idx[name]] = math.log(v0)
        lo[idx[name]], hi[idx[name]] = math.log(K.SIGMA_MIN_KMS), math.log(K.SIGMA_MAX_KMS)

    # One width for both permitted broad lines: Halpha and Hbeta come from the
    # same broad-line region. Bounding it is essential -- an unbounded
    # 2500 km/s "broad" component is indistinguishable from a slowly varying
    # continuum over a 370 A window and ran to its bound on the broad-line AGN.
    p0[idx[K.NAME_BROAD]] = math.log(700.0)
    lo[idx[K.NAME_BROAD]], hi[idx[K.NAME_BROAD]] = (
        math.log(K.BROAD_MIN_KMS), math.log(K.BROAD_MAX_KMS))

    if use_abs:
        p0[idx[K.NAME_ABS_SCALE]] = 0.3
        lo[idx[K.NAME_ABS_SCALE]], hi[idx[K.NAME_ABS_SCALE]] = 0.0, 1.5
        p0[idx[K.NAME_ABS_RATIO]] = math.log(1.8)
        lo[idx[K.NAME_ABS_RATIO]] = math.log(K.ABS_RATIO_MIN)
        hi[idx[K.NAME_ABS_RATIO]] = math.log(K.ABS_RATIO_MAX)
        p0[idx[K.NAME_GBAND]] = 0.15
        lo[idx[K.NAME_GBAND]], hi[idx[K.NAME_GBAND]] = 0.0, 1.0

    for k in free_keys:
        p0[idx[f"amp_{k}"]] = 5.0
        lo[idx[f"amp_{k}"]] = 0.0
        if k == "Hb_b" and not config.allow_broad_b:
            hi[idx[f"amp_{k}"]] = 1e-9
            p0[idx[f"amp_{k}"]] = 0.0

    # NOTE (preserved quirk): the reference implementation initialises and
    # bounds only c0 and c1, i.e. the constant term of the first complex and the
    # *slope* of the first complex -- not the constant term of the second. All
    # other continuum coefficients start at 0 and are unbounded. It is kept
    # verbatim so the port reproduces the baseline; see docs/limitations.md L7.
    for i in range(len(complexes)):
        p0[idx[f"c{i}"]] = float(np.median(case.flux[windows[i]]))
        lo[idx[f"c{i}"]] = -50.0
        hi[idx[f"c{i}"]] = 500.0

    return names, idx, p0, lo, hi


def at_bounds(
    p: np.ndarray,
    names: list[str],
    lo: np.ndarray,
    hi: np.ndarray,
    rel_tol: float = 1e-9,
) -> list[str]:
    """Names of parameters sitting on a bound (to a relative tolerance).

    Reported because a bound-pinned amplitude is exactly what makes the
    covariance fragile (docs/limitations.md L4), and because a reader should be
    able to tell "the model put this at zero" from "the model estimated a small
    number".
    """
    pinned = []
    for n, i in ((n, i) for i, n in enumerate(names)):
        span = hi[i] - lo[i] if np.isfinite(hi[i] - lo[i]) else 1.0
        tol = rel_tol * max(span, abs(p[i]), 1e-12)
        if np.isfinite(lo[i]) and abs(p[i] - lo[i]) <= tol:
            pinned.append(n)
        elif np.isfinite(hi[i]) and abs(p[i] - hi[i]) <= tol:
            pinned.append(n)
    return pinned


def fit_case(
    case: Case,
    z_init: float | None = None,
    config: ModelConfig | None = None,
    keep_jacobian: bool = False,
) -> SpectrumResult:
    """Fit one spectrum at (near) a given redshift seed.

    Returns a :class:`SpectrumResult`. On failure the result carries
    ``error`` and empty ``lines``; it never raises for a data-level problem.
    """
    config = config or ModelConfig()
    complexes = build_complexes(config)

    if z_init is None:
        seeds = redshift_peaks(case, n=1)
        z_init = seeds[0] if seeds else 0.05
    z0 = float(z_init)

    windows = _select_windows(case, complexes, z0)
    if windows is None:
        return SpectrumResult(case_id=case.case_id, redshift=z0,
                              redshift_err=float("nan"), lines={},
                              model=config, model_selected=config.label(),
                              error="window too small")

    names, idx, p0, lo, hi = build_parameters(case, complexes, windows, config, z0)
    use_abs = config.use_absorption

    def resid(p):
        return residuals(p, idx, case, complexes, windows, config.ndeg, use_abs)

    try:
        sol = least_squares(
            resid, p0, bounds=(lo, hi),
            x_scale=X_SCALE, max_nfev=MAX_NFEV, ftol=FTOL, xtol=XTOL,
        )
    except Exception as exc:                             # noqa: BLE001
        return SpectrumResult(case_id=case.case_id, redshift=z0,
                              redshift_err=float("nan"), lines={},
                              model=config, model_selected=config.label(),
                              error=str(exc))

    p = sol.x
    ndata = int(sum(w.sum() for w in windows))
    dof = max(ndata - len(p), 1)
    chi2 = float(np.sum(sol.fun ** 2))
    cov = covariance(sol.jac, chi2, dof)

    lines: dict[str, LineMeasurement] = {}
    for line in K.LINE_NAMES:
        f, ferr = sigma_of_flux(
            lambda pv, _l=line: line_total(pv, idx, case, complexes, _l), p, cov)
        amp_broad = broad_amplitude(p, idx, line)
        flag = classify(f, ferr, amp_broad)
        parts = [
            ComponentMeasurement(
                key=c["key"],
                rest=c["rest"],
                amplitude=c["amplitude"],
                sigma_v_kms=c["sigma_v_kms"],
                sigma_int=c["sigma_int"],
                sigma_lsf=c["sigma_lsf"],
                sigma_obs=c["sigma_obs"],
                flux=c["flux"],
            )
            for c in components(p, idx, case, complexes, line)
        ]
        lines[line] = LineMeasurement(
            line=line,
            flux=float(f),
            flux_err=float(ferr),
            flag=flag,
            sigma_kms=float(math.exp(p[idx[SIGMA_PARAM_FOR_LINE[line]]])),
            broad_amplitude=float(amp_broad),
            components=parts,
        )

    return SpectrumResult(
        case_id=case.case_id,
        redshift=float(p[idx[K.NAME_Z]]),
        redshift_err=float(math.sqrt(max(cov[idx[K.NAME_Z], idx[K.NAME_Z]], 0.0))),
        lines=lines,
        chi2=chi2,
        dof=dof,
        n_par=len(p),
        covariance_scale=covariance_scale(chi2, dof),
        model=config,
        model_selected=config.label(),
        parameters={n: float(v) for n, v in zip(names, p)},
        at_bounds=at_bounds(p, names, lo, hi),
        jacobian=np.array(sol.jac, float) if keep_jacobian else None,
    )
