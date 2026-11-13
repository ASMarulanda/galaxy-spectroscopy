"""Covariance of the fitted parameters, and propagation to a line flux.

    Cov(p) = pinv(JtJ) * min(chi2/dof, COV_SCALE_CAP)
    sigma_F^2 = grad(F)^T Cov(p) grad(F)

Two deliberate choices, both measured, both documented as limitations:

**The scaling cap.** For a least-squares problem with the residuals weighted by
sqrt(ivar), ``JtJ`` is the inverse covariance. Scaling it by ``chi2/dof`` is the
standard correction for underestimated errors — but on real spectra chi2/dof is
driven by *model mismatch* (absorption lines, sky residuals), not by photon
noise, so the full correction over-corrects. The cap accepts up to 50 % extra
variance and no more. Measured on this dataset it binds on 2 of 8 cases
(docs/methods.md §6).

**The pseudoinverse.** ``pinv`` truncates directions whose eigenvalue falls
below ``rcond * lambda_max`` (numpy default ``rcond=1e-15``), which is what keeps
a rank-deficient normal matrix from producing infinities. It is also why
``sigma_F`` is not env-ironment-stable when an amplitude sits on its bound: on
S08 the direction carrying 37 % of ``sigma_F(Hbeta)^2`` has an eigenvalue ~12
orders of magnitude below ``lambda_max``, so a different optimiser path moves it.
See docs/limitations.md L4 and validation/experiments/.

The gradient is a forward finite difference of the *same* flux function the
pipeline reports, so a change in the flux formula cannot silently leave the
uncertainty behind.
"""

from __future__ import annotations

import math

import numpy as np

from . import constants as K


def covariance(jac: np.ndarray, chi2: float, dof: int, cap: float = K.COV_SCALE_CAP) -> np.ndarray:
    """Scaled covariance of the fitted parameters.

    ``jac`` is the Jacobian of the weighted residuals, so ``jac.T @ jac`` is the
    information matrix.
    """
    scale = min(chi2 / dof, cap)
    return np.linalg.pinv(jac.T @ jac) * scale


def covariance_scale(chi2: float, dof: int, cap: float = K.COV_SCALE_CAP) -> float:
    """The factor actually applied, exposed so results can report it."""
    return min(chi2 / dof, cap)


def gradient(flux_of_p, p: np.ndarray, rel_step: float = 1e-6) -> np.ndarray:
    """Forward-difference gradient of a scalar flux function of the parameters."""
    f0 = flux_of_p(p)
    g = np.zeros(len(p))
    for j in range(len(p)):
        step = rel_step * max(abs(p[j]), 1e-3)
        pj = p.copy()
        pj[j] += step
        g[j] = (flux_of_p(pj) - f0) / step
    return g


def sigma_of_flux(flux_of_p, p: np.ndarray, cov: np.ndarray) -> tuple[float, float]:
    """Return ``(F, sigma_F)`` for a scalar flux function at parameters ``p``."""
    f = flux_of_p(p)
    g = gradient(flux_of_p, p)
    var = float(g @ cov @ g)
    return f, math.sqrt(max(var, 0.0))
