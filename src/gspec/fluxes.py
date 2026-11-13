"""Integrated line fluxes, and the per-component breakdown behind them.

The scientific core of the package. A reported flux is the integral of the
*observed* line profile:

    F = A * sqrt(sigma_int^2 + sigma_LSF^2) * sqrt(2 pi)

Three properties, each of which a reader should be able to defend:

1. **Quadrature, not addition.** Intrinsic and instrumental profiles are both
   Gaussian, so their convolution has the sum of variances.
2. **The integral is of the observed profile.** Convolution preserves area, so
   this is the true flux of the component; deconvolving the LSF first would be
   wrong.
3. **Sum over components.** A reported line is the sum of its narrow and broad
   components, so Halpha is total narrow + broad by definition.
"""

from __future__ import annotations

import math

import numpy as np

from . import constants as K
from .io import Case

SQRT2PI = math.sqrt(2.0 * math.pi)

#: Narrow-emission parameter that carries each reported line's velocity width.
#: Only used to report ``sigma_kms``; it does not enter any flux.
SIGMA_PARAM_FOR_LINE = {
    "Halpha": K.NAME_NARROW_A,
    "NII_6583": K.NAME_NARROW_A,
    "Hbeta": K.NAME_NARROW_B,
    "OIII_5007": K.NAME_NARROW_O,
}

#: Components whose amplitude, if non-zero, means the line has a broad part.
BROAD_KEYS = ("Ha_b", "Hb_b")


def integrated_flux(amplitude: float, sigma_int: float, sigma_lsf: float) -> float:
    """Integral of a Gaussian of height ``amplitude`` and observed width sqrt(int^2+lsf^2)."""
    return amplitude * math.sqrt(sigma_int ** 2 + sigma_lsf ** 2) * SQRT2PI


def components(
    p: np.ndarray,
    idx: dict[str, int],
    case: Case,
    complexes,
    line: str,
) -> list[dict]:
    """Every free component contributing to ``line``, with its own integral.

    Returned for diagnostics and plots: the sum of these ``flux`` values is what
    the pipeline reports, so a reader can see *which* component moved when a
    number changes.
    """
    out = []
    for cx in complexes:
        for e in cx.emis:
            if e.line != line or e.tie is not None:
                continue
            mu = e.rest * (1.0 + p[idx[K.NAME_Z]])
            sig_int = mu * math.exp(p[idx[e.sigma]]) / K.C_KMS
            j = int(np.argmin(np.abs(case.wave - mu)))
            sig_lsf = float(case.wdisp[j] * case.dlam[j])
            amp = p[idx[f"amp_{e.key}"]]
            out.append(dict(
                key=e.key,
                rest=e.rest,
                mu=mu,
                amplitude=float(amp),
                sigma_v_kms=float(math.exp(p[idx[e.sigma]])),
                sigma_int=sig_int,
                sigma_lsf=sig_lsf,
                sigma_obs=math.sqrt(sig_int ** 2 + sig_lsf ** 2),
                flux=integrated_flux(amp, sig_int, sig_lsf),
            ))
    return out


def line_total(p: np.ndarray, idx: dict[str, int], case: Case, complexes, line: str) -> float:
    """Total integrated flux of one reported line: the sum over its components."""
    return float(sum(c["flux"] for c in components(p, idx, case, complexes, line)))


def broad_amplitude(p: np.ndarray, idx: dict[str, int], line: str) -> float:
    """Amplitude of the broad component of ``line``, or 0.0 when there is none."""
    line_of = {"Ha_b": "Halpha", "Hb_b": "Hbeta"}
    amp = 0.0
    for key in BROAD_KEYS:
        if f"amp_{key}" in idx and line_of[key] == line:
            amp = float(p[idx[f"amp_{key}"]])
    return amp
