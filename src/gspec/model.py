"""Composing the full spectral model, and the weighted residual vector.

The model is a pure function of (parameters, case, configuration). Nothing here
optimises, and nothing here decides which pixels are usable — the caller passes
already-masked windows in.

    f(lambda) = C(lambda) * prod_i [1 - d_i G_i(lambda)]  +  sum_j A_j G_j(lambda)

The residual returned to the optimiser is weighted by ``sqrt(ivar)``, so a
least-squares solve on it is exactly a chi-square minimisation.
"""

from __future__ import annotations

import math

import numpy as np

from . import constants as K
from .continuum import continuum_at
from .io import Case
from .lines import Complex
from .lsf import gaussian, sigma_pix


def amplitudes(cx: Complex, p: np.ndarray, idx: dict[str, int]) -> dict[str, float]:
    """Free amplitudes for one complex, then the tied ones resolved from them.

    Ties are applied after the free amplitudes are read so that a tie can point
    at a component defined later in the table.
    """
    amps = {e.key: p[idx[f"amp_{e.key}"]] for e in cx.emis if e.tie is None}
    for e in cx.emis:
        if e.tie is not None:
            src, factor = e.tie
            amps[e.key] = amps[src] * factor
    return amps


def model_segments(
    p: np.ndarray,
    idx: dict[str, int],
    case: Case,
    complexes: list[Complex],
    windows: list[np.ndarray],
    ndeg: int,
    use_absorption: bool,
) -> list[np.ndarray]:
    """One model array per complex, evaluated on that complex's window."""
    wave = case.wave
    segs = [
        continuum_at(wave[windows[ci]], complexes[ci], ci, p, idx, ndeg,
                     use_absorption, case.wdisp[windows[ci]], case.dlam[windows[ci]])
        for ci in range(len(complexes))
    ]
    for ci, cx in enumerate(complexes):
        w = wave[windows[ci]]
        wd, dl = case.wdisp[windows[ci]], case.dlam[windows[ci]]
        amps = amplitudes(cx, p, idx)
        for e in cx.emis:
            mu = e.rest * (1.0 + p[idx[K.NAME_Z]])
            sig = sigma_pix(mu, math.exp(p[idx[e.sigma]]), wd, dl)
            segs[ci] = segs[ci] + amps[e.key] * gaussian(w, mu, sig)
    return segs


def residuals(
    p: np.ndarray,
    idx: dict[str, int],
    case: Case,
    complexes: list[Complex],
    windows: list[np.ndarray],
    ndeg: int,
    use_absorption: bool,
) -> np.ndarray:
    """Weighted residual vector, concatenated over the complexes."""
    segs = model_segments(p, idx, case, complexes, windows, ndeg, use_absorption)
    return np.concatenate([
        (case.flux[windows[ci]] - segs[ci]) * np.sqrt(case.ivar[windows[ci]])
        for ci in range(len(complexes))
    ])


def model_on_pixels(
    p: np.ndarray,
    idx: dict[str, int],
    case: Case,
    complexes: list[Complex],
    windows: list[np.ndarray],
    ndeg: int,
    use_absorption: bool,
) -> np.ndarray:
    """Model flux density evaluated on an arbitrary pixel selection.

    Used by the plots and by the diagnostics, not by the optimiser: it evaluates
    the same parameters on all *clean* pixels of the case, which is what a
    reader wants to see overlaid on the data.
    """
    segs = model_segments(p, idx, case, complexes, windows, ndeg, use_absorption)
    out = np.full_like(case.flux, np.nan)
    for ci in range(len(complexes)):
        out[windows[ci]] = segs[ci]
    return out
