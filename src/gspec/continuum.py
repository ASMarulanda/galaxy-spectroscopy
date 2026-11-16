"""The continuum model: a local polynomial, attenuated by stellar absorption.

The measured spectrum is

    f(lambda) = C(lambda) * prod_i [1 - d_i G_i(lambda)]  +  sum_j A_j G_j(lambda)

This module owns ``C`` and the absorption product. It knows nothing about
emission lines and nothing about the optimiser.

Why a polynomial
----------------
This is a *local approximation*, not a stellar-population fit. It is adequate
over a 370 A window, and its error is absorbed into the reported uncertainty and
into the BIC model-selection step. The synthetic generator deliberately uses a
different basis (three power laws times a 12-line absorption table) so that
assuming the generator's model cannot pass.

The identifiability decision
----------------------------
The absorption width is not free: it is ``sigma_narrow,B * exp(abs_ratio)`` with
``abs_ratio >= 1.15``. That is a physical prior — stellar absorption is broader
than ionised-gas emission, because stars trace the whole galaxy while H II
regions are kinematically cold — and it is what makes the Hbeta depth
identifiable from a single Balmer line. Leaving the width free let the fit trade
depth against a spurious broad component and biased Hbeta by +10..+80 %.

.. warning::
   On the ten synthetic cases in this repository the BIC selection rejects the
   absorption term altogether (`docs/methods.md` §4.2). The model is kept because
   it is selected on 5 of 6 real spectra; the claim that it *improves* Hbeta is
   not supported and is not made anywhere in this repository.
"""

from __future__ import annotations

import math

import numpy as np

from . import constants as K
from .lines import Complex
from .lsf import gaussian, sigma_pix


def polynomial(
    w: np.ndarray,
    cx: Complex,
    p: np.ndarray,
    idx: dict[str, int],
    ndeg: int,
    ci: int,
) -> np.ndarray:
    """Rest-frame polynomial continuum for complex ``ci`` on pixels ``w``.

    ``x`` is a normalised rest-frame coordinate so that the coefficients are
    comparable between the two windows regardless of their width:

        x = (lambda/(1+z) - centre) / width
    """
    x = (w / (1.0 + p[idx[K.NAME_Z]]) - cx.centre) / max(cx.width, 1.0)
    base = ci * (ndeg + 1)
    cont = np.zeros_like(w)
    for d in range(ndeg + 1):
        cont = cont + p[idx[f"c{base + d}"]] * x ** d
    return cont


def apply_absorption(
    cont: np.ndarray,
    w: np.ndarray,
    cx: Complex,
    p: np.ndarray,
    idx: dict[str, int],
    wdisp: np.ndarray,
    dlam: np.ndarray,
) -> np.ndarray:
    """Apply the multiplicative absorption product to ``cont``, in place of it.

    The absorption width is tied to the narrow emission width of the Hbeta
    complex, widened by ``exp(abs_ratio)``. The G-band uses its own scale
    because it is a molecular band, not a Balmer line, and its depth is not
    tied to the Balmer pattern.
    """
    sv = math.exp(p[idx[K.NAME_NARROW_B]]) * math.exp(p[idx[K.NAME_ABS_RATIO]])
    for mu_rest, weight, wmult in cx.absorbers:
        mu = mu_rest * (1.0 + p[idx[K.NAME_Z]])
        sig = wmult * sigma_pix(mu, sv, wdisp, dlam)
        scale = p[idx[K.NAME_GBAND]] if wmult > 1.5 else p[idx[K.NAME_ABS_SCALE]]
        cont = cont * (1.0 - scale * weight * gaussian(w, mu, sig))
    return cont


def continuum_at(
    w: np.ndarray,
    cx: Complex,
    ci: int,
    p: np.ndarray,
    idx: dict[str, int],
    ndeg: int,
    use_absorption: bool,
    wdisp: np.ndarray,
    dlam: np.ndarray,
) -> np.ndarray:
    """Polynomial continuum with the absorption product applied, if enabled."""
    cont = polynomial(w, cx, p, idx, ndeg, ci)
    if use_absorption and cx.absorbers:
        cont = apply_absorption(cont, w, cx, p, idx, wdisp, dlam)
    return cont
