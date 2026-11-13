"""Instrumental line-spread function: one relation, stated once.

Both the intrinsic profile and the instrumental response are approximately
Gaussian, so the observed profile is their convolution and the **variances**
add:

    sigma_obs^2 = sigma_int^2 + sigma_LSF^2

The integral of a Gaussian is preserved under convolution, so the flux of the
*observed* profile is the true flux of the component. The LSF must therefore
never be deconvolved before integrating.

Two traps this module exists to prevent:

* adding widths linearly instead of in quadrature (+39.6 % on S01 Halpha),
* using a global delta-lambda instead of the local one (up to ~10 % across the
  3800-9200 A band, because delta-lambda = lambda * ln10 * dlog10(lambda)).
"""

from __future__ import annotations

import numpy as np

from . import constants as K


def sigma_pix(mu: float, sigma_v_kms: float, wdisp: np.ndarray, dlam: np.ndarray) -> np.ndarray:
    """Observed width in pixels for a component at observed wavelength ``mu``.

    Parameters
    ----------
    mu
        Observed (not rest) line centre, Angstrom.
    sigma_v_kms
        Intrinsic velocity dispersion, km/s.
    wdisp, dlam
        Instrumental sigma in pixels and the local pixel width in Angstrom,
        both arrays on the pixels of the window being modelled.

    Notes
    -----
    This two-line quadrature is the only routine the independent implementation
    shares with this package, because both express the same physical relation.
    """
    return np.sqrt((mu * sigma_v_kms / K.C_KMS) ** 2 + (wdisp * dlam) ** 2)


def sigma_int_angstrom(mu: float, sigma_v_kms: float) -> float:
    """Intrinsic width in Angstrom: a velocity width scaled by lambda_obs/c."""
    return float(mu * sigma_v_kms / K.C_KMS)


def sigma_lsf_angstrom(wdisp_pix: float, dlam_angstrom: float) -> float:
    """Instrumental width in Angstrom at one pixel: wdisp[pixels] * dlambda."""
    return float(wdisp_pix * dlam_angstrom)


def gaussian(w: np.ndarray, mu: float, sigma: np.ndarray | float) -> np.ndarray:
    """Unit-height Gaussian, the profile shape used for every component."""
    return np.exp(-0.5 * ((w - mu) / sigma) ** 2)
