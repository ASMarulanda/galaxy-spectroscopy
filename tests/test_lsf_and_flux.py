"""Instrumental broadening and Gaussian integration.

Two relations carry the whole flux measurement:

    sigma_obs^2 = sigma_int^2 + sigma_LSF^2        (convolution of Gaussians)
    F = A * sigma_obs * sqrt(2 pi)                 (integral of the observed profile)

Both are checked here against an *analytic* Gaussian rather than against the
pipeline, so a change in the pipeline cannot make the test agree with itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from gspec import constants as K
from gspec.fluxes import integrated_flux
from gspec.lsf import gaussian, sigma_int_angstrom, sigma_lsf_angstrom, sigma_pix

SQRT2PI = math.sqrt(2.0 * math.pi)


def test_quadrature_identity():
    mu, sv, wd, dl = 6827.4, 85.0, 0.878, 1.5686
    expected = math.sqrt((mu * sv / K.C_KMS) ** 2 + (wd * dl) ** 2)
    assert sigma_pix(mu, sv, wd, dl) == pytest.approx(expected, rel=1e-15)


def test_linear_addition_would_be_wrong_by_forty_percent():
    """The failure mode, measured: adding widths inflates every flux."""
    mu, sv, wd, dl = 6827.4, 85.0, 0.878, 1.5686
    si, sl = sigma_int_angstrom(mu, sv), sigma_lsf_angstrom(wd, dl)
    quad = math.sqrt(si ** 2 + sl ** 2)
    linear = si + sl
    assert linear / quad == pytest.approx(1.396, abs=2e-3)
    # and the flux ratio is the width ratio, so the error is +39.6 %
    assert integrated_flux(1.0, si, sl) / (1.0 * quad * SQRT2PI) == pytest.approx(1.0)


def test_integral_of_an_analytic_gaussian_is_recovered():
    """F = A sigma sqrt(2 pi) against a numerical integral of the same profile."""
    mu, amp, sigma = 7000.0, 150.0, 2.3866
    w = np.linspace(mu - 40 * sigma, mu + 40 * sigma, 2_000_001)
    y = amp * np.exp(-0.5 * ((w - mu) / sigma) ** 2)
    numeric = np.trapz(y, w)          # numpy 1.26 is pinned; 2.x renames this
    analytic = amp * sigma * SQRT2PI
    assert numeric == pytest.approx(analytic, rel=1e-10)


def test_convolution_preserves_the_integral():
    """Observed profile integral == intrinsic flux, which is why the LSF is not deconvolved."""
    mu, flux, sv = 6827.4, 900.0, 85.0
    si = sigma_int_angstrom(mu, sv)
    sl = sigma_lsf_angstrom(0.878, 1.5686)
    so = math.sqrt(si ** 2 + sl ** 2)
    amp = flux / (so * SQRT2PI)          # as the generator does it
    w = np.linspace(mu - 30 * so, mu + 30 * so, 2_000_001)
    y = amp * np.exp(-0.5 * ((w - mu) / so) ** 2)
    numeric = np.trapz(y, w)          # numpy 1.26 is pinned; 2.x renames this
    assert numeric == pytest.approx(flux, rel=1e-9)
    assert integrated_flux(amp, si, sl) == pytest.approx(flux, rel=1e-12)


def test_peak_is_not_flux_and_the_ratio_is_the_width():
    """The density/integral confusion: the ratio is sqrt(2 pi) * sigma_obs."""
    si, sl = 1.9353, 1.3966
    so = math.sqrt(si ** 2 + sl ** 2)
    amp = 150.3894
    assert integrated_flux(amp, si, sl) / amp == pytest.approx(SQRT2PI * so, rel=1e-12)
    assert SQRT2PI * so == pytest.approx(5.98, abs=0.01)   # NOT "two orders of magnitude"


def test_lsf_scales_with_wavelength():
    """sigma_LSF in Angstrom grows across the band; a global value would be wrong."""
    dlam_3800, dlam_9200 = 3800.0 * math.log(10) * 9.99e-5, 9200.0 * math.log(10) * 9.99e-5
    s_lo, s_hi = sigma_lsf_angstrom(1.0, dlam_3800), sigma_lsf_angstrom(1.0, dlam_9200)
    assert s_hi / s_lo == pytest.approx(9200.0 / 3800.0, rel=1e-12)


def test_gaussian_is_unit_height():
    w = np.linspace(6000.0, 7000.0, 1001)
    g = gaussian(w, 6500.0, 2.0)
    assert g.max() == pytest.approx(1.0, rel=1e-15)
    assert g[np.argmin(np.abs(w - 6500.0))] == pytest.approx(1.0, rel=1e-12)
