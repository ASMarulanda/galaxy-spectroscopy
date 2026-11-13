"""Redshift estimation: a matched-filter candidate generator, then a chi-square polish.

Two stages, for one reason. A raw chi-square scan over the whole spectrum is
insensitive to the *identity* of a line and aliases onto the wrong member of the
Halpha/[N II] doublet at low S/N — a 943 km/s trap. A matched filter over the
known optical lines is insensitive to the slowly varying continuum and does not
alias; but a filter estimates a location coarsely. So the filter proposes
candidates and a local linear (non-negative least squares) chi-square scan
refines each one.

The linear scan exploits the structure of the problem: at a fixed trial redshift
the model is **linear in every amplitude**, so the amplitudes have a closed-form
non-negative solution and only ``z`` and the widths need a search.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import nnls

from . import constants as K
from .io import Case, clean
from .lines import Complex, default_complexes
from .lsf import gaussian, sigma_pix

#: Trial-redshift grid used by the filter and by the polish.
Z_MIN_DEFAULT = 0.008
Z_MAX_DEFAULT = 0.102

#: Relative strengths of the optical lines, used to weight the matched filter.
#: These are typical star-forming ratios, not a physical model of any object:
#: the filter only has to prefer the right redshift, not to fit anything.
ZFILTER = (
    (6564.61, 1.00), (6585.27, 0.35), (6549.86, 0.12),
    (4862.68, 0.35), (5008.24, 0.32), (4960.29, 0.11),
    (4341.69, 0.16), (6718.29, 0.13), (6732.67, 0.11),
    (3727.09, 0.45), (3729.88, 0.45),
)


def design_at(
    case: Case,
    z: float,
    complexes: list[Complex],
    sigma_fixed: float = 110.0,
    ndeg: int = 2,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]] | None:
    """Weighted design blocks for a linear amplitude solve at one trial redshift.

    Returns one ``(A, y, w)`` per complex, where ``A @ a`` is the model and ``w``
    are the ``sqrt(ivar)`` weights; ``None`` when a window has fewer than 40
    usable pixels (the fit is not identifiable there).
    """
    good = clean(case)
    wave, flux, ivar = case.wave, case.flux, case.ivar
    wd, dl = case.wdisp, case.dlam
    blocks = []
    for cx in complexes:
        m = good & (wave > cx.lo_rest * (1 + z)) & (wave < cx.hi_rest * (1 + z))
        if m.sum() < 40:
            return None
        w = wave[m]
        x = (w / (1.0 + z) - cx.centre) / max(cx.width, 1.0)
        cols = [np.ones_like(x)] + [x ** d for d in range(1, ndeg + 1)]
        for e in cx.emis:
            if e.tie is not None:
                continue
            mu = e.rest * (1.0 + z)
            sig = sigma_pix(mu, sigma_fixed, wd[m], dl[m])
            cols.append(gaussian(w, mu, sig))
        blocks.append((np.vstack(cols).T, flux[m], np.sqrt(ivar[m])))
    return blocks


def linear_redshift_scan(
    case: Case,
    zlo: float,
    zhi: float,
    nz: int,
    widths: tuple[float, ...] = (150.0,),
    ndeg: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Local chi-square scan of a linear (NNLS) amplitude fit over ``[zlo, zhi]``.

    The complexes used here are always the reference (m1a) ones, even when the
    surrounding fit runs the m1b variant. That is deliberate: the scan is a
    *redshift* estimator, and changing its windows for a variant would confound
    "which model fits the lines better" with "which seed the search found".
    """
    complexes = default_complexes()
    zs = np.linspace(zlo, zhi, nz)
    chi2 = np.full(nz, np.inf)
    for iz, z in enumerate(zs):
        best = np.inf
        for sv in widths:
            blocks = design_at(case, z, complexes, sigma_fixed=sv, ndeg=ndeg)
            if blocks is None:
                continue
            total = 0.0
            for A, y, w in blocks:
                Aw, yw = A * w[:, None], y * w
                try:
                    a, _ = nnls(Aw, yw)
                except Exception:                     # noqa: BLE001
                    total = np.inf
                    break
                total += float(np.sum((yw - Aw @ a) ** 2))
            best = min(best, total)
        chi2[iz] = best
    return zs, chi2


def matched_filter_scan(
    case: Case,
    zlo: float = Z_MIN_DEFAULT,
    zhi: float = Z_MAX_DEFAULT,
    nz: int = 1400,
    sig_a: float = 2.2,
) -> tuple[np.ndarray, np.ndarray]:
    """Weighted matched filter over the strong optical lines.

    For each trial redshift and each template line, the local continuum is
    estimated from sidebands 2-5 pixels away on both sides, and the filter
    response is

        S(z) += w_l * sum_i (f_i - median(side)) K_i ivar_i / sqrt(sum_i K_i^2 ivar_i)

    with ``K`` a Gaussian of sigma ``sig_a`` pixels. The sideband subtraction is
    what makes the score insensitive to the continuum shape.
    """
    m = clean(case)
    wave, flux, ivar = case.wave[m], case.flux[m], case.ivar[m]
    zs = np.linspace(zlo, zhi, nz)
    score = np.full(nz, -np.inf)
    if len(wave) < 200:
        return zs, score
    for iz, z in enumerate(zs):
        s = 0.0
        for l0, wgt in ZFILTER:
            lc = l0 * (1.0 + z)
            if lc < wave[0] + 25 or lc > wave[-1] - 25:
                continue
            k = np.abs(wave - lc) < 4 * sig_a
            if k.sum() < 3:
                continue
            blue = (wave > lc - 5 * sig_a) & (wave < lc - 2 * sig_a)
            red = (wave > lc + 2 * sig_a) & (wave < lc + 5 * sig_a)
            side = np.concatenate([flux[blue], flux[red]])
            if len(side) < 4:
                continue
            y = flux[k] - np.median(side)
            kern = np.exp(-0.5 * ((wave[k] - lc) / sig_a) ** 2)
            den = math.sqrt(float(np.sum(kern ** 2 * ivar[k])))
            if den > 0:
                s += wgt * float(np.sum(y * kern * ivar[k])) / den
        score[iz] = s
    return zs, score


def redshift_peaks(
    case: Case,
    n: int = 2,
    sep: float = 0.0020,
    n_cand: int = 5,
    refine: bool = True,
) -> list[float]:
    """Candidate redshifts: matched-filter peaks, then a local chi-square polish.

    ``sep`` keeps candidates at least 0.0020 apart in ``z`` (~600 km/s), which
    prevents two branches of the same peak from being carried forward as if they
    were independent starts.
    """
    zs, score = matched_filter_scan(case)
    order = np.argsort(score)[::-1]
    cands: list[float] = []
    for i in order:
        if not np.isfinite(score[i]):
            continue
        z = float(zs[i])
        if all(abs(z - q) > sep for q in cands):
            cands.append(z)
        if len(cands) >= n_cand:
            break
    if not refine:
        return cands[:n]
    polished: list[float] = []
    for z0 in cands:
        lo, hi = max(0.005, z0 - 0.004), z0 + 0.004
        zz, cc = linear_redshift_scan(case, lo, hi, 41, widths=(110.0, 250.0))
        if not np.all(np.isinf(cc)):
            z0 = float(zz[int(np.argmin(cc))])
        if all(abs(z0 - q) > sep for q in polished):
            polished.append(z0)
    if not polished:
        return cands[:n]
    return polished[:n]
