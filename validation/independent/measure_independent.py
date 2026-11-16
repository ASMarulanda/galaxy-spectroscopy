"""Independent emission-line measurement pipeline.

Scientific purpose
------------------
Provide a second, methodologically different route to the same 80-row output
contract, so that the benchmark's tolerances can be tested against something
other than the code that produced them.

Independence contract (see HANDOVER/15 §6)
------------------------------------------
* This module imports NO project code. Not ``common``, not ``methods``,
  not ``verify``. Only numpy / scipy / astropy.
* Nothing here calls ``fit_case``, ``cont_of``, ``model_all``, ``line_total``,
  ``flux_err``, ``clean``, ``matched_filter_scan``, ``linear_redshift_scan``
  or ``_bic``.

What is deliberately different from the reference implementation
----------------------------------------------------------------
==================  ===========================  ================================
stage               reference                    this implementation
==================  ===========================  ================================
continuum           polynomial per complex       penalised smoothing spline fitted
                                                 to the whole spectrum, iteratively
                                                 sigma-clipped
redshift            matched filter over 11       FFT cross-correlation in
                    lines + NNLS chi2 polish     log-wavelength space, parabolic
                                                 peak refinement
line model          one joint 25-parameter      separate per-complex fits;
                    bounded least squares        linear parameters (no log
                    over both complexes          parameterisation)
optimiser           scipy.optimize.least_squares scipy.optimize.minimize
                                                 (L-BFGS-B)
model selection     BIC over 7 nested models     nested F-test at p < 0.01
absorption shape    depth+width tied, ratio      fixed Balmer/G-band shape scaled
                    free in [1.15, 6]            by one depth parameter
flux                analytic A*sigma*sqrt(2pi)   numerical trapezoid of the fitted
                                                 component on a dense grid
uncertainty         covariance * min(chi2/dof,   residual bootstrap (parametric
                    1.5)                         noise), 40 draws
==================  ===========================  ================================

Known shared assumptions (stated, not hidden)
---------------------------------------------
The rest wavelengths, the [N II] 6583/6548 = 2.94 and [O III] 5007/4959 = 2.98
intensity ratios, and the fact that stellar absorption is broader than ionised
gas emission are properties of the physics, not of the reference code. Any
correct implementation must use them. They are written out here rather than
imported, so a change to the reference cannot silently propagate.

Author caveat
-------------
Both implementations were written by the same author. This is independence of
*code and method*, not of mind. The external catalogue comparison (step 6 of
the validation plan) is the only check that is independent of the author.
"""
from __future__ import annotations

import glob
import math
import os
import sys

import numpy as np
from astropy.io import fits
from scipy import stats
from scipy.interpolate import UnivariateSpline
from scipy.optimize import minimize

C_KMS = 299792.458
SQRT2PI = math.sqrt(2.0 * math.pi)

HALPHA = 6564.61
HBETA = 4862.68
OIII_5007 = 5008.24
OIII_4959 = 4960.29
NII_6583 = 6585.27
NII_6548 = 6549.86
SII_6716 = 6718.29
SII_6731 = 6732.67
HGAMMA = 4341.69
GBAND = 4305.00

R_NII = 2.94          # I(6583)/I(6548), atomic physics
R_OIII = 2.98         # I(5007)/I(4959), atomic physics

FATAL_BITS = (1 << 18) | (1 << 22) | (1 << 24)   # FULLREJECT, NOSKY, NODATA

N_BOOT = 40
BROAD_MAX = 1200.0      # km/s; above this a "broad line" is a pseudo-continuum
F_TEST_P = 0.01


# ---------------------------------------------------------------- input
def read_case(path):
    """Minimal FITS reader. Deliberately does not use the project's common.py."""
    with fits.open(path) as h:
        prim = h[0].header
        d = h[1].data
        out = dict(
            case_id=str(prim.get("CASEID", "")).strip() or os.path.basename(path)[:-5],
            flux=np.asarray(d["flux"], float),
            loglam=np.asarray(d["loglam"], float),
            ivar=np.asarray(d["ivar"], float),
            and_mask=np.asarray(d["and_mask"], np.int32),
            wdisp=np.asarray(d["wdisp"], float),
        )
    out["wave"] = 10.0 ** out["loglam"]
    out["good"] = (np.isfinite(out["flux"]) & np.isfinite(out["ivar"])
                   & (out["ivar"] > 0)
                   & ((out["and_mask"] & FATAL_BITS) == 0))
    return out


# ---------------------------------------------------------------- continuum
def fit_continuum(case, z, n_iter=3, clip=3.0, mask_halfwidth=0.004):
    """Penalised cubic smoothing spline fitted to LINE-FREE pixels only.

    The emission regions must be excluded before the fit: a cubic smoothing
    spline has enough freedom to follow a strong narrow line, and a continuum
    that passes through the Halpha peak destroys the measurement. Lines are
    masked using the CCF redshift, with a half-width of 0.4 % in wavelength so
    that broad components are excluded too.
    """
    wave, flux, ivar = case["wave"], case["flux"], case["ivar"]
    good = case["good"].copy()

    linemask = np.zeros_like(good)
    for lam in (HALPHA, NII_6583, NII_6548, SII_6716, SII_6731,
                HBETA, OIII_5007, OIII_4959, HGAMMA):
        mu = lam * (1.0 + z)
        linemask |= (wave > mu * (1 - mask_halfwidth)) & (wave < mu * (1 + mask_halfwidth))
    good &= ~linemask

    for _ in range(n_iter):
        idx = np.flatnonzero(good)
        if idx.size < 100:
            break
        x, y = wave[idx], flux[idx]
        w = np.sqrt(ivar[idx])
        # weighted residual budget; ~N_free for a model that fits to the noise
        s = float(idx.size) * 1.2
        try:
            spl = UnivariateSpline(x, y, w=w, s=s, k=3, ext=3)
        except Exception:                                   # noqa: BLE001
            return None
        resid = (flux - spl(wave)) * np.sqrt(ivar)
        ok = np.isfinite(resid) & (np.abs(np.where(np.isfinite(resid), resid, 0.0)) < clip)
        good = good & ok
    return spl


# ---------------------------------------------------------------- redshift
CCF_LINES = [(HALPHA, 1.00), (NII_6583, 0.35), (NII_6548, 0.12),
             (HBETA, 0.35), (OIII_5007, 0.32), (OIII_4959, 0.11),
             (SII_6716, 0.13), (SII_6731, 0.11), (HGAMMA, 0.16)]


def _running_median(y, w, width=201):
    n = len(y)
    half = width // 2
    out = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        out[i] = np.median(y[lo:hi])
    return out


def ccf_redshift(case, zlo=0.005, zhi=0.13, sig_lsf_pix=2.0):
    """FFT cross-correlation of the spectrum against a rest-frame line template.

    In log-wavelength space a redshift is a pure translation, so the search is a
    correlation rather than a template refit at every trial z.
    """
    loglam, flux, ivar, good = (case["loglam"], case["flux"],
                                case["ivar"], case["good"])
    dx = float(np.median(np.diff(loglam)))
    n = len(loglam)

    base = _running_median(flux, good)
    y = np.where(good, flux - base, 0.0) * np.sqrt(np.where(good, ivar, 0.0))

    # rest-frame template on the same uniform log-wavelength grid
    x = loglam
    sig_dex = sig_lsf_pix * dx          # the width is in PIXELS, not in dex
    t = np.zeros(n)
    for lam, wt in CCF_LINES:
        lx = math.log10(lam)
        t += wt * np.exp(-0.5 * ((x - lx) / sig_dex) ** 2)

    # Remove the mean of both series: without this the correlation is dominated
    # by the DC term and the peak lands in the low-overlap (wrap-around) region.
    y = y - y[good].mean()
    t = t - t.mean()

    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    cc = np.fft.irfft(np.fft.rfft(y, nfft) * np.conj(np.fft.rfft(t, nfft)), nfft)

    lags = np.arange(nfft)
    shift_dex = lags * dx                      # positive lag -> data is redshifted
    zgrid = 10.0 ** shift_dex - 1.0
    sel = (zgrid >= zlo) & (zgrid <= zhi)
    if not np.any(sel):
        return 0.05
    cand = np.flatnonzero(sel)
    k = cand[int(np.argmax(cc[cand]))]

    # parabolic refinement on the correlation peak
    if 0 < k < nfft - 1:
        y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        z = 10.0 ** ((k + delta) * dx) - 1.0
    else:
        z = float(zgrid[k])
    return float(z)


# ---------------------------------------------------------------- line model
class Complex:
    def __init__(self, lo, hi, comps, absorbers):
        self.lo, self.hi = lo, hi
        self.comps = comps            # list of dicts: name, rest, width_key, tie
        self.absorbers = absorbers    # list of (rest, weight, width_mult)


def ha_complex(allow_broad):
    comps = [
        dict(name="NII6548", rest=NII_6548, width="n", tie=("NII6583", 1.0 / R_NII)),
        dict(name="Ha_n", rest=HALPHA, width="n"),
        dict(name="NII6583", rest=NII_6583, width="n"),
        dict(name="SII6716", rest=SII_6716, width="n"),
        dict(name="SII6731", rest=SII_6731, width="n"),
    ]
    if allow_broad:
        comps.insert(2, dict(name="Ha_b", rest=HALPHA, width="b"))
    return Complex(6420.0, 6790.0, comps, [])


def hb_complex(allow_broad, absorption):
    comps = [
        dict(name="Hb_n", rest=HBETA, width="n"),
        dict(name="OIII4959", rest=OIII_4959, width="o", tie=("OIII5007", 1.0 / R_OIII)),
        dict(name="OIII5007", rest=OIII_5007, width="o"),
    ]
    if allow_broad:
        comps.insert(1, dict(name="Hb_b", rest=HBETA, width="b"))
    absorbers = ([(GBAND, 0.45, 3.0), (HGAMMA, 1.00, 1.0), (HBETA, 1.00, 1.0)]
                 if absorption else [])
    return Complex(4250.0, 5090.0, comps, absorbers)


def _gauss(w, mu, sig):
    return np.exp(-0.5 * ((w - mu) / sig) ** 2)


def _sig_pix(mu, sv_kms, wdisp, dlam):
    return np.sqrt((mu * sv_kms / C_KMS) ** 2 + (wdisp * dlam) ** 2)


def fit_complex(case, cx, z, cont_spline, broad_hint=None):
    """Fit one complex with a common narrow/broad width and free amplitudes.

    Parameters are linear (amplitudes in flux units, widths in km/s), which is a
    different parameterisation from the reference's log-widths.
    """
    wave, flux, ivar = case["wave"], case["flux"], case["ivar"]
    good = case["good"]
    wd, dl = case["wdisp"], np.gradient(case["wave"])

    m = good & (wave > cx.lo * (1 + z)) & (wave < cx.hi * (1 + z))
    if m.sum() < 60:
        return None
    w = wave[m]
    y = flux[m]
    base = cont_spline(w)
    dv = wd[m] * dl[m]
    iv = ivar[m]
    sw = np.sqrt(iv)

    free = [c for c in cx.comps if "tie" not in c]
    names = [c["name"] for c in free]
    namp = len(free)

    def unpack(p):
        sn, sb, so = p[0], p[1], p[2]
        off = 3
        depth = 0.0
        if cx.absorbers:
            depth = p[3]
            off = 4
        amps = {names[i]: p[off + i] for i in range(namp)}
        for c in cx.comps:
            if "tie" in c:
                src, fac = c["tie"]
                amps[c["name"]] = amps[src] * fac
        return sn, sb, so, depth, amps

    def model(p):
        sn, sb, so, depth, amps = unpack(p)
        cont = base.copy()
        for mu0, wt, wm in cx.absorbers:
            mu = mu0 * (1 + z)
            sig = wm * _sig_pix(mu, sn * 2.5, wd[m], dl[m])
            cont = cont * (1.0 - depth * wt * _gauss(w, mu, sig))
        tot = cont.copy()
        for c in cx.comps:
            sw_k = {"n": sn, "b": sb, "o": so}[c["width"]]
            mu = c["rest"] * (1 + z)
            sig = _sig_pix(mu, sw_k, wd[m], dl[m])
            tot = tot + amps[c["name"]] * _gauss(w, mu, sig)
        return tot

    def chi2(p):
        return float(np.sum(((y - model(p)) * sw) ** 2))

    # starting point: amplitudes from a narrow-band estimate, widths from the
    # observed line extent
    amax = max(float(np.max(np.abs(y - base))), 1e-3)
    p0 = [110.0, 600.0, 110.0] + ([0.10] if cx.absorbers else []) + [amax] * namp
    b_lo, b_hi = 250.0, BROAD_MAX
    if broad_hint is not None:
        # Broad Halpha and broad Hbeta come from the same broad-line region, so
        # their widths must agree. Passing the Halpha width as a tight bound is
        # both physically motivated and what makes the broad Hbeta identifiable
        # on spectra where it is weak.
        b_lo = max(250.0, 0.8 * broad_hint)
        b_hi = min(BROAD_MAX, 1.25 * broad_hint)
        p0[1] = min(max(broad_hint, b_lo), b_hi)
    bounds = ([(40.0, 400.0), (b_lo, b_hi), (40.0, 400.0)]
              + ([(0.0, 0.9)] if cx.absorbers else [])
              + [(0.0, None)] * namp)

    res = minimize(chi2, p0, method="L-BFGS-B", bounds=bounds,
                   options=dict(maxiter=400, ftol=1e-12, gtol=1e-10))
    if not res.success and not np.isfinite(res.fun):
        return None

    p = res.x
    ndata = int(m.sum())
    k_par = len(p)
    return dict(p=p, chi2=float(res.fun), ndata=ndata, k=k_par,
                dof=max(ndata - k_par, 1), z=z, mask=m, names=names,
                free=free, unpack=unpack, model=model, wave=w, flux=y,
                cont=base, ivar=iv, cx=cx)


def line_key(name):
    return {"Ha_n": "Halpha", "Ha_b": "Halpha", "Hb_n": "Hbeta",
            "Hb_b": "Hbeta", "OIII5007": "OIII_5007",
            "NII6583": "NII_6583"}.get(name)


def flux_of_line(fit, line, case):
    """Total integrated flux of a line: sum over its kinematic components.

    The integral is computed numerically on a dense grid, including the
    instrumental width measured from the file, so this does not rely on the
    analytic Gaussian formula anywhere in the project.
    """
    cx = fit["cx"]
    sn, sb, so, depth, amps = fit["unpack"](fit["p"])
    total = 0.0
    for c in cx.comps:
        if line_key(c["name"]) != line or "tie" in c:
            continue
        sv = {"n": sn, "b": sb, "o": so}[c["width"]]
        mu = c["rest"] * (1 + fit["z"])
        sig_int = mu * sv / C_KMS
        lsf = local_lsf_sigma(case, mu)
        sig = math.sqrt(sig_int ** 2 + lsf ** 2)
        grid = np.linspace(mu - 8 * sig, mu + 8 * sig, 2001)
        total += float(np.trapz(
            amps[c["name"]] * np.exp(-0.5 * ((grid - mu) / sig) ** 2), grid))
    return total


def local_lsf_sigma(case, mu):
    wave, wdisp = case["wave"], case["wdisp"]
    j = int(np.argmin(np.abs(wave - mu)))
    dlam = wave[j] * 2.3025851 * np.median(np.diff(case["loglam"]))
    return float(wdisp[j] * dlam)


# ---------------------------------------------------------------- model choice
def broad_is_degenerate(fit):
    """True when the broad component has run to its upper bound.

    Over a 340 A window a Gaussian that broad is indistinguishable from a
    slowly varying continuum, and it absorbs neighbouring line flux. Phase 2
    measured exactly this failure on the broad-line AGN R05.
    """
    try:
        _, sb, _, _, _ = fit["unpack"](fit["p"])
    except Exception:                                       # noqa: BLE001
        return True
    return bool(sb > 0.95 * BROAD_MAX)


def f_test(red, full):
    """F-test for nesting; returns True when the extra parameters are justified."""
    if red is None or full is None:
        return False
    dk = full["k"] - red["k"]
    if dk <= 0:
        return False
    denom_df = max(full["ndata"] - full["k"], 1)
    f = ((red["chi2"] - full["chi2"]) / dk) / (full["chi2"] / denom_df)
    if f <= 0:
        return False
    p = stats.f.sf(f, dk, denom_df)
    return bool(p < F_TEST_P)


# ---------------------------------------------------------------- bootstrap
def bootstrap_sigma(case, fit, line, n_boot=N_BOOT, seed=0):
    """Parametric residual bootstrap over the fitted parameters."""
    rng = np.random.default_rng(seed)
    base_flux = flux_of_line(fit, line, case)
    w = fit["wave"]
    ymod = fit["model"](fit["p"])
    sig = 1.0 / np.sqrt(fit["ivar"])

    vals = []
    for _ in range(n_boot):
        yb = ymod + rng.normal(0.0, 1.0, size=len(w)) * sig
        # refit with the same parameterisation on the resampled data
        sub = dict(case)
        sub["flux"] = case["flux"].copy()
        sub["flux"][fit["mask"]] = yb
        f2 = fit_complex(sub, fit["cx"], fit["z"], fit["_spline"])
        if f2 is None:
            continue
        vals.append(flux_of_line(f2, line, case))
    if len(vals) < 5:
        return base_flux, float("nan")
    return float(np.mean(vals)), float(np.std(vals, ddof=1))


# ---------------------------------------------------------------- driver
def measure_case(case):
    # redshift first: the continuum fit needs to know where the lines are
    z = ccf_redshift(case)
    spline = fit_continuum(case, z)
    if spline is None:
        return None

    # Refine z on the Halpha complex. A coarse bracket is followed by a fine
    # scan whose chi2(z) also yields the 1-sigma redshift uncertainty from the
    # Delta-chi2 = 1 interval.
    z_err = float("nan")
    for half, nstep in ((0.0015, 15), (0.0006, 61)):
        zs, c2 = [], []
        for dz in np.linspace(-half, half, nstep):
            f = fit_complex(case, ha_complex(True), z + dz, spline)
            if f is None:
                continue
            zs.append(z + dz)
            c2.append(f["chi2"])
        if not zs:
            continue
        zs, c2 = np.asarray(zs), np.asarray(c2)
        k = int(np.argmin(c2))
        z = float(zs[k])
        inside = zs[np.abs(c2 - c2[k]) <= 1.0]
        if inside.size > 1:
            z_err = float(0.5 * (inside.max() - inside.min()))
    if not np.isfinite(z_err) or z_err <= 0:
        z_err = 0.0006 / 61      # conservative: half a fine-scan step

    out = dict(z=z, z_err=z_err, lines={})

    # ---- Halpha complex: narrow-only vs narrow+broad
    cands = []
    for allow_broad in (False, True):
        f = fit_complex(case, ha_complex(allow_broad), z, spline)
        if f is not None:
            f["_spline"] = spline
            cands.append(f)
    ha_fit = cands[0] if cands else None
    if len(cands) == 2 and f_test(cands[0], cands[1]) and not broad_is_degenerate(cands[1]):
        ha_fit = cands[1]

    # ---- Hbeta complex: with / without absorption, narrow-only vs +broad
    broad_hint = None
    if ha_fit is not None:
        _, sb_ha, _, _, _ = ha_fit["unpack"](ha_fit["p"])
        if 250.0 < sb_ha < 0.95 * BROAD_MAX:
            broad_hint = float(sb_ha)

    hb_cands = {}
    for with_abs in (True, False):
        for allow_broad in (False, True):
            f = fit_complex(case, hb_complex(allow_broad, with_abs), z, spline,
                            broad_hint=broad_hint if allow_broad else None)
            if f is not None:
                f["_spline"] = spline
                hb_cands[(with_abs, allow_broad)] = f
    hb_fit = None
    if hb_cands:
        base_abs = hb_cands.get((True, False)) or hb_cands.get((False, False))
        for key in ((True, True), (False, True)):
            if (key in hb_cands and not broad_is_degenerate(hb_cands[key])
                    and f_test(base_abs, hb_cands[key])):
                base_abs = hb_cands[key]
        hb_fit = base_abs

    for fit, lines in ((ha_fit, ("Halpha", "NII_6583")),
                       (hb_fit, ("Hbeta", "OIII_5007"))):
        if fit is None:
            continue
        for line in lines:
            flux = flux_of_line(fit, line, case)
            _, ferr = bootstrap_sigma(case, fit, line)
            if not np.isfinite(flux) or not np.isfinite(ferr) or ferr <= 0:
                flag = "failed"
            elif abs(flux) / ferr < 3.0:
                flag = "nondetection"
            else:
                flag = "ok"
            out["lines"][line] = dict(flux=flux, flux_err=ferr, flag=flag)
    return out


def main(case_dir, out_csv):
    import csv
    files = sorted(glob.glob(os.path.join(case_dir, "*.fits")))
    rows = []
    for fn in files:
        case = read_case(fn)
        cid = case["case_id"]
        try:
            r = measure_case(case)
        except Exception as exc:                            # noqa: BLE001
            print(f"  {cid}: FAILED {exc}")
            r = None
        if r is None:
            r = dict(z=float("nan"), lines={})
        for line in ("Halpha", "Hbeta", "OIII_5007", "NII_6583"):
            rec = r["lines"].get(line)
            if rec is None:
                rows.append(dict(case_id=cid, line=line, redshift=float("nan"),
                                 redshift_err=float("nan"), flux=float("nan"),
                                 flux_err=float("nan"), flag="failed"))
            else:
                rows.append(dict(case_id=cid, line=line, redshift=r["z"],
                                 redshift_err=r.get("z_err", float("nan")),
                                 flux=rec["flux"], flux_err=rec["flux_err"],
                                 flag=rec["flag"]))
        ha = r["lines"].get("Halpha", {})
        print(f"  {cid:5} z={r['z']:.5f} "
              f"Ha={ha.get('flux', float('nan')):8.2f} "
              f"+-{ha.get('flux_err', float('nan')):7.2f} {ha.get('flag','?')}",
              flush=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case_id", "line", "redshift",
                                           "redshift_err", "flux", "flux_err", "flag"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    print(f"wrote {out_csv}: {len(rows)} rows")


if __name__ == "__main__":
    base = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    main(os.path.join(base, "data", "cases"), os.path.join(os.path.dirname(__file__), "out_ind.csv"))
