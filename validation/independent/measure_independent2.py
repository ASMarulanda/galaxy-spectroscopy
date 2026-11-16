"""Independent pipeline v2: JOINT continuum and stellar absorption.

Why v2 exists
-------------
v1 fitted the continuum once (a smoothing spline on line-free pixels) and then
held it fixed while fitting the lines. On the broad-line cases that is exactly
the wrong decomposition: the spline is free enough to follow part of the strong
stellar Hbeta absorption, so the data no longer support a separate broad Hbeta
component and the fit pushes the residual into the narrow line (S05 -39 %,
S07 -48 %; see HANDOVER/16 §4.4).

v2 removes the two-stage structure. Continuum and absorption are fitted
SIMULTANEOUSLY with the emission lines inside each complex window, so the fit
must choose where the flux belongs instead of inheriting that choice from a
previous step. Two consequences follow, and both matter:

* the continuum is a free-form radial-basis function, not a polynomial, so the
  estimator stays methodologically distinct from the reference;
* the bootstrap can resample the data and refit the whole window model, which
  propagates continuum uncertainty automatically rather than only line-window
  noise.

Independence
------------
No project module is imported. This file imports only numpy/scipy, plus the FITS
reader and CCF redshift from ``measure_independent`` (v1), which are themselves
independent of the reference solution. It never calls any function from
``methods.py``.

Model, per complex window
-------------------------
    f(lambda) = [ sum_k c_k B_k(lambda) ] * prod_i [1 - d w_i G_i(lambda) ]
                + sum_j A_j G_j(lambda)

B_k    radial basis functions spanning the window (continuum, free)
G_i    stellar absorption lines, width = sigma_star * wmult_i (d, sigma_star free)
G_j    emission components, widths from one narrow and one broad velocity
       dispersion shared across the complex
"""
from __future__ import annotations

import glob
import math
import os
import sys

import numpy as np
from scipy import stats
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from measure_independent import (  # noqa: E402  (v1: I/O and CCF only)
    C_KMS, HALPHA, HBETA, NII_6583, NII_6548, OIII_5007, OIII_4959,
    SII_6716, SII_6731, HGAMMA, GBAND, R_NII, R_OIII,
    read_case, ccf_redshift, local_lsf_sigma, line_key,
)

BROAD_MAX = 1200.0
N_RBF = 6           # continuum basis functions per window
F_TEST_P = 0.01


# ------------------------------------------------------------------ model
def rbf_basis(w, n=N_RBF):
    """Partition-of-unity style Gaussian basis across the window."""
    lo, hi = float(w.min()), float(w.max())
    centres = np.linspace(lo, hi, n)
    width = (hi - lo) / max(n - 1, 1)
    return np.exp(-0.5 * ((w[:, None] - centres[None, :]) / width) ** 2)


class Complex:
    def __init__(self, lo, hi, comps, absorbers):
        self.lo, self.hi = lo, hi
        self.comps = comps
        self.absorbers = absorbers


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


def hb_complex(allow_broad, absorption=True):
    comps = [
        # Hgamma emission sits inside this window at 15-35 % of Hbeta and is
        # superimposed on the Hgamma absorption line; leaving it out biases the
        # local continuum.
        dict(name="Hg_em", rest=HGAMMA, width="n"),
        dict(name="Hb_n", rest=HBETA, width="n"),
        dict(name="OIII4959", rest=OIII_4959, width="o", tie=("OIII5007", 1.0 / R_OIII)),
        dict(name="OIII5007", rest=OIII_5007, width="o"),
    ]
    if allow_broad:
        comps.insert(1, dict(name="Hb_b", rest=HBETA, width="b"))
    absorbers = ([(GBAND, 0.45, 3.0), (HGAMMA, 1.00, 1.0), (HBETA, 1.00, 1.0)]
                 if absorption else [])
    return Complex(4250.0, 5090.0, comps, absorbers)


def _sig_pix(mu, sv, wd, dl):
    return np.sqrt((mu * sv / C_KMS) ** 2 + (wd * dl) ** 2)


def fit_window(case, cx, z, cont_init, broad_hint=None, maxiter=600, n_rbf=None):
    """Joint fit of continuum, stellar absorption and emission in one window.

    cont_init : initial continuum values on the window pixels (from the spline)
    """
    wave, flux, ivar = case["wave"], case["flux"], case["ivar"]
    good = case["good"]
    wd_all, dl_all = case["wdisp"], np.gradient(case["wave"])

    m = good & (wave > cx.lo * (1 + z)) & (wave < cx.hi * (1 + z))
    if m.sum() < 80:
        return None
    w = wave[m]
    y = flux[m]
    iv = ivar[m]
    cont_w = np.asarray(cont_init, float)[m] if np.ndim(cont_init) else float(cont_init)
    sw = np.sqrt(iv)
    wd, dl = wd_all[m], dl_all[m]

    NB = int(n_rbf or N_RBF)
    B = rbf_basis(w, NB)                  # (n, NB)
    free = [c for c in cx.comps if "tie" not in c]
    names = [c["name"] for c in free]
    namp = len(names)

    # ---- parameter layout
    # [0:N_RBF] continuum coefficients (log-scaled to stay positive)
    # [N_RBF] absorption depth
    # [N_RBF+1] sigma_star (km/s)
    # [N_RBF+2] narrow sigma (km/s)
    # [N_RBF+3] broad sigma (km/s)      (unused when no broad component)
    # then amplitudes
    P = NB
    i_depth, i_sstar, i_narrow, i_broad = P, P + 1, P + 2, P + 3
    i_amp = P + 4
    npar = i_amp + namp

    has_abs = bool(cx.absorbers)
    has_broad = any(c["width"] == "b" for c in cx.comps)

    def unpack(p):
        cont = B @ np.exp(p[:NB])
        depth = p[i_depth] if has_abs else 0.0
        sstar = p[i_sstar] if has_abs else 200.0
        sn = p[i_narrow]
        sb = p[i_broad] if has_broad else 600.0
        amps = {names[i]: p[i_amp + i] for i in range(namp)}
        for c in cx.comps:
            if "tie" in c:
                src, fac = c["tie"]
                amps[c["name"]] = amps[src] * fac
        return cont, depth, sstar, sn, sb, amps

    def atten_of(depth, sstar):
        a = np.ones_like(w)
        for mu0, wt, wm in cx.absorbers:
            mu = mu0 * (1 + z)
            sig = wm * _sig_pix(mu, sstar, wd, dl)
            a = a * (1.0 - depth * wt * np.exp(-0.5 * ((w - mu) / sig) ** 2))
        return a

    def model(p):
        cont, depth, sstar, sn, sb, amps = unpack(p)
        tot = cont * atten_of(depth, sstar)
        for c in cx.comps:
            sv = {"n": sn, "b": sb, "o": sn}[c["width"]]
            mu = c["rest"] * (1 + z)
            sig = _sig_pix(mu, sv, wd, dl)
            tot = tot + amps[c["name"]] * np.exp(-0.5 * ((w - mu) / sig) ** 2)
        return tot

    def chi2(p):
        r = (y - model(p)) * sw
        return float(np.sum(r * r))

    # ---- initial values
    coef0 = np.linalg.lstsq(B, np.maximum(cont_w, 1e-3), rcond=None)[0]
    p0 = np.zeros(npar)
    p0[:NB] = np.log(np.maximum(coef0, 1e-3))
    p0[i_depth] = 0.10
    p0[i_sstar] = 180.0
    p0[i_narrow] = 110.0
    p0[i_broad] = broad_hint if broad_hint else 600.0
    amax = max(float(np.max(np.abs(y - cont_w))), 1e-3)
    p0[i_amp:i_amp + namp] = amax

    lo = np.full(npar, -np.inf)
    hi = np.full(npar, np.inf)
    lo[:NB], hi[:NB] = np.log(1e-4), np.log(1e4)
    lo[i_depth], hi[i_depth] = 0.0, 0.60
    lo[i_sstar], hi[i_sstar] = 100.0, 400.0
    lo[i_narrow], hi[i_narrow] = 40.0, 400.0
    if has_broad:
        b_lo, b_hi = 250.0, BROAD_MAX
        if broad_hint:
            b_lo = max(250.0, 0.8 * broad_hint)
            b_hi = min(BROAD_MAX, 1.25 * broad_hint)
            p0[i_broad] = min(max(broad_hint, b_lo), b_hi)
        lo[i_broad], hi[i_broad] = b_lo, b_hi
    lo[i_amp:], hi[i_amp:] = 0.0, np.inf

    bounds = list(zip(lo.tolist(), hi.tolist()))
    res = minimize(chi2, p0, method="L-BFGS-B", bounds=bounds,
                   options=dict(maxiter=maxiter, ftol=1e-11, gtol=1e-9))
    if not np.isfinite(res.fun):
        return None

    return dict(p=res.x, chi2=float(res.fun), npar=npar, ndata=int(m.sum()),
                dof=max(int(m.sum()) - npar, 1), z=z, cx=cx,
                unpack=unpack, model=model, wave=w, mask=m, has_broad=has_broad)


# ------------------------------------------------------------------ flux
def flux_of_line(fit, line, case):
    cont, depth, sstar, sn, sb, amps = fit["unpack"](fit["p"])
    total = 0.0
    for c in fit["cx"].comps:
        if line_key(c["name"]) != line or "tie" in c:
            continue
        sv = {"n": sn, "b": sb, "o": sn}[c["width"]]
        mu = c["rest"] * (1 + fit["z"])
        sig_int = mu * sv / C_KMS
        lsf = local_lsf_sigma(case, mu)
        sig = math.sqrt(sig_int ** 2 + lsf ** 2)
        grid = np.linspace(mu - 8 * sig, mu + 8 * sig, 2001)
        total += float(np.trapz(
            amps[c["name"]] * np.exp(-0.5 * ((grid - mu) / sig) ** 2), grid))
    return total


def broad_degenerate(fit):
    if not fit.get("has_broad"):
        return False
    _, _, _, _, sb, _ = fit["unpack"](fit["p"])
    return bool(sb > 0.95 * BROAD_MAX)


def f_test(red, full):
    if red is None or full is None:
        return False
    dk = full["npar"] - red["npar"]
    if dk <= 0:
        return False
    dof = max(full["ndata"] - full["npar"], 1)
    f = ((red["chi2"] - full["chi2"]) / dk) / (full["chi2"] / dof)
    if f <= 0:
        return False
    return bool(stats.f.sf(f, dk, dof) < F_TEST_P)


# ------------------------------------------------------------------ driver
def measure_case(case, n_boot=25, seed=0, verbose=False):
    from measure_independent import fit_continuum

    z = ccf_redshift(case)
    spl = fit_continuum(case, z)
    if spl is None:
        return None

    # refine z on the Halpha window, then the Delta-chi2 = 1 interval
    ha_models = {}
    for allow_broad in (False, True):
        f = fit_window(case, ha_complex(allow_broad), z, spl(case["wave"]))
        if f is not None:
            ha_models[allow_broad] = f
    if not ha_models:
        return None

    def ha_select(zz):
        models = {}
        for ab in (False, True):
            f = fit_window(case, ha_complex(ab), zz, spl(case["wave"]))
            if f is not None:
                models[ab] = f
        if not models:
            return None
        best = models.get(False)
        if True in models and f_test(models.get(False), models[True]) \
                and not broad_degenerate(models[True]):
            best = models[True]
        return best

    z_err = float("nan")
    for half, nstep in ((0.0015, 11), (0.0006, 21)):
        zs, c2 = [], []
        for dz in np.linspace(-half, half, nstep):
            f = ha_select(z + dz)
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
        z_err = 1e-5

    cont = spl(case["wave"])
    ha_fit = ha_select(z)

    broad_hint = None
    if ha_fit is not None and ha_fit["has_broad"]:
        _, _, _, _, sb, _ = ha_fit["unpack"](ha_fit["p"])
        if 250 < sb < 0.95 * BROAD_MAX:
            broad_hint = float(sb)

    hb_models = {}
    for with_abs in (True, False):
        for ab in (False, True):
            f = fit_window(case, hb_complex(ab, with_abs), z, cont,
                           broad_hint=broad_hint if ab else None)
            if f is not None:
                hb_models[(with_abs, ab)] = f
    hb_fit = None
    if hb_models:
        base = None
        for key in ((True, False), (False, False)):
            if key in hb_models:
                base = hb_models[key]
                break
        for key in ((True, True), (False, True)):
            if key in hb_models and f_test(base, hb_models[key]) \
                    and not broad_degenerate(hb_models[key]):
                base = hb_models[key]
        hb_fit = base

    out = dict(z=z, z_err=z_err, lines={})

    for fit, lines in ((ha_fit, ("Halpha", "NII_6583")),
                       (hb_fit, ("Hbeta", "OIII_5007"))):
        if fit is None:
            continue
        for line in lines:
            flux = flux_of_line(fit, line, case)
            ferr = bootstrap_sigma(case, fit, line, z, spl, n_boot, seed)
            if not np.isfinite(flux) or not np.isfinite(ferr) or ferr <= 0:
                flag = "failed"
            elif abs(flux) / ferr < 3.0:
                flag = "nondetection"
            else:
                flag = "ok"
            out["lines"][line] = dict(flux=flux, flux_err=ferr, flag=flag)
    return out


def bootstrap_sigma(case, fit, line, z, spl, n_boot=25, seed=0):
    """Resample the whole spectrum and refit the JOINT window model.

    Because continuum, absorption and lines are fitted together, perturbing the
    data and refitting propagates continuum and absorption uncertainty as well
    as line-window noise.
    """
    rng = np.random.default_rng(seed + 17)
    sig = 1.0 / np.sqrt(np.where(case["ivar"] > 0, case["ivar"], np.inf))
    vals = [flux_of_line(fit, line, case)]
    for _ in range(n_boot):
        sub = dict(case)
        sub["flux"] = case["flux"] + rng.normal(0.0, 1.0, case["flux"].size) * sig
        f2 = fit_window(sub, fit["cx"], z, spl(case["wave"]),
                        broad_hint=None, maxiter=150)
        if f2 is None:
            continue
        vals.append(flux_of_line(f2, line, sub))
    if len(vals) < 5:
        return float("nan")
    return float(np.std(vals, ddof=1))


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
            print(f"  {cid}: FAILED {exc}", flush=True)
            r = None
        if r is None:
            r = dict(z=float("nan"), z_err=float("nan"), lines={})
        for line in ("Halpha", "Hbeta", "OIII_5007", "NII_6583"):
            rec = r["lines"].get(line)
            if rec is None:
                rows.append(dict(case_id=cid, line=line, redshift=float("nan"),
                                 redshift_err=float("nan"), flux=float("nan"),
                                 flux_err=float("nan"), flag="failed"))
            else:
                rows.append(dict(case_id=cid, line=line, redshift=r["z"],
                                 redshift_err=r["z_err"],
                                 flux=rec["flux"], flux_err=rec["flux_err"],
                                 flag=rec["flag"]))
        hb = r["lines"].get("Hbeta", {})
        print(f"  {cid:5} z={r['z']:.5f} Hb={hb.get('flux', float('nan')):8.2f} "
              f"+-{hb.get('flux_err', float('nan')):6.2f} {hb.get('flag', '?')}",
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
    main(os.path.join(base, "data", "cases"),
         os.path.join(os.path.dirname(__file__), "out_ind2.csv"))
