#!/usr/bin/env python3
"""A guided walk through the pipeline on one real case (S01).

Prints every intermediate value the documentation talks about, in the order the
pipeline computes them. Nothing here is a reimplementation: each step calls the
public module that does the work, and names it.

    python tools/walkthrough_s01.py            # default: data/cases/S01.fits
    python tools/walkthrough_s01.py S03        # any case id

Read alongside the v0.2 methods document; the section numbers below match it.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

import gspec
from gspec import constants as K
from gspec.covariance import covariance, covariance_scale, sigma_of_flux
from gspec.fit import build_complexes
from gspec.fluxes import line_total
from gspec.io import clean
from gspec.lines import LINE_OF_KEY
from gspec.pipeline import measure
from gspec.redshift import matched_filter_scan, redshift_peaks

ROOT = Path(__file__).resolve().parents[1]
CASE = sys.argv[1] if len(sys.argv) > 1 else "S01"
PATH = ROOT / "data" / "cases" / f"{CASE}.fits"


def rule(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def main() -> int:
    from astropy.io import fits

    # ---------------------------------------------------------------- step 1
    rule(f"STEP 1 — loading the FITS file  ({PATH.name})")
    print("module: gspec.io.read_case()   (methods.md §1)")
    with fits.open(PATH) as hdul:
        hdul.info()
        prim = hdul[0].header
        cols = [(c.name, c.format, c.unit) for c in hdul[1].columns]
    for key in ("CASEID", "BUNIT", "VACUUM", "PROVEN", "NWAVE"):
        if key in prim:
            print(f"  PRIMARY {key:7} = {prim[key]!r}")
    print("  COADD columns:")
    for name, fmt, unit in cols:
        print(f"    {name:9} {fmt:3} {unit or ''}")

    case = gspec.read_case(PATH)

    # ---------------------------------------------------------------- step 2
    rule("STEP 2 — wavelength, flux density and inverse variance  (methods.md §1-2)")
    print("module: gspec.io.Case (wave = 10**loglam, dlam from the native grid)")
    dlog = float(np.median(np.diff(case.loglam)))
    dlam_ha = float(case.dlam[np.argmin(np.abs(case.wave - K.HALPHA * (1 + 0.04)))])
    good = clean(case)
    print(f"  pixels            {len(case)}")
    print(f"  wavelength        {case.wave[0]:.1f} – {case.wave[-1]:.1f} Angstrom (vacuum)")
    print(f"  dlog10(lambda)    {dlog:.4e} dex  -> dlambda at Halpha = {dlam_ha:.4f} A/px")
    print(f"  flux density      median {np.median(case.flux):.3f} "
          f"(unit: {case.bunit})   <- a DENSITY, not an integral")
    print(f"  inverse variance  median sigma = {np.median(1/np.sqrt(case.ivar[good])):.4f}")
    print(f"  usable pixels     {good.sum()} / {len(case)}  ({100*(1-good.mean()):.2f} % rejected)")
    print(f"  instrumental LSF  wdisp at Halpha = "
          f"{case.wdisp[np.argmin(np.abs(case.wave - K.HALPHA*(1+0.04)))]:.3f} px")

    # ---------------------------------------------------------------- step 3
    rule("STEP 3 — redshift estimation  (methods.md §2)")
    print("modules: gspec.redshift.matched_filter_scan() -> redshift_peaks()")
    zs, score = matched_filter_scan(case)
    top = np.argsort(score)[::-1][:5]
    print("  matched-filter peaks (trial grid 0.008-0.102, 1400 points):")
    for i in top:
        print(f"    z = {zs[i]:.5f}   score = {score[i]:8.2f}")
    cands = redshift_peaks(case, n=2)
    print(f"  polished candidates carried into the fit: "
          f"{[round(c, 6) for c in cands]}")
    result = measure(case)
    print(f"  reported redshift  z = {result.redshift:.8f}  +- {result.redshift_err:.2e}")
    truth = json.loads((ROOT / "data" / "truth" / "truth_synthetic.json").read_text())
    if CASE in truth:
        dz = (result.redshift - truth[CASE]["z"]) * K.C_KMS / (1 + truth[CASE]["z"])
        print(f"  analytic truth     z = {truth[CASE]['z']:.8f}  ->  error {dz:+.2f} km/s")

    # ---------------------------------------------------------------- step 4
    rule("STEP 4 — continuum and emission lines  (methods.md §3-4)")
    print("modules: gspec.select.fit_with_model_selection() -> fit_case() -> "
          "model.build_segments()")
    print(f"  model selected     {result.model_selected}")
    print(f"  parameters         {result.n_par}   dof {result.dof}   "
          f"chi2/dof {result.chi2dof:.4f}")
    print(f"  covariance scale   min(chi2/dof, {K.COV_SCALE_CAP}) = "
          f"{result.covariance_scale:.4f}")
    print(f"  pinned at a bound  {result.at_bounds or 'none'}")
    complexes = build_complexes(result.model)
    for cx in complexes:
        keys = [e.key for e in cx.emis]
        print(f"  {cx.name} complex: {cx.lo_rest:.0f}-{cx.hi_rest:.0f} A rest, "
              f"components {keys}")
    print("  free amplitudes (the ties are resolved, not fitted):")
    for k in ("Ha_n", "Ha_b", "NII6583", "Hg_em", "Hb_n", "Hb_b", "OIII5007"):
        if f"amp_{k}" in result.parameters:
            print(f"    amp_{k:9} = {result.parameters[f'amp_{k}']:10.4f}"
                  f"   -> line {LINE_OF_KEY.get(k)}")

    # ---------------------------------------------------------------- step 5
    rule("STEP 5 — integrated line fluxes  (methods.md §5)")
    print("module: gspec.fluxes.components() / integrated_flux()")
    print("  F = A * sqrt(sigma_int^2 + sigma_LSF^2) * sqrt(2 pi), summed per line")
    for line in ("Halpha", "Hbeta", "OIII_5007", "NII_6583"):
        print(f"  {line}:")
        for c in result.lines[line].components:
            print(f"    {c.key:10} A={c.amplitude:9.4f}  sigma_int={c.sigma_int:7.4f} A  "
                  f"sigma_LSF={c.sigma_lsf:6.4f} A  sigma_obs={c.sigma_obs:7.4f} A  "
                  f"F={c.flux:9.3f}")
        print(f"    {'TOTAL':10} {result.lines[line].flux:44.3f}"
              f"   (unit {K.UNITS})")
        if CASE in truth:
            t = truth[CASE]["lines"][line]["flux"]
            print(f"    {'truth':10} {t:44.3f}   error "
                  f"{100*(result.lines[line].flux/t-1):+.2f} %")

    # ---------------------------------------------------------------- step 6
    rule("STEP 6 — uncertainties  (methods.md §6)")
    print("modules: gspec.covariance.covariance() / sigma_of_flux()")
    print("  Cov = pinv(J'J) * min(chi2/dof, cap);  sigma_F^2 = grad' Cov grad")
    print(f"  {'line':11} {'F':>10} {'sigma_F':>9} {'sigma/F':>9} {'|F|/sigma':>10}  flag")
    for line, m in result.lines.items():
        print(f"  {line:11} {m.flux:10.3f} {m.flux_err:9.4f} "
              f"{m.flux_err/m.flux:9.2%} {abs(m.flux)/m.flux_err:10.2f}  {m.flag}")
    print("  the gradient is a forward difference of the same function that was")
    print("  integrated, so a change in the flux formula cannot leave sigma behind.")

    # ---------------------------------------------------------------- step 7
    rule("STEP 7 — validation against analytic truth  (the v0.2 validation report §2)")
    if CASE not in truth:
        print(f"  {CASE} is a real spectrum: no analytic truth. It is validated on")
        print("  redshift and the archival band instead (the v0.2 validation report §5).")
        return 0
    band = {"Halpha": 0.07, "Hbeta": 0.20, "OIII_5007": 0.06, "NII_6583": 0.10}
    print(f"  {'line':11} {'measured':>10} {'truth':>10} {'error':>9} {'band':>7} "
          f"{'3 sigma':>9}  verdict")
    for line, m in result.lines.items():
        t = truth[CASE]["lines"][line]["flux"]
        err = m.flux - t
        tol = max(band[line] * t, 3.0 * m.flux_err)
        verdict = "inside band" if abs(err) <= band[line] * t else (
            "inside via 3 sigma" if abs(err) <= tol else "OUTSIDE")
        print(f"  {line:11} {m.flux:10.3f} {t:10.3f} {100*err/t:+8.2f}% "
              f"{100*band[line]:6.0f}% {3*m.flux_err:9.3f}  {verdict}")
    print("\n  A row inside the band is an accurate measurement; a row inside only")
    print("  via 3 sigma is an honest non-measurement. The distinction is the")
    print("  whole point of reporting uncertainties, and it is why this repository")
    print("  quotes medians and not pass counts.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
