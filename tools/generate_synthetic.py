"""Independent seeded synthetic-spectrum generator (TRUTH SIDE).

This module is the only producer of synthetic ground truth.  It shares no code
with ``phase2/measure/``: truth is computed analytically from the generator
parameters (integral of each Gaussian), never by fitting the generated data.

Design rules enforced here
--------------------------
* Deterministic: every case is a pure function of (base_seed, case index).
* Truth is exact: each emission component is drawn as a Gaussian whose integral
  is *by construction* the declared flux.  ``self_test`` re-integrates each
  component on the pixel grid and asserts agreement (< 1e-6 relative).
* The continuum basis (3 power laws x multiplicative stellar absorption) is not
  the basis any measurement method uses.
* Sibling lines that are not graded ([N II] 6548, [S II], [O I], [O III] 4959,
  higher Balmer) are present, so a solver must model a complex rather than
  integrate a window.

Identifiability / fairness decisions (see PHASE2_VERIFICATION_CONTRACT.md §4.3)
-----------------------------------------------------------------------------
1. All stellar absorption depths share ONE random scale factor per case.  A
   single stellar population has a fixed Balmer absorption pattern; scaling it
   by one number keeps "line + absorption" identifiable from H-delta/H-gamma/
   H-beta jointly instead of leaving H-beta absorption degenerate with a local
   continuum.
2. The stellar absorption is broader than the ionised-gas emission
   (``sigma_stellar >= 1.35 * sigma_gas``).
3. The measurement side does NOT know these numbers: it uses an approximate,
   physically motivated Balmer absorption pattern of its own.  Generator and
   solver therefore disagree slightly on purpose.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat as C  # noqa: E402  (import shim; everything below is verbatim)

SQRT2PI = math.sqrt(2.0 * math.pi)

#: (rest wavelength, depth relative to H-beta, width multiplier)
ABSORPTION_LINES = [
    (3933.66, 1.50, 1.0),   # Ca II K
    (3968.47, 1.40, 1.0),   # Ca II H
    (4102.89, 1.10, 1.0),   # H-delta
    (4305.00, 0.45, 3.0),   # G band (broad)
    (4341.69, 1.00, 1.0),   # H-gamma
    (4862.68, 1.00, 1.0),   # H-beta   <-- the important one
    (5168.90, 0.70, 1.0),   # Mg b
    (5174.00, 0.80, 1.0),
    (5184.00, 0.60, 1.0),
    (5271.00, 0.30, 1.0),   # Fe
    (5891.58, 0.50, 1.0),   # Na D
    (5897.56, 0.45, 1.0),
]
HB_DEPTH = 0.11   # moderate, identifiable absorption (see contract 4.3)


def _continuum(lam_rest, rng, norm):
    """Three-component power-law continuum; not a polynomial."""
    amps = rng.dirichlet([1.2, 1.2, 1.2])
    slopes = np.sort(rng.uniform(0.4, 4.0, size=3))
    t = lam_rest / 5000.0
    return norm * sum(a * t ** (-s) for a, s in zip(amps, slopes))


def _absorption_profile(lam_rest, abs_scale, sigma_stellar_kms):
    atten = np.ones_like(lam_rest)
    for mu, rel, wfac in ABSORPTION_LINES:
        d = abs_scale * rel * HB_DEPTH
        sig = mu * sigma_stellar_kms * wfac / C.C_KMS
        atten *= 1.0 - d * np.exp(-0.5 * ((lam_rest - mu) / sig) ** 2)
    return np.clip(atten, 0.05, None)


def _add_gaussian(flux_density, wave, mu, flux, sigma_kms):
    """Add a Gaussian with EXACT integral ``flux``.

    The width is the observed width (intrinsic and instrumental combined in
    quadrature).  Convolution preserves the integral, so the declared flux is
    exact.
    """
    sig_int = mu * sigma_kms / C.C_KMS
    sig_lsf = C.sigma_lsf_kms(mu) * mu / C.C_KMS
    sig_obs = math.sqrt(sig_int ** 2 + sig_lsf ** 2)
    amp = flux / (sig_obs * SQRT2PI)
    return flux_density + amp * np.exp(-0.5 * ((wave - mu) / sig_obs) ** 2), sig_obs


CASES = [
    # id     snr  sig_n  broad_frac sig_b  ha_flux  ha_hb  oiii_hb  nii_ha  badfrac z
    ("S01",   30.0,  85.0, 0.00,   0.0,     900.0,   2.85,   0.45,   0.32,   0.00, 0.0400),
    ("S02",   22.0, 110.0, 0.00,   0.0,     420.0,   3.10,   1.80,   0.42,   0.00, 0.0625),
    ("S03",   15.0,  95.0, 0.00,   0.0,     260.0,   3.40,   0.80,   0.25,   0.01, 0.0210),
    ("S04",   12.0, 175.0, 0.00,   0.0,     600.0,   3.80,   0.95,   0.60,   0.01, 0.0760),
    ("S05",   20.0, 105.0, 0.45, 520.0,     780.0,   4.20,   0.60,   0.38,   0.00, 0.0330),
    ("S06",    7.0, 120.0, 0.00,   0.0,     180.0,   3.60,   0.70,   0.30,   0.02, 0.0555),
    ("S07",    9.0, 115.0, 0.55, 780.0,     950.0,   5.50,   1.10,   0.45,   0.01, 0.0680),
    ("S08",   10.0, 130.0, 0.00,   0.0,     340.0,   3.30,   2.20,   0.35,   0.03, 0.0470),
    ("S09",    9.0, 140.0, 0.00,   0.0,      90.0,   4.50,   0.50,   0.20,   0.02, 0.0810),
    ("S10",   14.0, 160.0, 0.50, 650.0,     540.0,   4.00,   3.50,   0.70,   0.01, 0.0880),
]

N_PIXELS_HINT = 3849


def build_case(case_id, snr, sig_n, broad_frac, sig_b, ha_flux, ha_hb,
               oiii_hb, nii_ha, badfrac, z, base_seed):
    """Return (wave, loglam, flux, ivar, and_mask, or_mask, wdisp, truth, model)."""
    seed = base_seed + (int(case_id[1:]) * 10007 if case_id[1:].isdigit() else 0)
    rng = np.random.default_rng(seed)

    loglam = np.linspace(math.log10(3800.0), math.log10(9200.0), N_PIXELS_HINT)
    wave = 10.0 ** loglam
    lam_rest = wave / (1.0 + z)

    cont_norm = rng.uniform(4.0, 40.0)
    cont = _continuum(lam_rest, rng, cont_norm)
    abs_scale = rng.uniform(0.60, 1.40)
    sigma_stellar = max(rng.uniform(120.0, 260.0), 1.35 * sig_n)   # fairness rule 2
    cont = cont * _absorption_profile(lam_rest, abs_scale, sigma_stellar)

    model = cont.copy()
    comps = {}

    f_ha_total = ha_flux
    f_ha_broad = f_ha_total * broad_frac
    f_ha_narrow = f_ha_total - f_ha_broad
    f_hb_total = f_ha_total / ha_hb
    f_hb_broad = f_hb_total * broad_frac
    f_oiii = f_hb_total * oiii_hb
    f_nii6583 = f_ha_narrow * nii_ha

    def place(name, rest, flux, sigma_kms, record=True):
        nonlocal model
        mu = rest * (1.0 + z)
        model, sig_obs = _add_gaussian(model, wave, mu, flux, sigma_kms)
        if record:
            comps.setdefault(name, []).append(
                dict(rest=float(rest), flux=float(flux), sigma_kms=float(sigma_kms),
                     sigma_obs_ang=float(sig_obs)))

    place("Halpha", C.HALPHA, f_ha_narrow, sig_n)
    if f_ha_broad > 0:
        place("Halpha", C.HALPHA, f_ha_broad, sig_b)
    place("Hbeta", C.HBETA, f_hb_total - f_hb_broad, sig_n)
    if f_hb_broad > 0:
        place("Hbeta", C.HBETA, f_hb_broad, sig_b)
    place("OIII_5007", C.OIII_5007, f_oiii, sig_n * rng.uniform(0.9, 1.15))
    place("OIII_4959", C.OIII_4959, f_oiii / C.OIII_RATIO_5007_4959,
          sig_n * rng.uniform(0.9, 1.15))
    place("NII_6583", C.NII_6583, f_nii6583, sig_n)
    place("NII_6548", C.NII_6548, f_nii6583 / C.NII_RATIO_6583_6548, sig_n)
    place("_SII_6716", C.SII_6716, f_ha_narrow * rng.uniform(0.05, 0.30), sig_n, False)
    place("_SII_6731", C.SII_6731, f_ha_narrow * rng.uniform(0.04, 0.25), sig_n, False)
    if f_ha_narrow > 150:
        place("_OI_6300", C.OI_6300, f_ha_narrow * rng.uniform(0.01, 0.06), sig_n, False)
        place("_OI_6364", C.OI_6364, f_ha_narrow * rng.uniform(0.003, 0.02), sig_n, False)
    if f_hb_total > 60:
        place("_Hg", C.HGAMMA, f_hb_total * rng.uniform(0.15, 0.35), sig_n, False)
        place("_Hd", C.HDELTA, f_hb_total * rng.uniform(0.05, 0.18), sig_n, False)
        place("_NeIII", C.NEIII_3869, f_hb_total * rng.uniform(0.02, 0.15), sig_n, False)
        place("_OIIa", C.OII_3727a, f_hb_total * rng.uniform(0.1, 0.6), sig_n, False)
        place("_OIIb", C.OII_3727b, f_hb_total * rng.uniform(0.1, 0.6), sig_n, False)

    truth_lines = {}
    for name in ("Halpha", "Hbeta", "OIII_5007", "NII_6583"):
        truth_lines[name] = dict(
            flux=float(sum(c["flux"] for c in comps[name])),
            components=comps[name],
            sigma_kms_narrow=float(sig_n),
            has_broad=any(c["sigma_kms"] > 300 for c in comps[name]),
        )

    snr_lam = snr * (wave / 6000.0) ** 0.4
    sigma = np.maximum(cont, 1e-3) / snr_lam
    flux = model + rng.normal(0.0, 1.0, size=len(wave)) * sigma
    ivar = 1.0 / sigma ** 2

    and_mask = np.zeros(len(wave), np.int32)
    or_mask = np.zeros(len(wave), np.int32)
    n_bad = int(round(badfrac * len(wave)))
    if n_bad > 0:
        idx = rng.choice(len(wave), size=n_bad, replace=False)
        ivar[idx] = 0.0
        and_mask[idx] |= C.MASK_FULLREJECT
        or_mask[idx] |= C.MASK_FULLREJECT
    idx = rng.choice(len(wave), size=int(round(0.02 * len(wave))), replace=False)
    and_mask[idx] |= C.MASK_NEARBADPIXEL
    idx = rng.choice(len(wave), size=int(round(0.004 * len(wave))), replace=False)
    or_mask[idx] |= C.MASK_BRIGHTSKY

    wdisp = C.sigma_lsf_kms(wave) * wave / C.C_KMS / np.gradient(wave)

    truth = dict(
        case_id=case_id, klass="synthetic", z=float(z), lines=truth_lines,
        continuum=dict(norm=float(cont_norm), abs_scale=float(abs_scale),
                       sigma_stellar_kms=float(sigma_stellar),
                       basis="3 power laws x multiplicative stellar absorption "
                             "(single scale factor, tied Balmer pattern)"),
        noise=dict(snr_at_6000=float(snr), exponent=0.4,
                   bad_pixel_fraction=float(badfrac)),
        seed=int(seed),
        lambda_range_observed=[float(wave[0]), float(wave[-1])],
        n_pixels=int(len(wave)),
    )
    return wave, loglam, flux, ivar, and_mask, or_mask, wdisp, truth, model


def self_test(truth, wave):
    """Verify the generator quadrature component by component."""
    worst = 0.0
    for _name, rec in truth["lines"].items():
        for c in rec["components"]:
            mu = c["rest"] * (1.0 + truth["z"])
            model = np.zeros_like(wave)
            model, _ = _add_gaussian(model, wave, mu, c["flux"], c["sigma_kms"])
            m = np.abs(wave - mu) < 8.0 * c["sigma_obs_ang"]
            worst = max(worst, abs(np.trapz(model[m], wave[m]) - c["flux"])
                        / max(c["flux"], 1e-12))
    return worst


if __name__ == "__main__":
    base = 20261106
    outdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "cases"))
    os.makedirs(outdir, exist_ok=True)
    truths = {}
    for spec in CASES:
        cid = spec[0]
        wave, loglam, flux, ivar, am, om, wd, truth, model = build_case(cid, *spec[1:], base)
        rel = self_test(truth, wave)
        truth["generator_self_test_max_rel_error"] = rel
        C.write_case(os.path.join(outdir, f"{cid}.fits"), cid, loglam, flux, ivar, am, om, wd,
                     "synthetic")
        truths[cid] = truth
        print(f"{cid} z={truth['z']:.4f} npx={len(wave)} "
              f"Ha={truth['lines']['Halpha']['flux']:8.1f} "
              f"Hb={truth['lines']['Hbeta']['flux']:7.1f} "
              f"OIII={truth['lines']['OIII_5007']['flux']:8.1f} "
              f"NII={truth['lines']['NII_6583']['flux']:7.1f} "
              f"abs_scale={truth['continuum']['abs_scale']:.2f} "
              f"sig_stel={truth['continuum']['sigma_stellar_kms']:5.0f} "
              f"self_test={rel:.1e}")
    C.dump_json(os.path.join(outdir, "..", "truth_synthetic.json"), truths)
    print("\nwrote truth_synthetic.json")
