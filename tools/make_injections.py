"""Build the real-data cases: redacted SDSS spectra + differential injections.

TRUTH SIDE.  Two products:

* ``R01..R06``  -- redacted copies of six Phase 1 SDSS DR18 spectra.
* ``I01..I04``  -- the same spectra with a synthetic line of EXACTLY known
  integrated flux added on top, at H-alpha, H-beta and [O III] 5007.

Redaction (why the case files are not the raw archive files)
------------------------------------------------------------
The archive files carry (a) a ``SPZLINE`` extension holding the SDSS pipeline's
own line fluxes and (b) ``SPECOBJ.Z`` holding the catalogue redshift.  Shipping
either would let an agent answer without measuring.  The case files therefore
keep only the calibrated spectrum (COADD, six columns) plus the unit/calibration
metadata, and carry an opaque ``CASEID``.  Plate/MJD/fibre/RA/Dec are removed so
the object cannot be looked up from memory either.  Provenance is retained in
this module and in the truth manifest, which lives in the verifier container.

The injected flux is exact by construction: a Gaussian of width
``sigma_obs = sqrt(sigma_int^2 + sigma_lsf^2)`` normalised so that its integral
equals the declared value.  Because the same noise realisation is used for the
original and the injected copy, the *difference* of two correct measurements is
the injected line alone -- an exact, catalogue-free ground truth on real data.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat as C  # noqa: E402  (import shim; everything below is verbatim)

SQRT2PI = math.sqrt(2.0 * math.pi)
PHASE1 = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "phase1"))
CASES = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "cases"))
TRUTH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

#: case_id -> (phase1 filename, role, fallback injection sigma_v in km/s or None).
#:
#: R05 has NO injection copy.  Phase 2 measured that a +40 % H-alpha increment
#: into this broad-line AGN is not reproducible: three independent method
#: configurations recovered +55 %, +99 % and +100 % of it, the full fit sits at
#: chi2/dof = 2.89, and the narrow-line amplitude moves by 44 % when the
#: increment is 16 % of the line.  A case whose truth cannot be recovered by any
#: valid method is not a verification case (contract §11 A7); broad-line
#: recovery is recorded as an open limitation instead.
SELECTION = [
    ("R01", "spec-0408-51821-0446.fits", "clean star-forming, high S/N, weak blend", 130.0),
    ("R02", "spec-0406-51817-0490.fits", "strong lines, low metallicity", 130.0),
    ("R03", "spec-0379-51789-0398.fits", "high redshift, clean narrow lines", 130.0),
    ("R04", "spec-0418-51817-0546.fits", "high stellar dispersion, blended", 130.0),
    ("R05", "spec-0968-52412-0011.fits", "broad-line AGN, severe Halpha/NII blend", None),
    ("R06", "spec-0294-51986-0036.fits", "low S/N, Hbeta formally undetected", None),
]


def read_phase1(fn):
    with fits.open(os.path.join(PHASE1, "data", fn)) as f:
        d = f[1].data
        so = f[2].data[0]
        out = dict(
            flux=np.asarray(d["flux"], float),
            loglam=np.asarray(d["loglam"], float),
            ivar=np.asarray(d["ivar"], float),
            and_mask=np.asarray(d["and_mask"], np.int32),
            or_mask=np.asarray(d["or_mask"], np.int32),
            wdisp=np.asarray(d["wdisp"], float),
            z=float(so["Z"]), z_err=float(so["Z_ERR"]),
            plate=int(so["PLATE"]), mjd=int(so["MJD"]), fiber=int(so["FIBERID"]),
            specobjid=str(so["SPECOBJID"]).strip(),
            subclass=str(so["SUBCLASS"]).strip(),
            sn=float(so["SN_MEDIAN_ALL"]), vdisp=float(so["VDISP"]),
            spzline={str(r["LINENAME"]).strip(): float(r["LINEAREA"]) for r in f[3].data},
            spzsigma={str(r["LINENAME"]).strip(): float(r["LINESIGMA"]) for r in f[3].data},
        )
    out["wave"] = 10.0 ** out["loglam"]
    return out


def add_line(flux, wave, rest, z, total_flux, sigma_int_kms):
    """Add a Gaussian whose integral is exactly ``total_flux``.  Returns (flux, meta)."""
    mu = rest * (1.0 + z)
    sig_int = mu * sigma_int_kms / C.C_KMS
    sig_lsf = C.sigma_lsf_kms(mu) * mu / C.C_KMS
    sig = math.sqrt(sig_int ** 2 + sig_lsf ** 2)
    amp = total_flux / (sig * SQRT2PI)
    return flux + amp * np.exp(-0.5 * ((wave - mu) / sig) ** 2), dict(
        rest=float(rest), mu_observed=float(mu), flux=float(total_flux),
        sigma_int_kms=float(sigma_int_kms), sigma_obs_ang=float(sig),
    )


def main():
    os.makedirs(CASES, exist_ok=True)
    cand = {f"{r['plate']}/{r['fiber']}": r for r in
            json.load(open(os.path.join(PHASE1, "logs", "candidates10.json")))}
    truth = {}
    inj_truth = {}

    for cid, fn, role, inj_sigma in SELECTION:
        d = read_phase1(fn)
        key = f"{d['plate']}/{d['fiber']}"
        mpa = cand[key]
        C.write_case(os.path.join(CASES, f"{cid}.fits"), cid, d["loglam"], d["flux"],
                     d["ivar"], d["and_mask"], d["or_mask"], d["wdisp"],
                     "sdss-dr18-redacted")
        truth[cid] = dict(
            case_id=cid, klass="real",
            z=d["z"], z_err=d["z_err"], role=role,
            archive=dict(release="DR18", run2d="26", plate=d["plate"], mjd=d["mjd"],
                         fiber=d["fiber"], specobjid=d["specobjid"],
                         subclass=d["subclass"], sn_median=d["sn"],
                         vdisp_kms=d["vdisp"], file_sha256_source=fn),
            lines={
                # strict grading does NOT use these; they are a broad consistency band
                "Halpha": dict(pipe=d["spzline"].get("H_alpha"),
                               mpa=mpa["mpa_Ha"], mpa_err=mpa["mpa_Ha_err"]),
                "Hbeta": dict(pipe=d["spzline"].get("H_beta"),
                              mpa=mpa["mpa_Hb"], mpa_err=mpa["mpa_Hb_err"]),
                "OIII_5007": dict(pipe=d["spzline"].get("[O_III] 5007"),
                                  mpa=mpa["mpa_OIII5007"], mpa_err=mpa["mpa_OIII5007_err"]),
                "NII_6583": dict(pipe=d["spzline"].get("[N_II] 6583"),
                                 mpa=mpa["mpa_NII6583"], mpa_err=mpa["mpa_NII6583_err"]),
            },
            n_pixels=int(len(d["wave"])),
            lambda_range_observed=[float(d["wave"][0]), float(d["wave"][-1])],
        )
        # ---- injection copy -------------------------------------------------
        if inj_sigma is None:
            continue
        f_ha = d["spzline"].get("H_alpha", 0.0) or 0.0
        f_hb = d["spzline"].get("H_beta", 0.0) or 0.0
        f_o3 = d["spzline"].get("[O_III] 5007", 0.0) or 0.0
        if f_hb <= 5.0:                       # G03-style non-detection
            f_hb = max(f_hb, 0.0) + f_ha / 3.2
        if f_o3 <= 5.0:
            f_o3 = f_ha / 8.0
        d_ha = 0.40 * max(f_ha, 50.0)
        d_hb = 0.40 * f_hb
        d_o3 = 0.40 * f_o3

        # Inject at the object's OWN measured line width (SDSS pipeline
        # LINESIGMA).  A fixed width would create an artificial double-component
        # profile that no single-Gaussian model can fit, which Phase 2 measured
        # to cost 15 % on [O III] for a reason that has nothing to do with the
        # science being tested (failure F-I1).
        def own_sigma(name, fallback):
            v = d["spzsigma"].get(name)
            return float(v) if v and 20.0 < v < 1200.0 else fallback

        s_ha = own_sigma("H_alpha", inj_sigma)
        s_hb = own_sigma("H_beta", 0.85 * s_ha)
        s_o3 = own_sigma("[O_III] 5007", 0.85 * s_ha)

        flux = d["flux"].copy()
        flux, meta_ha = add_line(flux, d["wave"], C.HALPHA, d["z"], d_ha, s_ha)
        flux, meta_hb = add_line(flux, d["wave"], C.HBETA, d["z"], d_hb, s_hb)
        flux, meta_o3 = add_line(flux, d["wave"], C.OIII_5007, d["z"], d_o3, s_o3)
        # keep the [O III] doublet ratio physical
        flux, meta_o3959 = add_line(flux, d["wave"], C.OIII_4959, d["z"],
                                    d_o3 / C.OIII_RATIO_5007_4959, s_o3)

        iid = "I" + cid[1:]
        C.write_case(os.path.join(CASES, f"{iid}.fits"), iid, d["loglam"], flux,
                     d["ivar"], d["and_mask"], d["or_mask"], d["wdisp"],
                     "sdss-dr18-redacted")
        inj_truth[iid] = dict(
            case_id=iid, klass="injection", base_case=cid, z=d["z"], z_err=d["z_err"],
            injected={"Halpha": meta_ha, "Hbeta": meta_hb, "OIII_5007": meta_o3,
                      "OIII_4959": meta_o3959},
            n_pixels=int(len(d["wave"])),
            lambda_range_observed=[float(d["wave"][0]), float(d["wave"][-1])],
        )
        truth[iid] = dict(case_id=iid, klass="injection", z=d["z"], z_err=d["z_err"],
                          base_case=cid, n_pixels=int(len(d["wave"])),
                          lambda_range_observed=[float(d["wave"][0]), float(d["wave"][-1])])

    C.dump_json(os.path.join(TRUTH, "truth_real.json"), truth)
    C.dump_json(os.path.join(TRUTH, "truth_injections.json"), inj_truth)
    print(f"wrote {len(truth)} real cases, {len(inj_truth)} injection cases")
    for k, v in inj_truth.items():
        print(f"  {k} (base {v['base_case']}): Halpha +{v['injected']['Halpha']['flux']:.1f} "
              f"Hbeta +{v['injected']['Hbeta']['flux']:.1f} "
              f"OIII5007 +{v['injected']['OIII_5007']['flux']:.1f} "
              f"sigma_int={v['injected']['Halpha']['sigma_int_kms']:.0f} km/s")


if __name__ == "__main__":
    main()
