"""E-Abs: does the stellar-absorption term actually do anything?

The question
------------
The pipeline *contains* a multiplicative Balmer absorption model. Whether that
model improves the measurement is a separate, testable claim. Stage 1 found that
the BIC sweep rejects it on all ten synthetic spectra (docs/methods.md §4.2).
This experiment asks the follow-up questions the sweep cannot answer:

    1. If the term is forced ON, how much do the fluxes move?
    2. Does the fit recover the absorption the generator actually put there?
    3. Is the basis the problem? The measurement side uses three absorbers with
       weights (Hgamma 0.95, G band 0.60, Hbeta 1.00); the generator uses the
       same three lines inside this window but with (1.00, 0.45, 1.00).

Design
------
Four configurations, all fitted at the **same redshift seed** (the value the
default pipeline reports) so that the only thing that changes between them is
the model:

    bic      the shipped pipeline: absorption on the table, BIC decides
    abs_on   3-line absorption forced on, ndeg=3
    abs_off  no absorption term, ndeg=3
    abs_gen  absorption forced on with the generator's own relative weights

Nothing here is imported by the package. The default pipeline is untouched.

Run:  python validation/experiments/e_abs_absorption.py
"""

from __future__ import annotations

import json
import statistics as st
import sys
from pathlib import Path

import gspec
from gspec import constants as K
from gspec.pipeline import measure, measure_fixed_model
from gspec.results import ModelConfig
from gspec.select import bic

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = Path(__file__).resolve().parent / "results"

SYNTHETIC = tuple(f"S{i:02d}" for i in range(1, 11))
FLOORS = {"Halpha": 0.07, "Hbeta": 0.20, "OIII_5007": 0.06, "NII_6583": 0.10}

#: The generator's absorption lines that fall inside the Hbeta window
#: (4250-5090 A rest), with its relative weights. Verified against
#: tools/generate_synthetic.py ABSORPTION_LINES.
GENERATOR_MATCHED = (
    (4341.69, 1.00, 1.0),      # Hgamma
    (4305.00, 0.45, 3.0),      # G band
    (4862.68, 1.00, 1.0),      # Hbeta
)

CONFIGS = {
    "abs_on": ModelConfig(variant="exp", ndeg=3, use_absorption=True),
    "abs_off": ModelConfig(variant="exp", ndeg=3, use_absorption=False),
    "abs_gen": ModelConfig(variant="exp", ndeg=3, use_absorption=True,
                           absorbers_hb=GENERATOR_MATCHED),
}


def main() -> int:
    truth = json.loads((DATA / "truth" / "truth_synthetic.json").read_text())
    cases = {cid: gspec.read_case(DATA / "cases" / f"{cid}.fits") for cid in SYNTHETIC}

    rows = []
    for cid, case in cases.items():
        # 1. the shipped pipeline, so the comparison has a reference point
        shipped = measure(case)
        z0 = shipped.redshift

        results = {"bic": shipped}
        for name, cfg in CONFIGS.items():
            results[name] = measure_fixed_model(case, cfg, z_known=z0)

        true_abs = truth[cid]["continuum"]["abs_scale"]
        true_sig_n = truth[cid]["lines"]["Halpha"]["sigma_kms_narrow"]
        true_sig_star = truth[cid]["continuum"]["sigma_stellar_kms"]

        for name, r in results.items():
            rec = dict(
                case=cid, config=name, n_par=r.n_par, chi2=round(r.chi2, 2),
                chi2dof=round(r.chi2dof, 4), bic=round(bic(r), 2),
                use_absorption=r.model.use_absorption,
                z=round(r.redshift, 8),
                abs_scale_fit=(round(r.parameters["abs_scale"], 4)
                               if "abs_scale" in r.parameters else None),
                abs_ratio_fit=(round(__import__("math").exp(r.parameters["abs_ratio"]), 3)
                               if "abs_ratio" in r.parameters else None),
                abs_scale_true=round(true_abs, 4),
                sigma_star_over_gas_true=round(true_sig_star / true_sig_n, 3),
            )
            for line in FLOORS:
                m = r.lines[line]
                err = (m.flux - truth[cid]["lines"][line]["flux"]) / truth[cid]["lines"][line]["flux"]
                rec[f"{line}_err_pct"] = round(100 * err, 3)
                rec[f"{line}_flux"] = round(m.flux, 4)
            rows.append(rec)
        print(f"  {cid}: " + "  ".join(
            f"{n}={results[n].chi2:.1f}/{results[n].n_par}par"
            for n in ("bic", "abs_on", "abs_off", "abs_gen")), flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "e_abs_absorption.json").write_text(json.dumps(rows, indent=1))

    # ------------------------------------------------------------- summary
    print("\n=== median error per line and configuration (%) ===")
    header = f"{'config':9}" + "".join(f"{ln:>12}" for ln in FLOORS) + f"{'mean|err|':>11}"
    print(header)
    for name in ("bic", "abs_on", "abs_off", "abs_gen"):
        sel = [r for r in rows if r["config"] == name]
        meds = [st.median(r[f"{ln}_err_pct"] for r in sel) for ln in FLOORS]
        meanabs = st.mean(abs(r[f"{ln}_err_pct"]) for r in sel for ln in FLOORS)
        print(f"{name:9}" + "".join(f"{m:12.2f}" for m in meds) + f"{meanabs:11.2f}")

    print("\n=== worst |error| per line (%) ===")
    print(header)
    for name in ("bic", "abs_on", "abs_off", "abs_gen"):
        sel = [r for r in rows if r["config"] == name]
        worst = [max(abs(r[f"{ln}_err_pct"]) for r in sel) for ln in FLOORS]
        print(f"{name:9}" + "".join(f"{w:12.2f}" for w in worst))

    print("\n=== does forcing absorption on pay for itself? (delta-BIC: abs_on - abs_off) ===")
    print(f"{'case':5} {'dchi2':>8} {'dpar':>5} {'dBIC':>8}  verdict")
    wins = 0
    for cid in SYNTHETIC:
        on = next(r for r in rows if r["case"] == cid and r["config"] == "abs_on")
        off = next(r for r in rows if r["case"] == cid and r["config"] == "abs_off")
        dchi = on["chi2"] - off["chi2"]
        dpar = on["n_par"] - off["n_par"]
        dbic = on["bic"] - off["bic"]
        verdict = "absorption wins" if dbic < 0 else "absorption LOSES"
        wins += dbic < 0
        print(f"{cid:5} {dchi:8.1f} {dpar:5d} {dbic:8.1f}  {verdict}")
    print(f"\n  absorption is preferred on {wins} of {len(SYNTHETIC)} synthetic cases")

    print("\n=== does the fit recover the absorption the generator put there? ===")
    print(f"{'case':5} {'abs_scale true':>15} {'fitted':>8} {'ratio':>7} "
          f"{'sigma* /sig_gas true':>21} {'fitted':>7}")
    for cid in SYNTHETIC:
        on = next(r for r in rows if r["case"] == cid and r["config"] == "abs_on")
        ratio = on["abs_scale_fit"] / on["abs_scale_true"] if on["abs_scale_true"] else float("nan")
        print(f"{cid:5} {on['abs_scale_true']:15.3f} {on['abs_scale_fit']:8.3f} "
              f"{ratio:7.2f} {on['sigma_star_over_gas_true']:21.2f} {on['abs_ratio_fit']:7.2f}")

    print(f"\nwrote {OUT / 'e_abs_absorption.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
