#!/usr/bin/env python3
"""Figure 2 — Hβ depends on the fitting configuration.

    .venv/bin/python validation/figures/fig2_hbeta_configurations.py     # ~4 min

Input
    data/cases/S01..S10.fits  (synthetic, analytic truth known)
    data/truth/truth_synthetic.json
    validation/results/out_independent.csv  (the separate implementation)
Output
    validation/figures/fig2_hbeta_configurations.{png,pdf}

The four configurations, and what each one isolates
    1. ``BIC selection (default)``  the shipped pipeline: cubic continuum,
       absorption available, 7-candidate BIC sweep.
    2. ``no selection``             identical constants, sweep switched off. The
       difference between 1 and 2 is *the effect of model selection alone*.
    3. ``m1b``                      quadratic continuum, windows 8 Å wider, Hγ
       weight 1.05, no sweep. Same fitter, different constants.
    4. ``independent``              a separately written program (RBF continuum,
       bounded optimiser, bootstrap) — not a variant of this code.

Configurations 2 and 3 are therefore not independent draws from a model space;
they are three points in it plus one external method. The figure is a sensitivity
map, not a model comparison in the statistical sense.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import gspec                                          # noqa: E402
from gspec import plotting, variants                  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SYNT = [f"S{i:02d}" for i in range(1, 11)]


def main() -> None:
    truth = json.loads((ROOT / "data/truth/truth_synthetic.json").read_text())
    t = {c: truth[c]["lines"]["Hbeta"]["flux"] for c in SYNT}

    configs: dict[str, dict[str, float]] = {}

    # 1 — the shipped pipeline
    sel = {}
    for c in SYNT:
        r = gspec.measure(gspec.read_case(ROOT / f"data/cases/{c}.fits"))
        sel[c] = r.lines["Hbeta"].flux
        print(f"  {c}  selection    Hbeta = {sel[c]:8.2f}  ({r.model.label()})")
    configs["BIC selection (default)"] = sel

    # 2 — same constants, sweep off: isolates the effect of model selection
    nosel = {}
    for c in SYNT:
        r = gspec.measure(gspec.read_case(ROOT / f"data/cases/{c}.fits"),
                          model=variants.M1A, model_select=False)
        nosel[c] = r.lines["Hbeta"].flux
        print(f"  {c}  no selection Hbeta = {nosel[c]:8.2f}")
    configs["no selection"] = nosel

    # 3 — the alternative configuration
    m1b = {}
    for c in SYNT:
        r = gspec.measure(gspec.read_case(ROOT / f"data/cases/{c}.fits"),
                          model=variants.M1B)
        m1b[c] = r.lines["Hbeta"].flux
        print(f"  {c}  m1b          Hbeta = {m1b[c]:8.2f}")
    configs["m1b"] = m1b

    # 4 — the separate implementation, read from its committed output
    import csv as _csv
    ind = {}
    with (ROOT / "validation/results/out_independent.csv").open(newline="") as fh:
        for row in _csv.DictReader(fh):
            if row["line"] == "Hbeta":
                ind[row["case_id"]] = float(row["flux"])
    configs["independent"] = {c: ind[c] for c in SYNT}

    fig = plotting.fig_hbeta_configurations(t, configs)
    for p in plotting.save(fig, str(HERE / "fig2_hbeta_configurations")):
        print(f"wrote {p.relative_to(ROOT)}")

    print("\n  median signed error per configuration:")
    for lab, d in configs.items():
        errs = sorted(100 * (d[c] - t[c]) / t[c] for c in SYNT)
        med = 0.5 * (errs[4] + errs[5])
        print(f"    {lab:24s} {med:+7.1f} %   worst {max(errs, key=abs):+8.1f} %")


if __name__ == "__main__":
    main()
