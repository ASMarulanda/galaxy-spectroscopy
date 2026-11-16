#!/usr/bin/env python3
"""Figure 3 — are the reported uncertainties honest?

    .venv/bin/python validation/figures/fig3_pull_calibration.py     # ~2 s

Input
    validation/results/out_gspec.csv         (this pipeline, 10 synthetic cases)
    validation/results/out_independent.csv   (the separate implementation)
    data/truth/truth_synthetic.json          (analytic truth)
Output
    validation/figures/fig3_pull_calibration.{png,pdf}

Definition. For a measurement F with reported uncertainty σ_F and known true
value F_true, the pull is

    p = (F − F_true) / σ_F

If σ_F is the standard deviation of the measurement error, the pulls are drawn
from N(0,1): their RMS is ≈ 1, and 68 % lie within ±1. An RMS well above 1 means
the errors are *under*estimated; well below 1 means they are too conservative.
This is a necessary condition for calibration, not a sufficient one — see the
caveat printed on the figure.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from gspec import plotting                              # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SYNT = [f"S{i:02d}" for i in range(1, 11)]
LINES = ("Halpha", "Hbeta", "OIII_5007", "NII_6583")


def load_pulls(path: Path, truth: dict) -> list[dict]:
    import csv

    rows, pulls = [], []
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append((row["case_id"], row["line"],
                         float(row["flux"]), float(row["flux_err"])))
    for case, line, flux, err in rows:
        if case not in truth or line not in truth[case]["lines"]:
            continue
        if err <= 0 or not math.isfinite(err):
            continue
        t = truth[case]["lines"][line]["flux"]
        pulls.append(dict(case=case, line=line, pull=(flux - t) / err))
    return pulls


def main() -> None:
    truth = json.loads((ROOT / "data/truth/truth_synthetic.json").read_text())
    sets = {}
    for label, fname in (("gspec", "golden_reference.csv"),
                         ("independent", "out_independent.csv")):
        path = ROOT / "validation/results" / fname
        if not path.exists():
            print(f"missing {path}; run the pipeline first")
            continue
        sets[label] = load_pulls(path, truth)
        vals = [r["pull"] for r in sets[label]]
        rms = math.sqrt(sum(v * v for v in vals) / len(vals))
        inside = sum(abs(v) <= 1 for v in vals) / len(vals)
        print(f"  {label:12s} n={len(vals):3d}  RMS={rms:.3f}  "
              f"median|p|={sorted(abs(v) for v in vals)[len(vals)//2]:.3f}  "
              f"within 1σ={100*inside:.0f}%  max|p|={max(abs(v) for v in vals):.2f}")

    fig = plotting.fig_pulls(sets)
    for p in plotting.save(fig, str(HERE / "fig3_pull_calibration")):
        print(f"wrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
