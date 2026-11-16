#!/usr/bin/env python3
"""Figure 5 — line-ratio diagnostics on the six real spectra.

    .venv/bin/python validation/figures/fig5_diagnostics.py      # ~40 s

Input
    data/cases/R01..R06.fits   (real SDSS spectra)
Output
    validation/figures/fig5_diagnostics.{png,pdf}
    validation/results/tables/diagnostic_ratios.csv

What is plotted
    (a) the [N II] BPT diagram. Only objects with *both* [N II]/Hα and
        [O III]/Hβ measured above 3σ can be placed; the rest are named in the
        panel with the reason. The two reference curves are drawn as references,
        not as a classification of the data.
    (b) Hα/Hβ per object, with the case-B value marked and non-detections shown
        as crosses rather than as points at some plausible value.

The secondary axis of (b) converts the ratio to A_V under a foreground-screen,
case-B assumption. It is a scale, not a measurement, and the assumptions it needs
are listed in ``docs/figures.md`` §5.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import gspec                                       # noqa: E402
from gspec import diagnostics, plotting            # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REAL = [f"R{i:02d}" for i in range(1, 7)]


def main() -> None:
    results = {}
    for case_id in REAL:
        results[case_id] = gspec.measure(gspec.read_case(ROOT / f"data/cases/{case_id}.fits"))

    rows = diagnostics.diagnostic_table(results)

    print(f"{'case':5} {'BPT':>4}  {'[NII]/Ha':>16} {'[OIII]/Hb':>16} {'Ha/Hb':>16}  A_V")
    for r in rows:
        def cell(d):
            return f"{d['value']:7.3f}±{d['error']:5.3f}" if d["detected"] else "   not det."
        ex = r["extinction"]
        av = f"{ex['av']:.2f}" if ex["valid"] else (f"<{ex['av_upper']:.2f}" if ex["av_upper"] == ex["av_upper"] else "--")
        print(f"{r['case_id']:5} {'yes' if r['bpt_placeable'] else 'no':>4}  "
              f"{cell(r['nii_ha']):>16} {cell(r['oiii_hb']):>16} {cell(r['balmer']):>16}  {av}")

    fig = plotting.fig_diagnostics(rows)
    for p in plotting.save(fig, str(HERE / "fig5_diagnostics")):
        print(f"wrote {p.relative_to(ROOT)}")

    out = ROOT / "validation/results/tables/diagnostic_ratios.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case_id", "halpha_snr", "hbeta_snr", "oiii_snr", "nii_snr",
                    "nii_ha", "nii_ha_err", "oiii_hb", "oiii_hb_err",
                    "balmer", "balmer_err", "bpt_placeable",
                    "e_bv", "e_bv_err", "av", "av_err", "av_upper", "extinction_note"])
        for r in rows:
            ex = r["extinction"]
            w.writerow([
                r["case_id"], f"{r['halpha_snr']:.2f}", f"{r['hbeta_snr']:.2f}",
                f"{r['oiii_snr']:.2f}", f"{r['nii_snr']:.2f}",
                f"{r['nii_ha']['value']:.6g}", f"{r['nii_ha']['error']:.6g}",
                f"{r['oiii_hb']['value']:.6g}", f"{r['oiii_hb']['error']:.6g}",
                f"{r['balmer']['value']:.6g}", f"{r['balmer']['error']:.6g}",
                int(r["bpt_placeable"]),
                f"{ex['e_bv']:.6g}", f"{ex['e_bv_err']:.6g}",
                f"{ex['av']:.6g}", f"{ex['av_err']:.6g}", f"{ex['av_upper']:.6g}",
                ex["reason"].replace(",", ";"),
            ])
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
