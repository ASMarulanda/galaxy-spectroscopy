#!/usr/bin/env python3
"""Figure 4 — injection response: does the pipeline transfer added flux?

    .venv/bin/python validation/figures/fig4_injection_recovery.py    # ~1 min

Input
    data/cases/I01..I04.fits and their base spectra R01, R02(base R03), ...
    data/truth/truth_injections.json   (the *known* added flux per line)
Output
    validation/figures/fig4_injection_recovery.{png,pdf}

Method. Each injection spectrum is its base spectrum plus one Gaussian line of
known flux, added at 40 % of the base pipeline flux. The measured increment is

    ΔF_measured = F(I_case) − F(base_case)

and the figure plots it against the known ΔF. This is a *differential* test: the
continuum, the absorber strength and the neighbouring lines are identical in both
spectra, so a common error in the continuum level cancels. It therefore measures
whether the fitter responds correctly to more line flux — not whether it gets the
absolute scale right (see the caveat printed on the figure).

Error bars are sqrt(σ_I² + σ_base²) from the two fits. As a estimate of the
uncertainty on the difference this is conservative: the two spectra share
everything except the added line, so their errors are strongly correlated and the
true uncertainty on the increment is smaller than the quadrature sum.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import gspec                                       # noqa: E402
from gspec import plotting                         # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LINES = ("Halpha", "Hbeta", "OIII_5007")


def main() -> None:
    inj = json.loads((ROOT / "data/truth/truth_injections.json").read_text())
    rows: list[dict] = []

    for case_id in sorted(inj):
        base = inj[case_id]["base_case"]
        r_inj = gspec.measure(gspec.read_case(ROOT / f"data/cases/{case_id}.fits"))
        r_base = gspec.measure(gspec.read_case(ROOT / f"data/cases/{base}.fits"))
        for line in LINES:
            known = inj[case_id]["injected"][line]["flux"]
            meas = r_inj.lines[line].flux - r_base.lines[line].flux
            err = math.hypot(r_inj.lines[line].flux_err, r_base.lines[line].flux_err)
            rows.append(dict(case_id=case_id, base=base, line=line,
                             known_increment=known, measured_increment=meas,
                             increment_err=err))
            print(f"  {case_id}/{base} {line:10s} known {known:8.2f}  "
                  f"measured {meas:8.2f}  ratio {meas/known:5.3f}")

    fig = plotting.fig_injection(rows)
    for p in plotting.save(fig, str(HERE / "fig4_injection_recovery")):
        print(f"wrote {p.relative_to(ROOT)}")

    print("\n  median recovery fraction by line:")
    for line in LINES:
        fr = sorted(r["measured_increment"] / r["known_increment"]
                    for r in rows if r["line"] == line)
        med = 0.5 * (fr[1] + fr[2]) if len(fr) == 4 else fr[len(fr) // 2]
        print(f"    {line:10s} {100*med:6.1f} %   (n={len(fr)})")


if __name__ == "__main__":
    main()
