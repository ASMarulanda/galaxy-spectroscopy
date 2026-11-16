#!/usr/bin/env python3
"""Figure 1 — the S01 fit, decomposed.

    .venv/bin/python validation/figures/fig1_s01_fit.py [CASE_ID]

Input
    data/cases/S01.fits — one SDSS-format spectrum (synthetic S01 by default;
    any other case id works). Loaded by ``gspec.io.read_case``.
Output
    validation/figures/fig1_s01_fit.png
    validation/figures/fig1_s01_fit.pdf

What the figure shows, panel by panel
    (a) the whole spectrum with the two fitted windows marked;
    (b) the Hα complex zoomed, with each fitted Gaussian drawn *on top of* the
        continuum so its contribution is visible, and the total model over them;
    (c) the same for Hβ + [O III];
    (d) normalised residuals (data − model)/σ over every fitted pixel, which is
        the panel that decides whether the fit is acceptable at all.

Reading it: if the continuum (dashed) sits above the data in the line cores, the
Gaussians are doing the work; if the residuals in (d) show a coherent dip under
Hβ, the model is missing absorption rather than noise.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import gspec                                     # noqa: E402
from gspec import plotting                       # noqa: E402

HERE = Path(__file__).resolve().parent
CASE = sys.argv[1] if len(sys.argv) > 1 else "S01"


def main() -> None:
    path = HERE.parents[1] / "data" / "cases" / f"{CASE}.fits"
    case = gspec.read_case(path)
    result = gspec.measure(case)

    print(f"{CASE}: z = {result.redshift:.6f}  ({result.model.label()})")
    print(f"      chi2/dof = {result.chi2dof:.4f}   parameters = {result.n_par}")
    for ln in ("Halpha", "Hbeta", "OIII_5007", "NII_6583"):
        m = result.lines[ln]
        print(f"      {ln:10s} F = {m.flux:9.2f} +- {m.flux_err:5.2f}   {m.flag}")

    fig = plotting.fig_spectrum_fit(case, result)
    stem = HERE / f"fig1_{CASE.lower()}_fit"
    for p in plotting.save(fig, str(stem)):
        print(f"wrote {p.relative_to(HERE.parents[1])}")


if __name__ == "__main__":
    main()
