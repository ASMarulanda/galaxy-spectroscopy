"""The figures are only trustworthy if they are the model that was fitted.

A plot can be wrong in ways no unit test of the fitter would catch: it can draw a
decomposition that does not sum to the model, use a continuum evaluated at the
wrong redshift, or quietly label a non-detection as a measurement. These tests
attack exactly that.

The central one is `test_curves_reproduce_the_fitted_chi2`: the arrays returned by
`model_curves` are re-evaluated from `result.parameters` alone, and their residual
sum of squares must equal the `chi2` the optimiser reported. If the plotting layer
ever drifts from the fitter — a changed parameter order, a wrong window, a
forgotten absorption term — that identity breaks.
"""

from __future__ import annotations

import math

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gspec import diagnostics, plotting                      # noqa: E402
from conftest import REAL                                    # noqa: E402


@pytest.fixture(scope="module")
def s01(cases, reference):
    return cases["S01"], reference["S01"]


def test_decomposition_sums_to_the_total(s01):
    """continuum + Σ components == total, elementwise, where both are defined."""
    case, result = s01
    cur = plotting.model_curves(case, result)
    tot = cur["continuum"].copy()
    for _key, _rest, prof in cur["components"]:
        tot = tot + np.where(np.isfinite(prof), prof, 0.0)
    ok = np.isfinite(cur["total"])
    assert ok.sum() > 1000
    assert np.allclose(tot[ok], cur["total"][ok], rtol=0, atol=1e-12)


def test_curves_reproduce_the_fitted_chi2(s01):
    """The replotted model is the fitted model: same residuals, same chi²."""
    case, result = s01
    cur = plotting.model_curves(case, result)
    m = cur["in_window"] & np.isfinite(cur["total"])
    r = (case.flux[m] - cur["total"][m]) * np.sqrt(case.ivar[m])
    assert math.isclose(float(np.sum(r ** 2)), result.chi2, rel_tol=1e-8), (
        f"replotted chi2 {np.sum(r**2):.6f} != fitted {result.chi2:.6f}")


def test_components_are_gaussians_of_the_right_scale(s01):
    """Each component integrates to A·σ·√(2π), the flux the pipeline reports."""
    case, result = s01
    cur = plotting.model_curves(case, result)
    com = {k: p for k, _r, p in cur["components"]}
    ha = com["Ha_n"]
    m = cur["in_window"]
    integral = float(np.nansum(ha[m] * case.dlam[m]))
    line = np.nansum([c.flux for c in result.lines["Halpha"].components
                      if c.key in ("Ha_n", "Ha_b")])
    assert math.isclose(integral, line, rel_tol=2e-3), (
        f"integrated component {integral:.2f} vs reported {line:.2f}")


def test_figure_has_four_populated_panels(s01, tmp_path):
    case, result = s01
    fig = plotting.fig_spectrum_fit(case, result)
    assert len(fig.axes) == 4
    for ax in fig.axes:
        assert ax.lines or ax.collections, "an empty panel is not a figure"
    paths = plotting.save(fig, str(tmp_path / "f"))
    assert all(p.stat().st_size > 5000 for p in paths)


def test_figures_are_byte_reproducible(s01, tmp_path):
    """Documented property: no timestamps, so a regenerated figure can be diffed."""
    case, result = s01
    a = plotting.save(plotting.fig_spectrum_fit(case, result), str(tmp_path / "a"))[0]
    b = plotting.save(plotting.fig_spectrum_fit(case, result), str(tmp_path / "b"))[0]
    assert a.read_bytes() == b.read_bytes()


def test_nondetections_are_not_drawn_as_measurements(cases, reference):
    """R03 has Hα but no Hβ: it must not appear as a Balmer-decrement point."""
    rows = diagnostics.diagnostic_table({c: reference[c] for c in REAL})
    by_case = {r["case_id"]: r for r in rows}
    assert not by_case["R03"]["balmer"]["detected"]

    fig = plotting.fig_diagnostics(rows)
    ax_b = fig.axes[1]
    # the measured points are the only artists drawn with a circular marker;
    # cap lines, the case-B guide and the non-detection crosses are not
    markers = [ln for ln in ax_b.lines if ln.get_marker() == "o"]
    ys = [v for ln in markers for v in ln.get_ydata()]
    n_detected = sum(1 for r in rows if r["balmer"]["detected"])
    assert n_detected == 2
    assert len(ys) == n_detected, f"{len(ys)} points drawn for {n_detected} detections"
    assert all(np.isfinite(ys))
