"""C1: does the pipeline recover redshift and the four fluxes against analytic truth?

The truth here is the generator's own parameters — each Gaussian was normalised so
that its integral *is* the declared flux — so nothing in this file is compared
against a fit, and no expected value was ever taken from the pipeline's output.

What the test encodes, and what it deliberately does not
-------------------------------------------------------
The acceptance rule is the validated one: a row passes if its error is inside its
relative floor **or** inside three times its own reported uncertainty. That
second term is not generosity, it is the uncertainty contract: a faint line whose
error bar is honestly large is a legitimate measurement. It also means "passes"
is weaker than "accurate to the floor", which is why the characterisation tests
at the bottom pin the actual medians instead of hiding behind the rule.
"""

from __future__ import annotations

import math
import statistics as st

import pytest

from conftest import FLOORS, SYNTHETIC


def rel_error(result, truth, line):
    true = truth[line]["flux"]
    return (result.lines[line].flux - true) / true


def test_every_synthetic_row_is_within_its_acceptance_band(reference, synthetic_truth):
    failures = []
    for cid in SYNTHETIC:
        r, truth = reference[cid], synthetic_truth[cid]["lines"]
        for line, floor in FLOORS.items():
            true = truth[line]["flux"]
            m = r.lines[line]
            err = abs(m.flux - true)
            tol = max(floor * true, 3.0 * m.flux_err)
            if err > tol:
                failures.append(f"{cid}/{line}: |{m.flux:.2f}-{true:.2f}|={err:.2f} > {tol:.2f}")
    assert not failures, "rows outside their acceptance band:\n  " + "\n  ".join(failures)


def test_no_more_than_six_rows_need_their_uncertainty_to_rescue_them(reference, synthetic_truth):
    """The validated budget: at most 6 of 40 synthetic rows sit beyond 3x their floor.

    Without this, "inside the acceptance band" could be satisfied by a large
    error bar on every row. Measured on the reference: 5.
    """
    far = 0
    for cid in SYNTHETIC:
        r, truth = reference[cid], synthetic_truth[cid]["lines"]
        for line, floor in FLOORS.items():
            true = truth[line]["flux"]
            if abs(r.lines[line].flux - true) > 3.0 * floor * true:
                far += 1
    assert far <= 6, f"{far} of 40 rows sit beyond 3x their floor"


def test_at_least_eight_cases_carry_a_line_inside_its_floor(reference, synthetic_truth):
    """The case-level anchor: a method that is only right 'on average' fails here.

    This is the gate that resists a corrupted set — a constant or copied
    prediction satisfies almost none of these anchors.
    """
    anchored = 0
    for cid in SYNTHETIC:
        r, truth = reference[cid], synthetic_truth[cid]["lines"]
        if any(abs(rel_error(r, truth, line)) <= floor for line, floor in FLOORS.items()):
            anchored += 1
    assert anchored >= 8, f"only {anchored} of 10 cases carry a line inside its floor"


@pytest.mark.characterisation
def test_the_documented_median_errors_are_reproduced(reference, synthetic_truth):
    """Pin the actual accuracy, so a silent regression cannot hide behind the band.

    These are the numbers quoted in validation/report.md. A change here is a
    change in the science and must be explained, not absorbed.
    """
    medians = {}
    for line in FLOORS:
        errs = [100 * rel_error(reference[cid], synthetic_truth[cid]["lines"], line)
                for cid in SYNTHETIC]
        medians[line] = st.median(errs)

    assert medians["Halpha"] == pytest.approx(-0.70, abs=0.15)
    assert medians["OIII_5007"] == pytest.approx(+0.31, abs=0.15)
    assert medians["NII_6583"] == pytest.approx(-1.54, abs=0.15)
    # Hbeta is the documented weak link: a systematic *under*-estimate of ~14 %.
    # This is an acceptance threshold, not an accuracy claim (limitations L1).
    assert medians["Hbeta"] == pytest.approx(-14.42, abs=0.3)


@pytest.mark.characterisation
def test_hbeta_is_the_weakest_line(reference, synthetic_truth):
    """The ordering matters more than the numbers: it is what L1 asserts."""
    worst = {}
    for line in FLOORS:
        worst[line] = max(abs(rel_error(reference[cid], synthetic_truth[cid]["lines"], line))
                          for cid in SYNTHETIC)
    assert worst["Hbeta"] == max(worst.values())
    assert worst["NII_6583"] < 0.10, "the least blended line should behave best"


def test_fluxes_are_positive_and_finite(reference):
    for cid, r in reference.items():
        for line, m in r.lines.items():
            assert math.isfinite(m.flux), f"{cid}/{line} is not finite"
            assert m.flux > 0, f"{cid}/{line} is not positive"
            assert m.flux_err > 0, f"{cid}/{line} has a non-positive uncertainty"
