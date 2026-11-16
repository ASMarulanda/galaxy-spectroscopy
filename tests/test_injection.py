"""C2: differential recovery on real spectra.

This is the strongest test in the repository, because the truth it uses cannot be
derived from any measurement of the data:

    the copy and its base share **one noise realisation**, and the copy differs
    only by a Gaussian of exactly known integral, so

        F(copy) - F(base) = F_injected        exactly,

    whatever the stellar continuum does. The continuum, the absorption and the
    calibration all cancel. No catalogue and no model enters.

A measurement operator that is only approximately linear over the increment
recovers an approximately exact difference; the tolerance below is the validated
one (7 / 20 / 12 %, widened by 3 sigma of the difference and capped at twice the
floor, because the increment is known exactly and a 50 % miss cannot be excused
by a large claimed uncertainty).
"""

from __future__ import annotations

import math

import pytest

from conftest import INJECTED, INJECTION_FLOORS

CAP_FACTOR = 2.0


def _pair(reference, iid):
    base = "R" + iid[1:]
    return reference[iid], reference[base], base


def test_every_injection_increment_is_recovered(reference, injection_truth):
    failures, checked = [], 0
    for iid in INJECTED:
        copy, base, base_id = _pair(reference, iid)
        injected = injection_truth[iid]["injected"]
        for line, floor in INJECTION_FLOORS.items():
            if line not in injected:
                continue
            need = injected[line]["flux"]
            if need <= 0:
                continue
            measured = copy.lines[line].flux - base.lines[line].flux
            sig = math.hypot(copy.lines[line].flux_err, base.lines[line].flux_err)
            tol = min(max(floor * need, 3.0 * sig), CAP_FACTOR * floor * need)
            checked += 1
            if abs(measured - need) > tol:
                failures.append(
                    f"{iid}/{line}: recovered {measured:.2f} of {need:.2f} "
                    f"({100 * (measured / need - 1):+.1f} %), tol {tol:.2f}")
    assert checked >= 12, f"only {checked} comparisons were graded"
    assert not failures, "injection differential failures:\n  " + "\n  ".join(failures)


@pytest.mark.characterisation
def test_bright_halpha_increments_are_recovered_to_a_few_percent(reference, injection_truth):
    """Halpha on the bright injections: ~1-3 %. This is the headline validation."""
    worst = 0.0
    for iid in INJECTED:
        copy, base, _ = _pair(reference, iid)
        need = injection_truth[iid]["injected"]["Halpha"]["flux"]
        measured = copy.lines["Halpha"].flux - base.lines["Halpha"].flux
        worst = max(worst, abs(measured / need - 1.0))
    assert worst < 0.05, f"worst Halpha increment miss {100 * worst:.1f} %"


def test_the_copy_differs_from_its_base(reference):
    """An implementation that returns the base for the copy is not measuring.

    This is the anti-degeneracy check: if the added line were ignored, every
    difference would be identically zero and the differential test above would
    have nothing to grade.
    """
    for iid in INJECTED:
        copy, base, base_id = _pair(reference, iid)
        assert copy.lines["Halpha"].flux != pytest.approx(base.lines["Halpha"].flux, rel=1e-6), \
            f"{iid} is indistinguishable from {base_id}: the increment was lost"
