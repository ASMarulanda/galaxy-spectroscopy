"""C3: are the reported uncertainties honest?

An uncertainty contract is only meaningful if it is tested against the errors
actually made:

    pull = (value - true) / sigma          -> should scatter with unit RMS

Three things are checked, and each fails in a different way:

* **pull RMS in [0.5, 2.5]** — optimistic sigma gives RMS >> 1, defensive sigma
  gives RMS << 1.
* **at most 20 % of pulls beyond 5 sigma** — a heavy tail means a few rows are
  wildly wrong while the mean looks acceptable.
* **no row reports sigma above 1.5x the true flux** — without this ceiling the
  tolerance ``max(floor, 3 sigma)`` could be bought outright by inflating sigma.
  This is the rule that closed a measured exploit in the validation history.

The sample is synthetic rows plus the twelve differential comparisons, which is
what the validated calibration used (n = 52).
"""

from __future__ import annotations

import math
import statistics as st

import pytest

from conftest import FLOORS, INJECTED, INJECTION_FLOORS, SYNTHETIC

SIGMA_PER_FLUX = 1.5
PULL_RMS_LO, PULL_RMS_HI = 0.5, 2.5
HARD_Z, HARD_FRAC = 5.0, 0.20
MIN_POINTS = 46


def _pulls(reference, synthetic_truth, injection_truth):
    """Pull for every non-failed synthetic row and every differential comparison."""
    pulls = []
    for cid in SYNTHETIC:
        truth = synthetic_truth[cid]["lines"]
        for line in FLOORS:
            m = reference[cid].lines[line]
            if m.flag == "failed" or m.flux_err <= 0:
                continue
            pulls.append((m.flux - truth[line]["flux"]) / m.flux_err)
    for iid in INJECTED:
        base = "R" + iid[1:]
        injected = injection_truth[iid]["injected"]
        for line in INJECTION_FLOORS:
            if line not in injected:
                continue
            need = injected[line]["flux"]
            mc, mb = reference[iid].lines[line], reference[base].lines[line]
            sig = math.hypot(mc.flux_err, mb.flux_err)
            if sig <= 0:
                continue
            pulls.append((mc.flux - mb.flux - need) / sig)
    return pulls


def test_pull_distribution_is_calibrated(reference, synthetic_truth, injection_truth):
    pulls = _pulls(reference, synthetic_truth, injection_truth)
    assert len(pulls) >= MIN_POINTS, f"only {len(pulls)} graded points"
    rms = math.sqrt(sum(p * p for p in pulls) / len(pulls))
    assert PULL_RMS_LO <= rms <= PULL_RMS_HI, f"pull RMS {rms:.2f} outside [{PULL_RMS_LO}, {PULL_RMS_HI}]"
    beyond = sum(abs(p) > HARD_Z for p in pulls) / len(pulls)
    assert beyond <= HARD_FRAC, f"{100 * beyond:.0f} % of pulls beyond {HARD_Z} sigma"


@pytest.mark.characterisation
def test_measured_calibration(reference, synthetic_truth, injection_truth):
    """Pull RMS 1.08 over 52 points, median sigma/F 5.2 % (validated values)."""
    pulls = _pulls(reference, synthetic_truth, injection_truth)
    rms = math.sqrt(sum(p * p for p in pulls) / len(pulls))
    assert rms == pytest.approx(1.08, abs=0.25)
    assert st.median(abs(p) for p in pulls) == pytest.approx(0.85, abs=0.4)


def test_no_row_reports_a_defensive_uncertainty(reference, synthetic_truth):
    """sigma <= 1.5 x the true flux: the ceiling that closed the exploit.

    Both validated implementations stay below 1.18, and only on faint
    non-detections.
    """
    offenders = []
    for cid in SYNTHETIC:
        for line in FLOORS:
            true = synthetic_truth[cid]["lines"][line]["flux"]
            m = reference[cid].lines[line]
            if m.flux_err > SIGMA_PER_FLUX * true:
                offenders.append(f"{cid}/{line}: sigma={m.flux_err:.2f} > {SIGMA_PER_FLUX}*{true:.2f}")
    assert not offenders, "\n  ".join(offenders)


def test_well_measured_lines_have_small_relative_uncertainties(reference):
    """For rows flagged `ok`, sigma/F must be small — otherwise the flag lies.

    Deliberately restricted to `ok`: on a `nondetection` the flux is near zero and
    sigma/F is meaningless (it reaches 18 on the faintest rows, which is exactly
    why those rows are not flagged ok).
    """
    ok_ratios = []
    for cid in SYNTHETIC:
        for line in FLOORS:
            m = reference[cid].lines[line]
            if m.flag == "ok":
                ok_ratios.append(m.flux_err / m.flux)
    assert len(ok_ratios) >= 25, f"only {len(ok_ratios)} ok rows on synthetic data"
    assert min(ok_ratios) < 0.02, "no well-measured line has a sub-2 % uncertainty"
    assert max(ok_ratios) < 0.50, "an `ok` row carries a 50 %+ relative uncertainty"


def test_flags_are_consistent_with_the_reported_numbers(reference):
    """`ok` <=> |F|/sigma >= 3, `nondetection` <=> |F|/sigma < 3.

    The flag must describe the numbers in the same row, never a hidden threshold.
    """
    for cid, r in reference.items():
        for line, m in r.lines.items():
            q = abs(m.flux) / m.flux_err
            if m.flag == "ok":
                assert q >= 3.0, f"{cid}/{line}: ok but |F|/sigma={q:.2f}"
            elif m.flag == "nondetection":
                assert q < 3.0, f"{cid}/{line}: nondetection but |F|/sigma={q:.2f}"
            elif m.flag == "unresolved_blend":
                assert m.broad_amplitude > 0 and m.flux_err > 0.9 * abs(m.flux)
            elif m.flag == "failed":
                raise AssertionError(f"{cid}/{line} failed on a well-behaved spectrum")
            else:
                raise AssertionError(f"unknown flag {m.flag!r}")
