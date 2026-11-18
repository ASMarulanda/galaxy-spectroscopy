"""The bound-aware flag policy, pinned on its own.

Two slow tests apply this policy to real fits: ``test_core_identities`` to a refit
of S08, and ``test_regression`` to the 80 rows of the golden baseline. Both need the
fitter, both take minutes, and neither can make the classifier misbehave on demand —
so the branch that must *fail* would otherwise never be exercised.

This file pins the policy itself against hand-built results, in milliseconds. The
rule is the one in ``conftest``:

* a flag is assigned from the reported flux and its uncertainty, so two fits can
  only be required to agree on it when the line's complex is free of active bounds
  in **both** fits;
* where a bound is pinned, each flag must instead be the one that its own flux,
  uncertainty and broad amplitude imply — so the flag may differ between fits, but
  never from its own evidence;
* the exclusion may not swallow the comparison: a coverage floor is enforced.

Nothing here is a portability claim. It says which comparisons are meaningful, not
that the numbers reproduce everywhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from conftest import (assert_flag_coverage, assert_flag_policy, bound_in_complex,
                      flags_comparable, reported_flag)
from gspec.flags import classify

LINE = "Halpha"


@dataclass
class _Line:
    """The four fields the policy and the classifier read from a line."""

    flux: float
    flux_err: float
    flag: str
    broad_amplitude: float = 0.0


@dataclass
class _Result:
    """The two things the policy reads from a fit: its lines and its bounds."""

    lines: dict = field(default_factory=dict)
    at_bounds: list = field(default_factory=list)


def _consistent(flux: float, flux_err: float, broad: float = 0.0) -> _Line:
    """A line whose flag is what its own numbers imply."""
    return _Line(flux, flux_err, classify(flux, flux_err, broad), broad)


def _fit(lines: dict, at_bounds=()) -> _Result:
    return _Result(lines=dict(lines), at_bounds=list(at_bounds))


# --------------------------------------------------------- unbound complexes
def test_unbound_complex_with_matching_flags_is_compared_exactly():
    """The strong claim still runs where it can hold."""
    a = _fit({LINE: _consistent(30.0, 1.0)})
    b = _fit({LINE: _consistent(30.5, 1.02)})

    assert flags_comparable(a, b, LINE)
    assert assert_flag_policy(a, b, LINE) == "exact"


def test_unbound_complex_with_mismatched_flags_fails():
    """A flag that moves with no bound to explain it is a defect, and is caught."""
    a = _fit({LINE: _consistent(30.0, 1.0)})          # ok
    b = _fit({LINE: _consistent(2.0, 1.0)})           # nondetection

    assert a.lines[LINE].flag != b.lines[LINE].flag
    with pytest.raises(AssertionError, match="no active bound"):
        assert_flag_policy(a, b, LINE)


# ------------------------------------------------------ bound-sensitive rows
def test_bound_sensitive_complex_tolerates_differing_flags_when_each_is_consistent():
    """S08 Hbeta on the CI runner: different flags, each true to its own numbers."""
    a = _fit({LINE: _consistent(30.0, 1.0)}, at_bounds=["amp_Ha_b"])
    b = _fit({LINE: _consistent(2.0, 1.0)}, at_bounds=["amp_Ha_b"])

    assert not flags_comparable(a, b, LINE)
    assert a.lines[LINE].flag != b.lines[LINE].flag, "the case under test needs a flip"
    assert assert_flag_policy(a, b, LINE) == "consistency"


def test_bound_sensitive_complex_with_an_internally_inconsistent_flag_fails():
    """Exemption from agreement is not exemption from its own evidence."""
    a = _fit({LINE: _consistent(30.0, 1.0)}, at_bounds=["amp_Ha_b"])
    b = _fit({LINE: _consistent(30.0, 1.0)}, at_bounds=["amp_Ha_b"])
    b.lines[LINE].flag = "nondetection"          # S/N is 30; the flag contradicts it

    with pytest.raises(AssertionError, match="imply"):
        assert_flag_policy(a, b, LINE)


def test_a_bound_in_either_fit_is_enough_to_withhold_exact_comparison():
    """The rule is 'both free', not 'the reference is free'."""
    a = _fit({LINE: _consistent(30.0, 1.0)})
    b = _fit({LINE: _consistent(30.0, 1.0)}, at_bounds=["amp_Ha_b"])

    assert bound_in_complex(a, LINE) is False
    assert bound_in_complex(b, LINE) is True
    assert not flags_comparable(a, b, LINE)


def test_a_bound_in_another_complex_does_not_disable_the_comparison():
    """The mapping is per complex, so an unrelated pin cannot widen the exemption."""
    other = _fit({"Hbeta": _consistent(30.0, 1.0)}, at_bounds=["amp_Ha_b"])
    twin = _fit({"Hbeta": _consistent(30.0, 1.0)}, at_bounds=["amp_Ha_b"])

    assert bound_in_complex(other, "Hbeta") is False
    assert flags_comparable(other, twin, "Hbeta")


# --------------------------------------------------------------- the floor
def test_coverage_floor_fails_when_the_exclusion_grows():
    with pytest.raises(AssertionError, match="bound exclusion"):
        assert_flag_coverage(compared=12, total=80, minimum=40)


def test_coverage_floor_passes_at_and_above_the_bound():
    assert_flag_coverage(compared=40, total=80, minimum=40)   # no raise
    assert_flag_coverage(compared=44, total=80, minimum=40)   # no raise


# ---------------------------------------------- fidelity to the real rule
def test_the_consistency_check_uses_the_published_classifier_not_a_sn_shortcut():
    """``reported_flag`` must be the production rule, all three criteria of it.

    `gspec.flags.classify` tests finiteness and a positive uncertainty, then the
    3-sigma convention, then the blend branch. The blend branch cannot be reached —
    `nondetection` is tested first, and |F|/sigma >= 3 implies sigma < 0.9|F| — which
    is exactly why a bound-pinned row can be judged on its own numbers, and why the
    broad amplitude cannot rescue a flag that its flux contradicts. Pinned
    independently by `test_unresolved_blend_is_unreachable_by_construction`.
    """
    for flux, flux_err in ((30.0, 1.0), (2.0, 1.0), (3.01, 1.0), (2.99, 1.0)):
        assert classify(flux, flux_err, 5.0) == classify(flux, flux_err, 0.0), (
            "a broad amplitude changed the flag where it is reachable")

    assert classify(float("nan"), 1.0) == "failed"
    assert classify(1.0, 0.0) == "failed"
    assert classify(1.0, -1.0) == "failed"

    # and the helper reproduces the rule through a result's own fields
    for line in (_Line(30.0, 1.0, "ok"), _Line(2.0, 1.0, "nondetection"),
                 _Line(1.0, 0.0, "failed")):
        assert reported_flag(_fit({LINE: line}), LINE) == line.flag
