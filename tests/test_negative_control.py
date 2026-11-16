"""The negative control: a wrong model must be rejected by the same tests.

`m2` removes the stellar-absorption term and changes nothing else. If the suite
accepts it, the suite is not measuring what it claims.

What `m2` actually fails — measured, not assumed
------------------------------------------------
On the ten synthetic spectra `m2` is *not* clearly worse: its Hbeta median error
is -11.1 % against the reference's -14.4 %. This is an uncomfortable but
important result, and it is why the repository makes no claim that the absorption
model improves Hbeta (`docs/methods.md` §4.2, limitations L3).

Where `m2` does fail is structural: it loses the I04 injection entirely (the fit
collapses to z = 0.0195 instead of 0.0826, a 18 910 km/s error), and it misses
the archival redshift on two real spectra. So the control is tested against the
property it demonstrably breaks — the differential — rather than against a
synthetic-accuracy claim it happens to survive.
"""

from __future__ import annotations

import math

import pytest

import gspec
from gspec import variants

from conftest import INJECTED, INJECTION_FLOORS

CAP_FACTOR = 2.0


@pytest.fixture(scope="module")
def control(cases) -> dict[str, gspec.SpectrumResult]:
    """The m2 (no-absorption) results on the injection cases and their bases."""
    model = variants.get("m2")
    want = list(INJECTED) + ["R" + iid[1:] for iid in INJECTED]
    return {cid: gspec.measure(cases[cid], model=model) for cid in want}


def test_control_loses_at_least_one_injection_increment(control, injection_truth, reference):
    """The control must fail where the reference passes."""
    broken = []
    for iid in INJECTED:
        base = "R" + iid[1:]
        injected = injection_truth[iid]["injected"]
        for line, floor in INJECTION_FLOORS.items():
            if line not in injected:
                continue
            need = injected[line]["flux"]
            if need <= 0:
                continue
            mc, mb = control[iid].lines[line], control[base].lines[line]
            measured = mc.flux - mb.flux
            sig = math.hypot(mc.flux_err, mb.flux_err)
            tol = min(max(floor * need, 3.0 * sig), CAP_FACTOR * floor * need)
            if abs(measured - need) > tol:
                broken.append(iid)
    assert broken, (
        "the no-absorption control passed every injection comparison; "
        "the negative control no longer discriminates and the suite is suspect"
    )


def test_control_collapses_the_i04_redshift(control, reference):
    """The documented failure mode: the seed fails and the whole case moves.

    I04 is the case whose injected increment is smallest relative to the line,
    and it is the case the control gets wrong by ~19 000 km/s.
    """
    dv = abs(control["I04"].redshift - reference["I04"].redshift) * 299792.458
    assert dv > 1000.0, f"I04 redshift moved only {dv:.0f} km/s; the control behaved like the reference"


def test_control_reports_a_different_model(control):
    for r in control.values():
        assert "no-abs" in r.model_selected, "the control claims to be fitting absorption"
        assert r.model.variant == "m2"
