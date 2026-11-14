"""The shipped alternative configuration (`m1b`) and the control (`m2`).

Both are exposed in `gspec.variants` and both are described in the v0.2 methods document,
so both need a test — otherwise a configuration that nobody exercises can drift
away from its documentation without anything noticing.

What `m1b` is *not*: an independent implementation. It runs the same fitter with
three constants changed (quadratic continuum, windows 8 Å wider, Hγ weight 1.05)
and no BIC sweep. The measurement here shows that this is a *different model*,
not a rerun of the same one — its Hβ differs from the reference by a median of
28 %, and in the opposite direction.
"""

from __future__ import annotations

import json

import pytest

import gspec
from gspec import variants
from gspec.pipeline import measure

LINES = ("Halpha", "Hbeta", "OIII_5007", "NII_6583")
BAND = {"Halpha": 0.07, "Hbeta": 0.20, "OIII_5007": 0.06, "NII_6583": 0.10}
SYNTHETIC = ("S01", "S03", "S07")


@pytest.fixture(scope="module")
def m1b(data_dir) -> dict:
    model = variants.get("m1b")
    return {cid: measure(gspec.read_case(data_dir / "cases" / f"{cid}.fits"), model=model)
            for cid in SYNTHETIC}


def test_m1b_produces_a_valid_measurement(m1b):
    for cid, r in m1b.items():
        assert r.model.variant == "m1b"
        assert r.model.ndeg == 2
        assert r.model.use_absorption is True
        assert 0.5 < r.chi2dof < 3.0
        assert len(r.lines) == 4
        for line, m in r.lines.items():
            assert m.flux > 0 and m.flux_err > 0
            assert m.flag in ("ok", "nondetection", "unresolved_blend")


def test_m1b_stays_inside_the_acceptance_band(data_dir, m1b):
    truth = json.loads((data_dir / "truth" / "truth_synthetic.json").read_text())
    for cid, r in m1b.items():
        for line in LINES:
            true = truth[cid]["lines"][line]["flux"]
            m = r.lines[line]
            assert abs(m.flux - true) <= max(BAND[line] * true, 3.0 * m.flux_err), \
                f"{cid}/{line}: {m.flux:.2f} vs {true:.2f} outside the band"


def test_m1b_is_a_different_model_not_a_rerun(data_dir, reference, m1b):
    """Its Hβ must differ from the reference — otherwise the variant is inert.

    This is the check that keeps `variants.py` honest: if changing the continuum
    order, the windows and the absorption weights changed nothing, the "three
    constants" claim in the docs would be false.
    """
    gaps = []
    for cid, r in m1b.items():
        a, b = reference[cid].lines["Hbeta"].flux, r.lines["Hbeta"].flux
        gaps.append(abs(a - b) / b)
    assert max(gaps) > 0.10, f"m1b and m1a agree on Hbeta to {max(gaps):.1%}; the variant is inert"


def test_the_control_has_no_absorption_term():
    """`m2` is the negative control, and this is the property that defines it."""
    assert variants.get("m2").use_absorption is False
    assert variants.get("m1a").use_absorption is True
    # and the lookup returns a copy, so a caller cannot corrupt the module state
    a = variants.get("m1a")
    a.ndeg = 99
    assert variants.get("m1a").ndeg == 3


def test_unknown_variant_is_rejected():
    with pytest.raises(KeyError):
        variants.get("m3")
