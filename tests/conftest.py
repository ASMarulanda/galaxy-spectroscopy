"""Shared fixtures.

The expensive thing in this suite is fitting: one full 20-case pass costs about a
minute on one core. It is therefore computed **once per session** and shared, so
the accuracy, redshift, uncertainty and regression tests all judge the *same*
numbers. That is not only faster — it means a disagreement between two tests can
never be caused by two different fits.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import gspec

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CASES = DATA / "cases"

#: Relative accuracy floors, per line, in the units of the reported flux.
#: These are the *acceptance* floors of the validation, not accuracy claims:
#: Hbeta's 20 % is an admission threshold (the v0.2 limitations document, L1).
FLOORS = {"Halpha": 0.07, "Hbeta": 0.20, "OIII_5007": 0.06, "NII_6583": 0.10}

#: Looser floors for a *difference* of two measurements: the continuum error is
#: carried twice, which is why [O III] is 12 % there and 6 % for an absolute flux.
INJECTION_FLOORS = {"Halpha": 0.07, "Hbeta": 0.20, "OIII_5007": 0.12}

SYNTHETIC = tuple(f"S{i:02d}" for i in range(1, 11))
REAL = tuple(f"R{i:02d}" for i in range(1, 7))
INJECTED = tuple(f"I{i:02d}" for i in range(1, 5))
ALL_CASES = tuple(sorted(SYNTHETIC + REAL + INJECTED))


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA


@pytest.fixture(scope="session")
def synthetic_truth() -> dict:
    return json.loads((DATA / "truth" / "truth_synthetic.json").read_text())


@pytest.fixture(scope="session")
def injection_truth() -> dict:
    return json.loads((DATA / "truth" / "truth_injections.json").read_text())


@pytest.fixture(scope="session")
def golden_rows() -> dict:
    """The validated baseline, as {(case_id, line): {field: float|str}}."""
    import csv

    path = ROOT / "validation" / "results" / "golden_reference.csv"
    with open(path, newline="") as fh:
        return {(r["case_id"], r["line"]): r for r in csv.DictReader(fh)}


@pytest.fixture(scope="session")
def cases() -> dict[str, gspec.Case]:
    return {cid: gspec.read_case(CASES / f"{cid}.fits") for cid in ALL_CASES}


@pytest.fixture(scope="session")
def reference(cases) -> dict[str, gspec.SpectrumResult]:
    """The reference model (m1a) on all 20 cases, fitted once."""
    return {cid: gspec.measure(case) for cid, case in cases.items()}
