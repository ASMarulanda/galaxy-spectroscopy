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
#: Hbeta's 20 % is an admission threshold (`docs/limitations.md`, L1).
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


# ------------------------------------------------------------------ V0.3 tools
# The retrieval and conversion tools live outside the package (tools/), because
# they are not part of the measured scientific API. Tests import them by path.

import sys as _sys

_sys.path.insert(0, str(ROOT / "tools"))


@pytest.fixture
def make_raw(tmp_path):
    """Build an SDSS-like archive file, by default the one R06 was made from.

    The arrays come from the shipped, redacted case ``R06.fits``, which DATA.md §3
    records as carrying the archive COADD columns unchanged. Rebuilding a raw file
    from it therefore gives a fixture whose conversion must reproduce R06 exactly —
    that identity is what the end-to-end test asserts. Each keyword argument breaks
    exactly one thing, so a test can prove the converter refuses it.
    """
    from astropy.io import fits
    import numpy as np

    def _make(name="spec-0294-51986-0036.fits", *, drop=(), loglam=None, bunit=None,
              vacuum=True, ivar=None, dflag=1, coeff0=None, coeff1=None,
              wdisp=None, n_rows=None):
        ref = fits.open(CASES / "R06.fits")[1].data
        n = n_rows or len(ref)
        data = {c: np.asarray(ref[c])[:n] for c in
                ("flux", "loglam", "ivar", "and_mask", "or_mask", "wdisp")}
        if loglam is not None:
            data["loglam"] = np.asarray(loglam, dtype=np.float32)
        if ivar is not None:
            data["ivar"] = np.asarray(ivar, dtype=np.float32)
        if wdisp is not None:
            data["wdisp"] = np.asarray(wdisp, dtype=np.float32)

        prim = fits.PrimaryHDU()
        h = prim.header
        h["BUNIT"] = "1E-17 erg/cm^2/s/Ang" if bunit is None else bunit
        if vacuum is not None:
            h["VACUUM"] = vacuum
        h["COEFF0"] = float(data["loglam"][0]) if coeff0 is None else coeff0
        h["COEFF1"] = 1e-4 if coeff1 is None else coeff1
        h["DC-FLAG"] = dflag
        h["MJD"] = 51986
        h["FIBERID"] = 36
        h["RUN2D"] = 26
        h["TELESCOP"] = "SDSS 2.5-M"

        cols = [fits.Column(name=c, format="J" if "mask" in c else "E",
                            array=data[c]) for c in data if c not in drop]
        if "sky" not in drop:
            cols.append(fits.Column(name="sky", format="E",
                                    array=np.zeros(n, dtype=np.float32)))
            cols.append(fits.Column(name="model", format="E",
                                    array=np.asarray(ref["flux"])[:n] * 0.5))
        coadd = fits.BinTableHDU.from_columns(cols, name="COADD")
        # a second extension that carries an "answer"; the converter must drop it
        spec = fits.BinTableHDU.from_columns(
            [fits.Column(name="z", format="E", array=np.array([0.081202], dtype=np.float32)),
             fits.Column(name="class", format="6A", array=np.array(["GALAXY"]))],
            name="SPECOBJ")
        path = tmp_path / name
        fits.HDUList([prim, coadd, spec]).writeto(path, overwrite=True)
        return path

    return _make
