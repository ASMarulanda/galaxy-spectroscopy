"""Does the shipped generator still produce the shipped data?

The synthetic cases are the only ones with *exact* truth, and that truth is only
meaningful if the spectra came from the code in `tools/`. This test regenerates
all ten from the recorded seed and compares every column of every FITS file bit
for bit.

It is a reproducibility test of the data, not of the fitter: if it fails, the
generator or the shim has changed, and every accuracy figure in
the v0.2 validation report refers to data this repository can no longer produce.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

BASE_SEED = 20261106           # the seed recorded in DATA.md §4
COLUMNS = ("flux", "loglam", "ivar", "and_mask", "or_mask", "wdisp")


@pytest.fixture(scope="module")
def generator():
    import generate_synthetic
    return generate_synthetic


@pytest.mark.parametrize("index", range(10))
def test_synthetic_case_reproduces_bit_for_bit(generator, data_dir, index):
    spec = generator.CASES[index]
    cid = spec[0]
    wave, loglam, flux, ivar, am, om, wd, truth, _model = generator.build_case(
        cid, *spec[1:], BASE_SEED)

    produced = {
        "flux": np.asarray(flux, np.float32),
        "loglam": np.asarray(loglam, np.float32),
        "ivar": np.asarray(ivar, np.float32),
        "and_mask": np.asarray(am, np.int32),
        "or_mask": np.asarray(om, np.int32),
        "wdisp": np.asarray(wd, np.float32),
    }
    with fits.open(data_dir / "cases" / f"{cid}.fits") as hdul:
        stored = {k: np.asarray(hdul[1].data[k]) for k in COLUMNS}

    for name in COLUMNS:
        np.testing.assert_array_equal(produced[name], stored[name],
                                      err_msg=f"{cid}: column {name} differs")


def test_generator_self_test_still_holds(generator):
    """The generator's own quadrature check: every declared flux is re-integrated."""
    for spec in generator.CASES:
        cid = spec[0]
        wave, _l, _f, _i, _m, _o, _w, truth, _model = generator.build_case(
            cid, *spec[1:], BASE_SEED)
        rel = generator.self_test(truth, wave)
        assert rel < 1e-6, f"{cid}: generator self-test rel error {rel:.2e}"
