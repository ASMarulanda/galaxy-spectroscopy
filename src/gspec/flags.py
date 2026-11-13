"""Quality flags.

The flag describes the measurement, not the object. It is defined by a single
consistency rule — the reported flux against the reported uncertainty — so that a
row can never be flagged ``ok`` while its own numbers say it is a non-detection.

``ok``
    ``|F| / sigma_F >= 3``. Reliable measurement.
``nondetection``
    ``|F| / sigma_F < 3``. The ``flux`` column still carries the best estimate,
    not an upper limit.
``unresolved_blend``
    A broad component is present *and* the error exceeds 90 % of the flux: the
    line cannot be separated from its neighbour, so ``flux`` is a total or a
    best estimate.

    **This branch is unreachable in the current implementation.** ``nondetection``
    is tested first, so reaching this line requires ``|F|/σ ≥ 3``, which implies
    ``σ ≤ |F|/3 < 0.9|F|`` — the two conditions cannot hold together. The value is
    kept because it is part of the published output contract and because the
    ordering is inherited from the original implementation; changing it would
    change outputs. Pinned by
    ``tests/test_core_identities.py::test_unresolved_blend_is_unreachable_by_construction``.
``failed``
    No finite estimate, or a non-positive uncertainty.

The 3-sigma threshold is a convention, not a detection theory. On the 20 spectra
of this repository 9 of 80 rows sit within a factor 1.5 of that boundary, so for
~11 % of rows the flag is a knife-edge — see docs/limitations.md L6.
"""

from __future__ import annotations

import math

from . import constants as K

FLAGS = ("ok", "nondetection", "unresolved_blend", "failed")


def classify(
    flux: float,
    flux_err: float,
    broad_amp: float = 0.0,
    detection_sigma: float = K.DETECTION_SIGMA,
    blend_fraction: float = K.BLEND_ERR_FRACTION,
) -> str:
    """Assign the flag from the reported numbers alone."""
    if not math.isfinite(flux) or not math.isfinite(flux_err) or flux_err <= 0:
        return "failed"
    if flux / flux_err < detection_sigma:
        return "nondetection"
    if broad_amp > 0.0 and flux_err > blend_fraction * abs(flux):
        return "unresolved_blend"
    return "ok"


def significance(flux: float, flux_err: float) -> float:
    """``|F| / sigma_F``, guarded against a zero uncertainty."""
    if not math.isfinite(flux_err) or flux_err <= 0:
        return float("nan")
    return abs(flux) / flux_err
