"""Reading a case file and deciding which pixels carry information.

This module owns exactly two ideas and nothing else:

1. the **FITS contract** — which HDU, which columns, which units, what precision;
2. the **mask policy** — which pixels may constrain a fit.

It contains no physics. Anything that touches a wavelength, a line or a
parameter belongs downstream.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import constants as K


@dataclass
class Case:
    """One calibrated spectrum, in float64, on its native pixel grid.

    Attributes
    ----------
    wave
        Observed wavelength in Angstrom (vacuum), reconstructed as
        ``10 ** loglam``. Never re-binned: every array stays on the pixels the
        data were sampled on.
    flux
        Flux **density**, 1e-17 erg s^-1 cm^-2 A^-1. Not an integrated flux.
    ivar
        ``1 / sigma^2`` per pixel, in the same units as ``flux**2``.
        ``ivar == 0`` means "no information", not "zero variance".
    wdisp
        Instrumental sigma in **pixels**. Must be converted with the *local*
        delta-lambda; a global value would bias widths by ~10 % across the band.
    """

    case_id: str
    wave: np.ndarray
    flux: np.ndarray
    ivar: np.ndarray
    wdisp: np.ndarray
    loglam: np.ndarray
    and_mask: np.ndarray
    or_mask: np.ndarray
    bunit: str = K.BUNIT
    provenance: str = ""
    path: str = ""

    @cached_property
    def dlam(self) -> np.ndarray:
        """Local pixel width dlambda/dpixel (Angstrom), from the native grid."""
        return np.gradient(self.wave)

    def __len__(self) -> int:
        return int(self.wave.size)


def read_case(path: str | Path) -> Case:
    """Read one case file. Returns float64 arrays on the native pixel grid.

    Deliberately dumb: it validates nothing, cleans nothing, and does not assume
    the array is sorted. The only statement it makes is the schema.
    """
    path = Path(path)
    with fits.open(path) as hdul:
        data = hdul[1].data
        prim = hdul[0].header
        out = Case(
            case_id=str(prim.get("CASEID", "")).strip(),
            bunit=str(prim.get("BUNIT", "")).strip(),
            provenance=str(prim.get("PROVEN", "")).strip(),
            flux=np.asarray(data["flux"], float),
            loglam=np.asarray(data["loglam"], float),
            ivar=np.asarray(data["ivar"], float),
            and_mask=np.asarray(data["and_mask"], np.int32),
            or_mask=np.asarray(data["or_mask"], np.int32),
            wdisp=np.asarray(data["wdisp"], float),
            wave=10.0 ** np.asarray(data["loglam"], float),
            path=str(path),
        )
    return out


def clean(case: Case, fatal_bits: int = K.FATAL_MASK_BITS) -> np.ndarray:
    """Boolean mask of pixels that may constrain a fit.

    A pixel is usable when it is finite, has positive inverse variance, and
    carries no fatal mask bit.

    Two consequences that are easy to get wrong:

    * ``ivar = 0`` is "no information", so including it would assert infinite
      weight rather than zero weight.
    * ``wdisp`` can be zero on rejected pixels. Because this mask is applied
      *before* any width is computed, a zero LSF never reaches the model and
      cannot produce a zero-width (delta-function) component.
    """
    ok = np.isfinite(case.flux) & np.isfinite(case.ivar) & (case.ivar > 0)
    ok &= (case.and_mask & fatal_bits) == 0
    return ok
