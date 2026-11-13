"""gspec — redshift estimation and emission-line analysis of optical galaxy spectra.

Quick start
-----------
    import gspec
    case = gspec.read_case("data/cases/S01.fits")
    result = gspec.measure(case)
    result.redshift, result.lines["Halpha"].flux

The public API is defined in :mod:`gspec.results` (the data structures),
:mod:`gspec.io` (reading a case) and :mod:`gspec.pipeline` (the measurement).
Everything else is an implementation module whose contract is documented in
``docs/methods.md``.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [
    "__version__",
]
