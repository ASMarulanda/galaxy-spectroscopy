"""gspec — redshift estimation and emission-line analysis of optical galaxy spectra.

Quick start
-----------
    import gspec
    case = gspec.read_case("data/cases/S01.fits")
    result = gspec.measure(case)
    result.redshift, result.lines["Halpha"].flux

Public API
----------
:func:`read_case` / :class:`Case`
    Reading one calibrated spectrum; the FITS contract and the mask policy.
:func:`measure` / :func:`measure_fixed_model`
    The pipeline, and its no-selection form for controlled experiments.
:class:`SpectrumResult`, :class:`LineMeasurement`, :class:`ComponentMeasurement`
    Structured results, carrying the fit diagnostics alongside the numbers.
:data:`VERSIONS`
    The exact dependency versions the numbers in this repository were produced
    under. Check these before comparing against a published value.
"""

from __future__ import annotations

__version__ = "0.1.0"

from . import constants, variants  # noqa: E402
from .io import Case, clean, read_case  # noqa: E402
from .pipeline import measure, measure_file, measure_fixed_model  # noqa: E402
from .results import (  # noqa: E402
    ComponentMeasurement,
    LineMeasurement,
    ModelConfig,
    SpectrumResult,
)

#: Dependency versions that define the numerics of this package. A result quoted
#: in the docs was produced under these; another environment may move the last
#: digits, and one specific uncertainty is known to be environment-sensitive
#: (docs/limitations.md L4).
VERSIONS = {
    "python": "3.12.3",
    "numpy": "1.26.4",
    "scipy": "1.11.4",
    "astropy": "7.0.0",
}

__all__ = [
    "__version__",
    "Case",
    "clean",
    "read_case",
    "measure",
    "measure_file",
    "measure_fixed_model",
    "SpectrumResult",
    "LineMeasurement",
    "ComponentMeasurement",
    "ModelConfig",
    "constants",
    "variants",
    "VERSIONS",
]
