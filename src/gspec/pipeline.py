"""The measurement pipeline: one call from a case to a structured result.

This module is the only place the stages are wired together, in the order they
must run:

    read a case  ->  find redshift candidates  ->  fit jointly  ->  simplify by
    BIC  ->  propagate uncertainty  ->  classify

Everything it calls is testable on its own; nothing downstream of
:func:`measure` needs to know how the stages are arranged.
"""

from __future__ import annotations

from pathlib import Path

from . import variants
from .fit import fit_case
from .io import Case, read_case
from .results import ModelConfig, SpectrumResult
from .select import fit_with_model_selection


def measure(
    case: Case,
    model: ModelConfig | None = None,
    n_starts: int = 2,
    model_select: bool = True,
    z_known: float | None = None,
    keep_jacobian: bool = False,
) -> SpectrumResult:
    """Measure one spectrum.

    Parameters
    ----------
    model
        A named configuration from :mod:`gspec.variants`, or ``None`` for the
        reference (``m1a``).
    n_starts
        Number of redshift seeds carried into the full fit. Two is the default;
        one is enough when the seed is supplied explicitly.
    model_select
        Run the BIC parsimony sweep. Ignored for variants other than ``m1a``.
    z_known
        Fit at this redshift seed and skip the search. Used by experiments that
        compare models under identical conditions.
    keep_jacobian
        Retain the Jacobian on the result. Needed only by covariance
        experiments; the default keeps results light.
    """
    cfg = model if model is not None else variants.get(variants.DEFAULT)
    return fit_with_model_selection(
        case, n_starts=n_starts, z_known=z_known,
        model_select=model_select, config=cfg,
    )


def measure_fixed_model(
    case: Case,
    model: ModelConfig,
    z_known: float | None = None,
    keep_jacobian: bool = False,
) -> SpectrumResult:
    """Fit one configuration with no model selection at all.

    This is the entry point for controlled experiments: the model is exactly what
    was asked for, so a difference between two results cannot be attributed to
    the selection step.
    """
    return fit_case(case, z_init=z_known, config=model, keep_jacobian=keep_jacobian)


def measure_file(path: str | Path, **kwargs) -> SpectrumResult:
    return measure(read_case(path), **kwargs)
