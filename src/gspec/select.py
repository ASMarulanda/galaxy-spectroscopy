"""BIC model selection over nested models, and the multi-start driver.

A broad component or a full cubic can soak up the flux of a weak line on a low
S/N spectrum. Rather than fixing a complexity in advance, the pipeline fits from
several redshift seeds and then accepts a simpler model whenever

    BIC = chi2 + k * ln(N)

says the data support it. Lower is better, so the extra parameters must buy
more than ``ln(N)`` each in chi-square.

This layer is also a scientific measurement in its own right: on the ten
synthetic spectra it rejects the stellar-absorption term entirely, and on five
of six real spectra it keeps it (`docs/methods.md` §4.2). The report of *which*
model won therefore travels with every result.
"""

from __future__ import annotations

import math

from . import constants as K
from .fit import fit_case
from .io import Case
from .redshift import redshift_peaks
from .results import ModelConfig, SpectrumResult

#: The parsimony sweep. Each entry is a candidate (allow_broad_b, use_abs, ndeg,
#: allow_broad_a). The order is the order they are tried; the first is the
#: richest model that still respects the physical prior.
CANDIDATES: tuple[tuple[bool, bool, int, bool], ...] = (
    (True, False, 3, True),
    (False, True, 3, True),
    (False, False, 3, True),
    (False, False, 2, True),
    (False, False, 1, True),
    (True, True, 3, False),
    (False, False, 3, False),
)


def bic(result: SpectrumResult) -> float:
    """Bayesian information criterion of a fitted result."""
    n = max(result.dof + result.n_par, 1)
    return result.chi2 + result.n_par * math.log(n)


def fit_with_model_selection(
    case: Case,
    n_starts: int = 2,
    z_known: float | None = None,
    model_select: bool = True,
    config: ModelConfig | None = None,
) -> SpectrumResult:
    """Fit from several redshift seeds, then optionally simplify by BIC.

    A fixed seed (``z_known``) skips the sweep entirely, which is what a
    controlled experiment wants when it compares models on identical conditions.
    """
    base = config or ModelConfig()

    if z_known is not None:
        return fit_case(case, z_init=z_known, config=base)

    starts = redshift_peaks(case, n=n_starts) or [0.05]
    best: SpectrumResult | None = None
    for z0 in starts:
        try:
            r = fit_case(case, z_init=z0, config=base)
        except Exception:                                 # noqa: BLE001
            continue
        if r.lines and (best is None or r.chi2 < best.chi2):
            best = r
    if best is None:
        return SpectrumResult(case_id=case.case_id, redshift=float("nan"),
                              redshift_err=float("nan"), lines={},
                              model=base, model_selected=base.label(),
                              error="all starts failed")

    # The sweep runs only for the reference variant: it is a statement about the
    # reference model family, not about a variant that already fixed its shape.
    if model_select and base.variant == "m1a":
        for allow_bb, use_abs, nd, allow_ba in CANDIDATES:
            candidate = ModelConfig(
                variant=base.variant, ndeg=nd, allow_broad_b=allow_bb,
                allow_broad_a=allow_ba, use_absorption=use_abs,
                pad=base.pad, hb_window=base.hb_window,
                gamma_weight=base.gamma_weight, absorbers_hb=base.absorbers_hb,
            )
            try:
                r = fit_case(case, z_init=best.redshift, config=candidate)
            except Exception:                             # noqa: BLE001
                continue
            if r.lines and bic(r) < bic(best):
                best = r

    best.n_starts_used = len(starts)
    return best
