"""Named model configurations: the reference, one alternative, one control.

These are *configuration objects*, not forked code. Every variant runs the same
fitter through the same modules; only the window table, the polynomial order and
the absorption switch change. That is what makes an ablation honest: if a
variant used its own copy of the fitter, "the models differ" would be
indistinguishable from "the code differs".

===========  ==========================================================  =====
config       what it is                                                  status
===========  ==========================================================  =====
``m1a``      cubic continuum, 3-line Balmer absorption, BIC selection     default
``m1b``      quadratic continuum, windows 8 A wider, Hgamma weight 1.05,
             no BIC sweep                                                alternative
``m2``       no stellar-absorption term at all                            negative control
===========  ==========================================================  =====

``m1b`` is **not** an independent implementation. It shares this package's
fitter and differs in three constants, so agreement between m1a and m1b shows
that the answer is insensitive to those constants — nothing more. Genuine
independence is a separately written program, published in v0.2.

``m2`` exists to be rejected. It removes the absorption term and must fail the
accuracy test; if it passes, the test is not measuring what it claims.
"""

from __future__ import annotations

from .results import ModelConfig

#: The graded reference: what the baseline numbers were produced with.
M1A = ModelConfig(variant="m1a", ndeg=3, use_absorption=True)

#: Alternative configuration: different continuum order, wider windows, Hgamma
#: weight 1.05, no model selection.
M1B = ModelConfig(
    variant="m1b", ndeg=2, use_absorption=True,
    pad=8.0, hb_window=(4235.0, 5105.0), gamma_weight=1.05,
)

#: Negative control: no stellar absorption.
M2 = ModelConfig(variant="m2", ndeg=3, use_absorption=False)

VARIANTS = {"m1a": M1A, "m1b": M1B, "m2": M2}

DEFAULT = "m1a"


def get(name: str) -> ModelConfig:
    """Look up a named configuration, returning a fresh copy."""
    try:
        base = VARIANTS[name]
    except KeyError:
        raise KeyError(f"unknown variant {name!r}; choose from {sorted(VARIANTS)}") from None
    return ModelConfig(
        variant=base.variant, ndeg=base.ndeg, allow_broad_b=base.allow_broad_b,
        allow_broad_a=base.allow_broad_a, use_absorption=base.use_absorption,
        pad=base.pad, hb_window=base.hb_window, gamma_weight=base.gamma_weight,
        absorbers_hb=base.absorbers_hb,
    )
