"""Structured results: what the pipeline returns, diagnostics included.

The diagnostics are not debugging aids. On this dataset the *selected model* is a
scientific result: the BIC sweep rejects the stellar-absorption term on all ten
synthetic spectra and keeps it on five of six real ones (the v0.2 methods document §4.2).
A reader who only sees ``flux`` cannot tell that apart from a fit that used
absorption everywhere, so the model and the bound state travel with the numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ComponentMeasurement:
    """One fitted emission component (a narrow or broad part of a line)."""

    key: str
    rest: float
    amplitude: float
    sigma_v_kms: float
    sigma_int: float
    sigma_lsf: float
    sigma_obs: float
    flux: float


@dataclass
class LineMeasurement:
    """One reported line: the sum over its components."""

    line: str
    flux: float
    flux_err: float
    flag: str
    sigma_kms: float = float("nan")
    broad_amplitude: float = 0.0
    components: list[ComponentMeasurement] = field(default_factory=list)

    @property
    def significance(self) -> float:
        """``|F| / sigma_F``; ``nan`` when the uncertainty is not positive."""
        if self.flux_err <= 0 or not np.isfinite(self.flux_err):
            return float("nan")
        return abs(self.flux) / self.flux_err


@dataclass
class ModelConfig:
    """Which model is being fitted.

    ``variant`` selects a named configuration (see :mod:`gspec.variants`); the
    remaining fields are the individual knobs that a controlled experiment needs
    to move one at a time.
    """

    variant: str = "m1a"
    ndeg: int = 3
    allow_broad_b: bool = True
    allow_broad_a: bool = True
    use_absorption: bool = True
    pad: float = 0.0
    hb_window: tuple[float, float] = (4250.0, 5090.0)
    gamma_weight: float = 0.95
    absorbers_hb: tuple[tuple[float, float, float], ...] | None = None

    def label(self) -> str:
        """Human-readable identity of the model, e.g. 'no-abs, ndeg=2, no broad Hb'."""
        parts = ["with-abs" if self.use_absorption else "no-abs", f"ndeg={self.ndeg}"]
        if not self.allow_broad_b:
            parts.append("no broad Hb")
        if not self.allow_broad_a:
            parts.append("no broad Ha")
        if self.absorbers_hb is not None:
            parts.append(f"abs-basis={len(self.absorbers_hb)}")
        return ", ".join(parts)


@dataclass
class SpectrumResult:
    """Everything the pipeline knows about one spectrum."""

    case_id: str
    redshift: float
    redshift_err: float
    lines: dict[str, LineMeasurement]
    chi2: float = float("nan")
    dof: int = 0
    n_par: int = 0
    covariance_scale: float = float("nan")
    model: ModelConfig = field(default_factory=ModelConfig)
    model_selected: str = ""
    n_starts_used: int = 0
    parameters: dict[str, float] = field(default_factory=dict)
    at_bounds: list[str] = field(default_factory=list)
    error: str | None = None
    jacobian: np.ndarray | None = None

    @property
    def chi2dof(self) -> float:
        return self.chi2 / self.dof if self.dof else float("nan")

    @property
    def ok_fraction(self) -> float:
        if not self.lines:
            return 0.0
        return sum(m.flag == "ok" for m in self.lines.values()) / len(self.lines)

    def row(self, line: str) -> dict:
        """One CSV-shaped row for ``line``."""
        m = self.lines[line]
        return dict(
            case_id=self.case_id,
            line=line,
            redshift=self.redshift,
            redshift_err=self.redshift_err,
            flux=m.flux,
            flux_err=m.flux_err,
            flag=m.flag,
        )
