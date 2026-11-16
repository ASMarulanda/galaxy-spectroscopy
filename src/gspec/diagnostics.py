"""Line-ratio diagnostics: BPT, Balmer decrement, and what they can support here.

Everything in this module is a *measured* ratio with a propagated uncertainty.
The only inferred physical quantity is the dust extinction, and it is returned
together with the assumptions that produced it, because those assumptions are
what a reader needs in order to discount the number.

Three things this module deliberately does not do:

* it does not classify a galaxy. With six objects, four of which have no
  measurable Hβ (see ``docs/figures.md`` §5), a BPT position exists for two of
  them and a classification for none;
* it does not invent a value for a non-detection. :class:`Ratio` carries
  ``detected=False`` and a reason, and the plotting code is expected to honour it;
* it does not assume a reddening law silently. `ccm89_k` is validated against
  published values in ``tests/test_diagnostics.py`` before any E(B−V) is quoted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .results import LineMeasurement, SpectrumResult

#: Case-B recombination value of Hα/Hβ for a T ≈ 10⁴ K, n_e ≈ 10²–10⁴ cm⁻³ gas
#: (Osterbrock & Ferland 2006, Table 4.2; the often-quoted 2.86).
HALPHA_HBETA_CASE_B = 2.86

#: Cardelli, Clayton & Mathis (1989) total-to-selective extinction.
R_V_CCM89 = 3.1

#: Wavelengths used by the diagnostics, vacuum Å.
LAMBDA_HA = 6564.61
LAMBDA_HB = 4862.68

DETECTION_SNR = 3.0


# --------------------------------------------------------------------------- #
# Ratios
# --------------------------------------------------------------------------- #
@dataclass
class Ratio:
    """A line ratio with a first-order uncertainty.

    ``detected`` is False when either line fails the 3σ criterion (or the
    denominator is non-positive, which makes the ratio meaningless). The caller
    must not plot the ``value`` of an undetected ratio as if it were a
    measurement.
    """

    name: str
    value: float = float("nan")
    error: float = float("nan")
    detected: bool = False
    reason: str = ""
    numerator: str = ""
    denominator: str = ""

    @property
    def significance(self) -> float:
        if not math.isfinite(self.error) or self.error <= 0:
            return float("nan")
        return abs(self.value) / self.error

    def as_dict(self) -> dict:
        return dict(name=self.name, value=self.value, error=self.error,
                    detected=self.detected, reason=self.reason,
                    denominator=self.denominator)


def _usable(m: LineMeasurement, min_snr: float) -> tuple[bool, str]:
    if m.flag == "failed":
        return False, f"{m.line} failed to fit"
    if not math.isfinite(m.flux) or not math.isfinite(m.flux_err) or m.flux_err <= 0:
        return False, f"{m.line} has no finite uncertainty"
    if abs(m.flux) / m.flux_err < min_snr:
        return False, f"{m.line} below {min_snr:g}σ"
    return True, ""


def ratio(num: LineMeasurement, den: LineMeasurement, name: str | None = None,
          min_snr: float = DETECTION_SNR) -> Ratio:
    """``num/den`` with first-order propagation of independent errors.

    σ_R² = (σ_num/den)² + (num·σ_den/den²)² — the standard small-error
    propagation, valid while both relative errors are ≪ 1. It says nothing about
    the covariance between the two lines, which the fit does produce but does not
    report per line, so the quoted σ is a *lower bound* on the true uncertainty
    (the lines share the continuum, so they are positively correlated and the
    real ratio is usually *better* determined than this formula implies).
    """
    name = name or f"{num.line}/{den.line}"
    r = Ratio(name=name, numerator=num.line, denominator=den.line)
    ok_n, why_n = _usable(num, min_snr)
    ok_d, why_d = _usable(den, min_snr)
    if not ok_n or not ok_d:
        r.reason = "; ".join(w for w in (why_n, why_d) if w)
        return r
    if den.flux == 0:
        r.reason = f"{den.line} is zero"
        return r
    value = num.flux / den.flux
    error = math.hypot(num.flux_err / den.flux, value * den.flux_err / den.flux)
    r.value, r.error = value, error
    r.detected = True
    return r


def balmer_decrement(result: SpectrumResult, min_snr: float = DETECTION_SNR) -> Ratio:
    """Hα/Hβ, the classic reddening indicator."""
    return ratio(result.lines["Halpha"], result.lines["Hbeta"],
                 name="Hα/Hβ", min_snr=min_snr)


def nii_ha(result: SpectrumResult, min_snr: float = DETECTION_SNR) -> Ratio:
    """[N II] λ6583/Hα — the abscissa of the [N II] BPT diagram."""
    return ratio(result.lines["NII_6583"], result.lines["Halpha"],
                 name="[N II]/Hα", min_snr=min_snr)


def oiii_hb(result: SpectrumResult, min_snr: float = DETECTION_SNR) -> Ratio:
    """[O III] λ5007/Hβ — the ordinate of every BPT diagram."""
    return ratio(result.lines["OIII_5007"], result.lines["Hbeta"],
                 name="[O III]/Hβ", min_snr=min_snr)


# --------------------------------------------------------------------------- #
# Dust extinction
# --------------------------------------------------------------------------- #
def ccm89_k(wavelength_angstrom: float, R_V: float = R_V_CCM89) -> float:
    """Cardelli, Clayton & Mathis (1989) extinction curve, k(λ) = A_λ/E(B−V).

    Valid for 0.3 ≤ x ≤ 3.3 µm⁻¹, i.e. 3030 Å ≤ λ ≤ 33333 Å — which covers Hα
    and Hβ comfortably and excludes the ultraviolet, where CCM89 has a separate
    polynomial that this project has no use for.

    The two branches are the paper's IR (x < 1.1) and optical (1.1 ≤ x ≤ 3.3)
    fits; the y-polynomial coefficients are Table 3 of the paper.
    """
    x = 1.0e4 / wavelength_angstrom  # µm⁻¹
    if not 0.3 <= x <= 3.3:
        raise ValueError(f"CCM89 branch not implemented at x={x:.3f} µm⁻¹ "
                         f"(λ={wavelength_angstrom:.0f} Å)")
    if x < 1.1:
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
    else:
        y = x - 1.82
        a = (1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4
             + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7)
        b = (1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4
             - 0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7)
    return R_V * a + b


@dataclass
class Extinction:
    """The inferred dust column, with the assumptions that produced it.

    ``valid`` means "a measurable reddening excess was detected". When the
    decrement is at or below the case-B value the honest answer is an *upper
    limit*, which is what ``av_upper`` carries.
    """

    e_bv: float = float("nan")
    e_bv_err: float = float("nan")
    av: float = float("nan")
    av_err: float = float("nan")
    av_upper: float = float("nan")
    valid: bool = False
    reason: str = ""
    assumptions: list[str] = field(default_factory=list)


def extinction_from_balmer(
    dec: Ratio,
    intrinsic: float = HALPHA_HBETA_CASE_B,
    R_V: float = R_V_CCM89,
    min_snr: float = DETECTION_SNR,
) -> Extinction:
    """E(B−V) and A_V from an observed Hα/Hβ, in the foreground-screen limit.

    E(B−V) = 2.5/(k(Hβ) − k(Hα)) · log₁₀[(Hα/Hβ)_obs / (Hα/Hβ)_int]

    Assumptions, all of which must survive scrutiny before the number is quoted:

    1. **Case B recombination**, T ≈ 10⁴ K, n_e ≈ 10²–10⁴ cm⁻³, so the intrinsic
       ratio is 2.86. A harder radiation field or collisional excitation deepens
       Hα and lowers the ratio.
    2. **A foreground screen** of dust (not a mixed geometry), so the reddening
       applies as a pure wavelength-dependent attenuation of the emergent flux.
    3. **The Balmer lines are optically thin and not absorbed by stellar
       photospheres.** For a post-starburst galaxy this is the weakest link: the
       underlying Balmer *absorption* depresses Hβ disproportionately, which
       inflates the apparent decrement exactly as dust does. In this pipeline the
       absorption term is fitted explicitly (``use_absorption``), but only for
       the real spectra where BIC keeps it — and the Hβ measurement is the least
       reliable number the pipeline produces (limitation L2).
    4. **CCM89 with R_V = 3.1.** A different reddening law changes E(B−V) by a
       few per cent, and R_V = 3.1 differs from the 3.1 ± 0.2 typical of
       starbursts to the same order.

    Nothing here is a measurement of dust: it is a *measurement of a line ratio*
    plus an inference that is only as good as assumptions 1–4.
    """
    assumptions = [
        "case-B Hα/Hβ = 2.86 (T ≈ 10⁴ K, n_e ≈ 10²–10⁴ cm⁻³)",
        "foreground dust screen, not a mixed geometry",
        "Balmer lines optically thin; stellar Balmer absorption negligible",
        f"CCM89 reddening law, R_V = {R_V:g}",
    ]
    ex = Extinction(assumptions=assumptions)
    if not dec.detected:
        ex.reason = f"Hα/Hβ not measurable: {dec.reason}"
        return ex

    dk = ccm89_k(LAMBDA_HB, R_V) - ccm89_k(LAMBDA_HA, R_V)
    scale = 2.5 / dk

    if dec.value <= intrinsic:
        # No reddening excess. A negative E(B−V) is not physical, so the only
        # defensible statement is an upper limit from the 1σ top of the ratio.
        ex.av_upper = R_V * scale * math.log10((dec.value + dec.error) / intrinsic)
        ex.reason = ("Hα/Hβ is at or below the case-B value: no measurable "
                     "reddening, only an upper limit")
        return ex

    ex.e_bv = scale * math.log10(dec.value / intrinsic)
    ex.e_bv_err = scale * dec.error / (dec.value * math.log(10.0))
    ex.av = R_V * ex.e_bv
    ex.av_err = R_V * ex.e_bv_err
    ex.valid = True
    return ex


# --------------------------------------------------------------------------- #
# BPT reference curves
# --------------------------------------------------------------------------- #
def kauffmann_line(x: np.ndarray | float) -> np.ndarray:
    """Kauffmann et al. (2003) empirical star-formation ceiling.

    log([O III]/Hβ) = 0.61/(log([N II]/Hα) − 0.05) + 1.3, valid below x = 0.05.
    """
    x = np.asarray(x, dtype=float)
    return 0.61 / (x - 0.05) + 1.3


def kewley_line(x: np.ndarray | float) -> np.ndarray:
    """Kewley et al. (2001) theoretical maximum-starburst line.

    log([O III]/Hβ) = 0.61/(log([N II]/Hα) − 0.47) + 1.19, valid below x = 0.47.
    """
    x = np.asarray(x, dtype=float)
    return 0.61 / (x - 0.47) + 1.19


def bpt_cloud(step: float = 0.002) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluation grids for the two curves, clipped to their domains of validity."""
    x_sf = np.arange(-1.80, 0.049, step)
    x_agn = np.arange(-1.80, 0.469, step)
    return x_sf, kauffmann_line(x_sf), x_agn, kewley_line(x_agn)


def diagnostic_table(results: dict[str, SpectrumResult]) -> list[dict]:
    """One row per object: the measured ratios, the inference, and the reasons."""
    rows = []
    for cid in sorted(results):
        r = results[cid]
        dec = balmer_decrement(r)
        n2 = nii_ha(r)
        o3 = oiii_hb(r)
        ex = extinction_from_balmer(dec)
        rows.append(dict(
            case_id=cid,
            halpha_snr=r.lines["Halpha"].significance,
            hbeta_snr=r.lines["Hbeta"].significance,
            oiii_snr=r.lines["OIII_5007"].significance,
            nii_snr=r.lines["NII_6583"].significance,
            balmer=dec.as_dict(),
            nii_ha=n2.as_dict(),
            oiii_hb=o3.as_dict(),
            bpt_placeable=bool(n2.detected and o3.detected),
            extinction=dict(e_bv=ex.e_bv, e_bv_err=ex.e_bv_err, av=ex.av,
                            av_err=ex.av_err, av_upper=ex.av_upper,
                            valid=ex.valid, reason=ex.reason),
        ))
    return rows
