"""Which lines exist, in which window, and how their amplitudes are tied.

A *complex* is a rest-frame wavelength window containing several emission
components that must be fitted together. There are two:

* the **Halpha complex** (6420-6790 A rest): [N II] 6548 (tied), Halpha narrow,
  Halpha broad, [N II] 6583, [S II] 6716, [S II] 6731;
* the **Hbeta complex** (4250-5090 A rest): Hgamma, Hbeta narrow, Hbeta broad,
  [O III] 4959 (tied), [O III] 5007.

Why two windows and not one: the complexes are ~1700 A apart at rest, and a
*shared* redshift across both is what breaks the redshift-width degeneracy that
a single window suffers from. Each window keeps its own polynomial continuum,
which absorbs the real shape difference without forcing a global model.

Why the ties: the [N II] and [O III] doublet ratios are set by transition
probabilities, not by the galaxy, so tying them removes two free parameters per
complex and stabilises the deblending.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import constants as K


@dataclass(frozen=True)
class Emission:
    """One emission component.

    Attributes
    ----------
    key
        Identifier used to build the parameter name ``amp_<key>``.
    rest
        Vacuum rest wavelength, Angstrom.
    sigma
        Name of the shared parameter holding ``log(sigma_v)`` for this component.
    tie
        ``(other_key, factor)`` when this component's amplitude is not free but
        is ``factor * amplitude(other_key)``.
    line
        Which reported line this component contributes to. ``None`` means it is
        modelled but never reported (it prevents it from contaminating the
        graded lines).
    """

    key: str
    rest: float
    sigma: str
    tie: tuple[str, float] | None = None
    line: str | None = None


@dataclass(frozen=True)
class Complex:
    """A rest-frame window and the components inside it."""

    name: str
    lo_rest: float
    hi_rest: float
    emis: tuple[Emission, ...]
    absorbers: tuple[tuple[float, float, float], ...] = field(default=())

    @property
    def width(self) -> float:
        return self.hi_rest - self.lo_rest

    @property
    def centre(self) -> float:
        return 0.5 * (self.lo_rest + self.hi_rest)


#: Broad G-band absorption: (rest A, depth relative to Hbeta, width multiplier).
GBAND = (4305.00, 0.60, 3.0)

#: The measurement side's Balmer absorption approximation. Deliberately NOT the
#: generator's 12-line table: the two are allowed to disagree, so a fit cannot
#: pass by assuming the model that produced the data.
HB_ABSORBERS_3LINE = ((4341.69, 0.95, 1.0), GBAND, (4862.68, 1.00, 1.0))


def alpha_complex(pad: float = 0.0) -> Complex:
    """The Halpha complex. ``pad`` widens the window by that many Angstrom."""
    return Complex(
        name="Ha",
        lo_rest=6420.0 - pad,
        hi_rest=6790.0 + pad,
        emis=(
            Emission("NII6548", K.NII_6548, "narrowA",
                     tie=("NII6583", 1.0 / K.NII_RATIO_6583_6548)),
            Emission("Ha_n", K.HALPHA, "narrowA", line="Halpha"),
            Emission("Ha_b", K.HALPHA, "broad", line="Halpha"),
            Emission("NII6583", K.NII_6583, "narrowA", line="NII_6583"),
            Emission("SII6716", K.SII_6716, "narrowA"),
            Emission("SII6731", K.SII_6731, "narrowA"),
        ),
        absorbers=(),
    )


def beta_complex(
    window: tuple[float, float] = (4250.0, 5090.0),
    gamma_weight: float = 0.95,
    absorbers: tuple[tuple[float, float, float], ...] | None = None,
) -> Complex:
    """The Hbeta complex. ``window``, ``gamma_weight`` and ``absorbers`` are the knobs.

    ``absorbers`` exists so a controlled experiment can swap the 3-line
    Balmer approximation for a richer basis without forking the model. It is
    never set by the default pipeline.
    """
    gband = (GBAND[0], GBAND[1], GBAND[2])
    if absorbers is None:
        absorbers = (
            (4341.69, gamma_weight, 1.0),
            gband,
            (4862.68, 1.00, 1.0),
        )
    return Complex(
        name="Hb",
        lo_rest=window[0],
        hi_rest=window[1],
        emis=(
            Emission("Hg_em", K.HGAMMA, "narrowB"),
            Emission("Hb_n", K.HBETA, "narrowB", line="Hbeta"),
            Emission("Hb_b", K.HBETA, "broad", line="Hbeta"),
            Emission("OIII4959", K.OIII_4959, "narrowO",
                     tie=("OIII5007", 1.0 / K.OIII_RATIO_5007_4959)),
            Emission("OIII5007", K.OIII_5007, "narrowO", line="OIII_5007"),
        ),
        absorbers=absorbers,
    )


def default_complexes(
    pad: float = 0.0,
    hb_window: tuple[float, float] = (4250.0, 5090.0),
    gamma_weight: float = 0.95,
) -> list[Complex]:
    """The reference configuration (m1a): cubic continuum, 3-line absorption."""
    return [alpha_complex(pad), beta_complex(hb_window, gamma_weight)]


#: Which reported line each emission key contributes to. Components absent from
#: this map are modelled but never reported, so a blend cannot leak a sibling's
#: flux into a graded line silently.
LINE_OF_KEY = {
    e.key: e.line
    for cx in default_complexes()
    for e in cx.emis
    if e.line is not None
}
