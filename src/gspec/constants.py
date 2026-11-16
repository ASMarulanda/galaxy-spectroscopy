"""Physical constants, rest wavelengths, atomic ratios and the mask policy.

Nothing here is computed. It is the vocabulary the rest of the package uses, and
it is the only module whose values a reader may want to change for a different
survey — hence the explicit statements of which convention each number follows.

Conventions
-----------
* Wavelengths are **vacuum** Angstrom. SDSS spectra are on a vacuum scale.
* Flux density is in 1e-17 erg s^-1 cm^-2 A^-1; integrated fluxes are in
  1e-17 erg s^-1 cm^-2. ``UNITS`` and ``BUNIT`` record both.
"""

from __future__ import annotations

#: Speed of light, km/s (CODATA). Used to convert velocity dispersions to widths.
C_KMS = 299792.458

# --------------------------------------------------------------------- lines
# Vacuum rest wavelengths (Angstrom).
HALPHA = 6564.61
HBETA = 4862.68
OIII_5007 = 5008.24
NII_6583 = 6585.27

NII_6548 = 6549.86
OIII_4959 = 4960.29
SII_6716 = 6718.29
SII_6731 = 6732.67
OI_6300 = 6300.30
OI_6364 = 6364.31
HGAMMA = 4341.69
HDELTA = 4102.89
HEPSILON = 3971.20
NEIII_3869 = 3869.86
OII_3727A = 3727.09
OII_3727B = 3729.88

#: The four lines this package reports, with canonical names.
TARGET_LINES = {
    "Halpha": HALPHA,
    "Hbeta": HBETA,
    "OIII_5007": OIII_5007,
    "NII_6583": NII_6583,
}
LINE_NAMES = tuple(TARGET_LINES)

#: Intensity ratios fixed by atomic physics, applied as exact ties between the
#: components of a doublet. Tying them removes two free parameters per complex.
#: Treating them as exact is a modelling choice; real gas with dust and shocks
#: varies by tens of percent (`docs/limitations.md`, L5).
NII_RATIO_6583_6548 = 2.94      # I(6583)/I(6548)
OIII_RATIO_5007_4959 = 2.98     # I(5007)/I(4959)

# ---------------------------------------------------------------- SDSS masks
#: SDSS SPZLINE-style mask bits (the subset this package ever inspects).
MASK_NOPLUG = 1 << 0
MASK_MANYBADCOLUMNS = 1 << 4
MASK_NEARBADPIXEL = 1 << 16
MASK_LOWFLAT = 1 << 17
MASK_FULLREJECT = 1 << 18
MASK_PARTIALREJECT = 1 << 19
MASK_SCATTEREDLIGHT = 1 << 20
MASK_NOSKY = 1 << 22
MASK_BRIGHTSKY = 1 << 23
MASK_NODATA = 1 << 24
MASK_COMBINEREJ = 1 << 25

#: Bits that make a pixel unusable. NEARBADPIXEL (16), LOWFLAT (17) and
#: BRIGHTSKY (23) are deliberately *not* fatal: they flag a caution, not a
#: rejection, and treating them as fatal would discard usable pixels.
FATAL_MASK_BITS = MASK_FULLREJECT | MASK_NOSKY | MASK_NODATA

# -------------------------------------------------------------------- units
UNITS = "1e-17 erg s-1 cm-2"
BUNIT = "1E-17 erg/cm^2/s/Ang"

# --------------------------------------------------------------- fit policy
#: Relative covariance scaling is capped here. On real spectra chi2/dof is
#: driven by model mismatch (absorption lines, sky residuals), not by photon
#: noise; uncapped scaling once produced 35 % errors on 40-sigma lines.
#: Measured on this dataset it binds on 2 of 8 cases (`docs/methods.md` §6).
COV_SCALE_CAP = 1.5

#: A line is "detected" when |F| / sigma_F exceeds this. A convention, not a
#: detection theory; it is the contract the flags are defined against.
DETECTION_SIGMA = 3.0

#: Above this error-to-flux ratio, a detected broad component means the line
#: cannot be separated from its neighbour.
BLEND_ERR_FRACTION = 0.9

#: Parameter names shared by every model configuration.
NAME_Z = "z"
NAME_NARROW_A = "narrowA"
NAME_NARROW_B = "narrowB"
NAME_NARROW_O = "narrowO"
NAME_BROAD = "broad"
NAME_ABS_SCALE = "abs_scale"
NAME_ABS_RATIO = "abs_ratio"
NAME_GBAND = "gband"

#: Bounds on velocity dispersions, in km/s (log-space internally).
SIGMA_MIN_KMS = 25.0
SIGMA_MAX_KMS = 600.0
BROAD_MIN_KMS = 250.0
BROAD_MAX_KMS = 1800.0

#: PHYSICAL PRIOR: stellar absorption must be at least this much broader than
#: the ionised gas. Stars fill the galaxy; H II regions are kinematically cold.
#: Without this tie, Hbeta absorption is degenerate with a local continuum and
#: the fit was biased by +10..+80 % (`docs/methods.md` §4).
ABS_RATIO_MIN = 1.15
ABS_RATIO_MAX = 6.0
