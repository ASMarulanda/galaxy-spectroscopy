"""Import shim for the verbatim provenance tools.

``tools/generate_synthetic.py`` and ``tools/make_injections.py`` are the programs
that *made the data*. They are copied here byte for byte from the original
project, with exactly one edit: the two lines that put a sibling ``common.py`` on
``sys.path`` now import this module instead.

That is deliberate. Rewriting them to the package's API would mean the generator
that produced `data/cases/` is no longer the generator sitting in `tools/` — and
the reproducibility claim ("S01–S10 come from this seed and this code") would rest
on a program that never ran. So the algorithm is untouched and only the names it
resolves are supplied here.

`tests/test_generator.py` regenerates the ten synthetic spectra through this shim
and checks them against the shipped FITS files bit for bit; if the shim were
wrong, that test would fail rather than silently produce different data.
"""

from __future__ import annotations

import json

import numpy as np
from astropy.io import fits

# Re-export the physical vocabulary the tools expect to find on `common`.
from gspec.constants import (  # noqa: F401
    BUNIT,
    C_KMS,
    HALPHA,
    HBETA,
    HDELTA,
    HEPSILON,
    HGAMMA,
    MASK_BRIGHTSKY,
    MASK_FULLREJECT,
    MASK_NEARBADPIXEL,
    NEIII_3869,
    NII_6548,
    NII_6583,
    NII_RATIO_6583_6548,
    OI_6300,
    OI_6364,
    OIII_4959,
    OIII_5007,
    OIII_RATIO_5007_4959,
    SII_6716,
    SII_6731,
)

#: The original module spelled the [O II] doublet with lower-case letters; kept so
#: the copied code needs no edit beyond the import line.
OII_3727a = 3727.09
OII_3727b = 3729.88


def sigma_lsf_kms(wave):
    """Empirical SDSS LSF sigma used by the generator: 55 km/s at 5000 A, slope 0.35.

    This is the *generator's* model, not the fitter's: the fitter reads the
    per-pixel ``wdisp`` column that this function was used to write. That
    asymmetry is intentional (the two sides are allowed to disagree).
    """
    return 55.0 * (np.asarray(wave, float) / 5000.0) ** 0.35


def write_case(path, case_id, loglam, flux, ivar, and_mask, or_mask, wdisp,
               provenance, extra_header=None):
    """Write one case file: PRIMARY (metadata + BUNIT) + COADD (six columns).

    Identical to the writer that produced `data/cases/`, including the float32
    storage of every column (DATA.md §5).
    """
    prim = fits.PrimaryHDU()
    h = prim.header
    h["CASEID"] = (case_id, "opaque case identifier")
    h["BUNIT"] = (BUNIT, "flux density unit of COADD.flux")
    h["VACUUM"] = (True, "wavelengths are vacuum")
    h["PROVEN"] = (provenance, "synthetic | sdss-dr18-redacted")
    h["NWAVE"] = (len(loglam), "number of spectral pixels")
    if extra_header:
        for k, v in extra_header.items():
            h[k] = v

    cols = [
        fits.Column(name="flux", format="E", unit=BUNIT, array=np.asarray(flux, np.float32)),
        fits.Column(name="loglam", format="E", unit="log10(Angstrom)",
                    array=np.asarray(loglam, np.float32)),
        fits.Column(name="ivar", format="E", unit="(1e-17 erg/s/cm2/A)^-2",
                    array=np.asarray(ivar, np.float32)),
        fits.Column(name="and_mask", format="J", array=np.asarray(and_mask, np.int32)),
        fits.Column(name="or_mask", format="J", array=np.asarray(or_mask, np.int32)),
        fits.Column(name="wdisp", format="E", unit="pixel", array=np.asarray(wdisp, np.float32)),
    ]
    coadd = fits.BinTableHDU.from_columns(cols, name="COADD")
    fits.HDUList([prim, coadd]).writeto(path, overwrite=True)


def read_case(path):
    """Read one case file, mirroring the writer above (used by make_injections)."""
    with fits.open(path) as f:
        d = f[1].data
        prim = f[0].header
        out = dict(
            case_id=str(prim.get("CASEID", "")).strip(),
            bunit=str(prim.get("BUNIT", "")).strip(),
            provenance=str(prim.get("PROVEN", "")).strip(),
            flux=np.asarray(d["flux"], float),
            loglam=np.asarray(d["loglam"], float),
            ivar=np.asarray(d["ivar"], float),
            and_mask=np.asarray(d["and_mask"], np.int32),
            or_mask=np.asarray(d["or_mask"], np.int32),
            wdisp=np.asarray(d["wdisp"], float),
        )
    out["wave"] = 10.0 ** out["loglam"]
    return out


def load_truth(path):
    with open(path) as fh:
        return json.load(fh)


def dump_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)


def gauss_int_flux(amp, sigma_ang):
    """Integrated flux of a Gaussian: A * sigma * sqrt(2 pi).

    Kept for the tools, but note it duplicates :func:`gspec.fluxes.integrated_flux`
    — there is exactly one place in the package that owns this physics, and the
    tools are the one documented exception because they must stay verbatim.
    """
    import math
    return float(amp) * float(sigma_ang) * math.sqrt(2.0 * math.pi)
