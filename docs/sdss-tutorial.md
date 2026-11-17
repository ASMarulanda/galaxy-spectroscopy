# Retrieving and measuring your own SDSS spectrum

A walkthrough of the V0.3 workflow: an object's identifiers, the archive, the FITS
file, the conversion, the measurement, and the provenance that ties them together.

Everything below was run to write this page. The numbers quoted are the numbers
that came out.

---

## 1. What you need

```
python3 -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
```

That installs `gspec` and its four pinned dependencies (`numpy 1.26.4`,
`scipy 1.11.4`, `astropy 7.0.0`, `matplotlib 3.6.3`). Network access is needed for
the download; nothing else in the repository is online.

## 2. What an SDSS optical spectrum is

A galaxy's optical spectrum is flux density against wavelength. Between 3800 and
9200 Å it is dominated by starlight — a pseudo-continuum crossed by absorption
features — with emission lines superimposed where gas is ionised, by young stars or
by an accretion flow, and recombines.

The SDSS Legacy spectrograph recorded 640 spectra at once on a drilled plug-plate,
and the survey's own pipeline reduced each fibre to a calibrated, sky-subtracted
one-dimensional spectrum. The file you will download is one such spectrum, from the
SDSS-I/II survey, reduced with `run2d=26` and served by the DR18 Science Archive
Server.

## 3. Identifying a spectrum

Three integers identify it, and they are unique together:

| identifier | what it is | range used by the tool |
|---|---|---|
| `plate` | which plug-plate was drilled | 1–9999 |
| `mjd` | which night it was observed | 40000–69999 |
| `fiber` | which hole on that plate | 1–1000 |

They build this path:

```
https://data.sdss.org/sas/dr18/prior-surveys/sdss2-dr8-sdss/spectro/redux/26/
        spectra/{plate:04d}/spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits
```

Plate and fibre are zero-padded to four digits in both the directory and the file
name. This matters: `spectra/294/…` returns **404**, `spectra/0294/…` returns 200.
(An earlier version of `DATA.md` showed the directory unpadded; that was wrong, and
V0.3 corrects it.)

The tool validates all three as integers inside the ranges above *before* it builds
a string, so nothing that is not an integer can reach a URL.

## 4. What is in the FITS file

```
 0  PRIMARY   header only: COEFF0, COEFF1, DC-FLAG, BUNIT, VACUUM, MJD, FIBERID, EXPTIME
 1  COADD     3828 rows x 8 columns: flux, loglam, ivar, and_mask, or_mask, wdisp, sky, model
 2  SPECOBJ   1 row x 126 columns: the catalogue entry (class, redshift, ...)
 3  SPZLINE   29 rows x 19 columns: the pipeline's own line fits
 4–9          plug-map exposures (B1/R1 ...)
```

**The wavelength grid is logarithmic.** With `DC-FLAG=1`,

$$\log_{10}\lambda_i = \mathrm{COEFF0} + (\mathrm{offset}+i)\,\mathrm{COEFF1},
\qquad \mathrm{COEFF1} = 10^{-4},$$

so a pixel covers 0.92 Å at 4000 Å and 2.03 Å at 8800 Å. The converter checks that
the stored `loglam` column actually follows this law, that the steps agree with
`COEFF1` to within 1e-5, and that the first pixel sits an integer number of pixels
from `COEFF0` (here, 15). A grid that does not satisfy this is refused.

**Flux is a density**, in units of `1E-17 erg/cm^2/s/Ang`. An emission-line flux is
the integral over the line profile, not a column.

**Noise is inverse variance**: `ivar = 1/σ²`, and `ivar = 0` means *no
information*, not *zero uncertainty*.

**`wdisp` is the instrumental resolution in pixels.** The file does not say so —
the column carries no `TUNIT` — and the distinction matters, because the pipeline
convolves its model with this width. It can be settled from the data: for roughly
constant resolving power, the width in Ångström grows with wavelength while the
width in pixels falls as 1/Δλ. Measured on this file, `wdisp` falls from 1.258 at
3984 Å to 0.696 at 8807 Å, and `wdisp × Δλ` stays within 1.10–1.41 Å, giving

| λ (Å) | 3984 | 5189 | 6757 | 8807 |
|---|---|---|---|---|
| R = λ / (2.355 σ) | 1466 | 1997 | 2168 | 2650 |

which is the published shape of the Legacy spectrograph's resolution curve.
Interpreted as Ångström the same numbers would give R ≈ 3200, inconsistent with the
instrument. The converter asserts the pixel interpretation and refuses a file whose
`wdisp` grows with wavelength.

**`SPECOBJ` and `SPZLINE` are the answers.** They carry the survey's own redshift
and line fluxes for this object. They are dropped, and the provenance record says
so.

The object used throughout this tutorial, SDSS 0294-51986-0036, is the source of
the case the repository calls **R06**; `DATA.md` §1.1 states the mapping and why it
is now published. That is an *identification* — it tells you which spectrum the case
file came from. It is not a validation of the fitted numbers, and §9 below is the
list of reasons a fitted Hβ on this object means very little.

## 5. The workflow, executed

```bash
python tools/download_sdss.py --plate 294 --mjd 51986 --fiber 36 --out data/raw
```

```
SDSS-0294-51986-0036: https://data.sdss.org/sas/dr18/prior-surveys/sdss2-dr8-sdss/spectro/redux/26/spectra/0294/spec-0294-51986-0036.fits
  fetched: 604800 bytes  sha256 12e069c2cceb5345...
  case: 3f35ce89277a08d9...  R(Halpha)~2097.8

1 fetched, 0 failed; manifest in data/raw/provenance/fetch_manifest.json
```

Then measure:

```bash
gspec measure --cases data/raw/cases --out data/raw/measurements.csv
```

```
SDSS-0294-51986-0036: z=0.08120  with-abs, ndeg=3, no broad Ha
wrote data/raw/measurements.csv (4 rows)
```

or both steps at once:

```bash
python tools/download_sdss.py --list objects.csv --out data/raw --measure
```

`objects.csv` needs columns named `plate`, `mjd` and `fiber` (the catalogue's
`fiberID` is accepted too), in any order.

### What came out

```
case            SDSS-0294-51986-0036
redshift        z = 0.081202   (+/- 0.000035)
model selected  with-abs, ndeg=3, no broad Ha

line         rest (A)        flux      err     S/N  flag
Halpha        6564.61       9.878    5.031     2.0  nondetection
Hbeta         4862.68      56.292  449.299     0.1  nondetection
OIII_5007     5008.24       1.628    2.561     0.6  nondetection
NII_6583      6585.27       8.851    1.798     4.9  ok
```

This object is faint, and three of four lines are non-detections. That is the
pipeline reporting honestly rather than failing: the flux beside a non-detection is
a limit, not a measurement. Only [N II] is recovered with a signal worth quoting.

## 6. What `gspec` extracts, and how

**The redshift.** One redshift is fitted jointly to both line windows, so the
positions of Hα, Hβ, [O III] and [N II] must all be consistent with a single
$z = \Delta\lambda/\lambda_{\rm rest}$. The quoted uncertainty is statistical and
comes from the covariance; it says nothing about systematics such as the wavelength
calibration of the survey.

**Line fluxes.** Each line is modelled as a Gaussian whose integral is the flux,

$$F = A\,\sigma_{\rm obs}\sqrt{2\pi}\,\Delta\lambda_{\rm pixel}\cdot\frac{1}{\Delta\lambda_{\rm pixel}}
= A\,\sigma_{\rm obs}\sqrt{2\pi},$$

with the instrumental width deconvolved from the observed one:

$$\sigma_{\rm intrinsic}^2 = \sigma_{\rm obs}^2 - \sigma_{\rm instr}^2,
\qquad \sigma_{\rm instr} = \texttt{wdisp}\times\Delta\lambda.$$

**Uncertainties.** The covariance of the model parameters is propagated through the
same arithmetic, so `flux_err` includes the correlation between amplitude and width
rather than treating them as independent.

**The continuum** is a polynomial in each window, and a stellar-absorption term is
included when BIC says it earns its parameters. This is why `model selected`
appears in the output: it is part of the result, not an implementation detail.

## 7. What the conversion does

`tools/fits_to_case.py` writes the six columns a case needs into a new FITS with a
case header. It copies; it does not transform:

| | |
|---|---|
| interpolation, resampling, smoothing | **none** — the pixels downloaded are the pixels fitted |
| renormalisation | **none** — the flux scale is the survey's |
| columns dropped | `sky`, `model` |
| extensions dropped | `SPECOBJ`, `SPZLINE`, the plug maps |
| conversions applied | float64 → float32, which is how the archive stores them already |

The case header carries `CASEID`, `BUNIT`, `VACUUM`, `PROVEN`, `NWAVE` plus
informational `PLATE`/`MJD`/`FIBERID`/`ORIGIN`/`SRCNAME` that `gspec` never reads.

If a file cannot be converted safely, the tool says which requirement failed and
exits 4 without writing anything. It does not bend the data to fit the pipeline,
and it does not modify the fitter to accept new input.

## 8. Reading the flags

| flag | means | how to read the flux |
|---|---|---|
| `ok` | detected, plausible width, S/N above the floor | a measurement |
| `nondetection` | the fit could not separate the line from the continuum | an upper limit |
| `unresolved_blend` | detected, but the width sits on the model's bound, so neighbouring lines cannot be separated | neither; treat the pair together |

Non-detections are not errors. A faint galaxy has faint lines, and a method that
reported `ok` for all of them would be lying about its own precision.

## 9. What remains limited

1. **One spectrum is not a sample.** Nothing here supports a statement about a
   population, and the pipeline fails on individual objects often enough that
   `docs/limitations.md` lists the cases.
2. **Hβ is the weakest of the four lines.** Across four fitting configurations its
   median error moved from −14 % to +20 % depending only on whether the
   absorption term was kept, with individual cases off by more than 90 %
   (`docs/limitations.md` L2).
3. **The absorption term is phenomenological.** The model includes it; the project
   has not shown that it separates stellar from nebular light (`docs/methods.md`
   §4.2).
4. **Uncertainties are fragile on some spectra.** Where parameters sit on their
   bounds, the formal errors depend on how the covariance handles them; the tested
   alternative made calibration worse (`docs/limitations.md` L4).
5. **The flux scale is the survey's**: no aperture correction, no extinction
   correction.

## 10. Provenance

Every download writes two records under `data/raw/provenance/`:

* `<name>.provenance.json` — the archive URL, the HTTP response headers the server
  sent, the bytes written, the SHA-256, and the identifiers **you** supplied;
* `<name>.conversion.json` — the input hash, the output hash, the columns and
  extensions dropped, and the verification report (pixel count, grid checks,
  inverse-variance fraction, resolution).

Plus one `fetch_manifest.json` per run, listing every object with its status, so a
CSV can always be traced to bytes:

```bash
sha256sum data/raw/raw/spec-0294-51986-0036.fits      # 12e069c2…
sha256sum data/raw/cases/spec-0294-51986-0036.fits    # 3f35ce89…
```

The fitter is deterministic, so the same case produces the same CSV byte for byte.

## 11. When it goes wrong

| symptom | cause | what to do |
|---|---|---|
| `HTTP 404` | one identifier is wrong, or the object was never observed | check the triple against a catalogue; 404 is not retried |
| `gave up after 3 attempts` | the archive is unreachable | try again later; `.part` files resume automatically |
| `incomplete download: N bytes of M` | the transfer broke | re-run: the `.part` is kept and the download resumes |
| `SHA-256 mismatch` | corrupt transfer or wrong object | the file is deleted; re-run, or check `--sha256` |
| `not compatible with the pipeline` | the file violates the instrument contract; the reasons are printed | read them — a different reduction or a different survey is usually the cause |
| `exists; pass --force to replace it` | an earlier download is present | the tool will not overwrite silently |

Exit codes: `0` success, `1` usage, `2` verification failed, `3` network, `4`
instrument contract violated.
