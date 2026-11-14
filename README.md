# galaxy-spectroscopy

**Redshift estimation and emission-line analysis of optical galaxy spectra.**

`gspec` measures, from a flux-calibrated optical spectrum of *unknown* redshift:

* the galaxy **redshift**,
* the continuum-subtracted **integrated fluxes** of Hα λ6565, Hβ λ4863,
  [O III] λ5008 and [N II] λ6585,
* a **1σ uncertainty** on each,
* a **quality flag** — by a joint fit of the Hα and Hβ line complexes at one
  common redshift.

It is a small, dependency-light scientific package: numpy, scipy and astropy, no
black boxes, and every number it reports is accompanied by the model that produced
it.

```
f(λ) = C(λ) · Π_i [1 − d_i·G_i(λ)]  +  Σ_j A_j·G_j(λ)

C(λ) = Σ_k c_k x^k                        one polynomial per line complex
       × multiplicative Balmer absorption, width tied to the gas width
F   = A · sqrt(σ_int² + σ_LSF²) · sqrt(2π)     integral of the observed profile
σ_F² = ∇Fᵀ · Cov(p) · ∇F                  Cov = pinv(JᵀJ)·min(χ²/dof, 1.5)
```

Redshift: a matched filter over eleven optical lines proposes candidates; an NNLS
χ² scan polishes each one; the two complexes are then fitted jointly with a shared
`z`, which is what breaks the redshift–width degeneracy that a single window
suffers from.

## Results on the ten synthetic spectra

Each synthetic spectrum was generated from known analytic parameters, so the truth
is arithmetic — not a fit.

| line | median error | worst case | acceptance floor | inside the floor alone |
|---|---|---|---|---|
| Hα | **−0.70 %** | +41.6 % | 7 % | 8 / 10 |
| Hβ | **−14.4 %** | −94.7 % | 20 % | 7 / 10 |
| [O III] 5007 | **+0.31 %** | −65.5 % | 6 % | 7 / 10 |
| [N II] 6583 | **−1.54 %** | −7.7 % | 10 % | 10 / 10 |

Redshift: median |Δv| **1.7 km/s**, maximum **18.5 km/s** over the ten cases.
Differential recovery on four real spectra with a line of exactly known flux
injected: **12 / 12** comparisons pass, Hα to 1–3 %.

Uncertainties: pull RMS **1.232** over the 40 synthetic line fluxes, where both the
value and its error are independently known.

## What this does not prove

Stated here rather than buried, because it is the most useful thing in the
repository:

* **Hβ is the weak link, and it is not a 14 % measurement.** Median −14.4 %, but
  the sign of that error depends on the fitting configuration: disabling the
  model-selection sweep moves Hβ to +19.6 % on the same spectra. A change of one
  switch is worth more than the whole acceptance band. Do not read the Hβ column
  as an accuracy claim.
* **The stellar-absorption term is not shown to help.** The BIC selection rejects
  it on all ten synthetic spectra, and forcing it on makes Hβ worse. It is kept
  because it is physically motivated and it is retained on five of six real
  spectra, not because it improved a metric.
* **σ is not reproducible across environments** where an amplitude sits on its
  parameter bound: measured, one row moves 19.70 → 51.98 between two machines.
* **Real-data fluxes rest on four injections and a factor-two band.** There is no
  exact flux truth for a real spectrum.
* **n = 20.** Nothing here is a statement about galaxies in general.

A full validation report, the independent-implementation comparison and five
figures ship in **v0.2**.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

The pins in `pyproject.toml` (numpy 1.26.4, scipy 1.11.4, astropy 7.0.0) are part
of the science: every number above was produced under them.

## Use

```python
import gspec

case = gspec.read_case("data/cases/S01.fits")
result = gspec.measure(case)

result.redshift                                  # 0.03999954
result.lines["Halpha"].flux                      # 901.331
result.lines["Halpha"].flux_err                  # 1.3898
result.lines["Halpha"].flag                      # 'ok'
result.model_selected                            # 'no-abs, ndeg=2, no broad Hb'
result.chi2dof, result.n_par, result.at_bounds   # diagnostics travel with the numbers
```

Batch interface:

```bash
gspec measure --cases data/cases --out line_measurements.csv
```

See one case end to end, with every intermediate printed in pipeline order:

```bash
python tools/walkthrough_s01.py            # or: python tools/walkthrough_s01.py S05
```

## The model, briefly

* **Two fixed windows** (rest 6420–6790 Å and 4250–5090 Å), each with its own
  polynomial continuum and its own line list: [N II] 6548/6583, Hα (narrow +
  broad), [S II] 6716/6731 on one side; Hγ, Hβ (narrow + broad), [O III]
  4959/5007 on the other.
* **Fixed doublet ratios** ([N II] 6583/6548 = 2.94, [O III] 5007/4959 = 2.98),
  so the doublets constrain the width instead of floating freely.
* **Instrumental broadening in quadrature**, recomputed per pixel from the log
  wavelength grid and the SDSS `wdisp` column.
* **Bounded least squares** with a BIC sweep over seven candidate models
  (continuum order, broad components, the absorption term).
* **Covariance from the Jacobian**, scale-capped at min(χ²/dof, 1.5) where
  χ²/dof > 1, and σ_F propagated through the integrated-flux function itself.

## Tests

```bash
pytest                     # 83 tests, ~3 min
```

They are written against analytic truth and exact injected increments, not
against the pipeline's own output. Notable ones: the synthetic spectra are
regenerated bit for bit from their documented seed; a deliberately broken control
configuration must fail the accuracy test; and the reported σ must be consistent
with the errors actually made.

## Data

20 spectra, 2.0 MB — 10 synthetic with analytic truth, 6 real SDSS-I/II Legacy
spectra with the catalogue answers removed, 4 injection copies of real spectra
with an exactly known increment.

The real spectra are SDSS DR18 (`run2d=26`), public domain. **SDSS acknowledgment
and the required citations are in `DATA.md` and must accompany any reuse.** Code is
MIT; original text, figures and synthetic data are CC BY 4.0; the SDSS spectra are
neither and keep their own terms.

Provenance, the redaction record and the licensing position: **`DATA.md`**.

## Layout

```
src/gspec/       constants → io → redshift → fit → fluxes → flags → pipeline
data/            20 case FITS + truth manifests + provenance record
tools/           the generator, the injection builder and a guided walkthrough
tests/           83 tests, one per scientific failure mode
validation/      the golden regression baseline
```

## Roadmap

* **v0.1** (this release) — the measurement pipeline, its tests, and the data.
* **v0.2** — validation report, independent-implementation comparison, five
  reproducible figures, and the full limitations document.
* **v0.3** — download a spectrum from the SDSS archive and measure it end to end,
  with integrity verification and recorded provenance.

## Citing

`CITATION.cff`. Any reuse must also acknowledge SDSS as described in `DATA.md`.

## Licence

MIT for the code (`LICENSE`); CC BY 4.0 for original text and synthetic data
(`LICENSE-content`). Third-party material is excluded from both grants — see
`DATA.md` §2.
