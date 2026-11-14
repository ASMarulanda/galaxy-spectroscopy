# galaxy-spectroscopy

**Redshift estimation and emission-line analysis of optical galaxy spectra.**

`gspec` measures, from a flux-calibrated optical galaxy spectrum of unknown
redshift: the redshift, the continuum-subtracted integrated fluxes of
Hα λ6565, Hβ λ4863, [O III] λ5008 and [N II] λ6585, a 1σ uncertainty on each, and
a quality flag — by a joint fit of the Hα and Hβ complexes at one common
redshift.

> **Status: scientific core validated and illustrated (stage 4 of 7).** The
> package is ported, tested, reproduces its baseline bit for bit, and ships five
> reproducible figures with the diagnostics they support. Notebooks and the
> machine-learning experiment are not built. Not yet published.

## Results on the 10 synthetic spectra

Each synthetic spectrum was generated from known analytic parameters, so the
truth is arithmetic, not a fit.

| line | median error | worst | acceptance floor | inside the floor alone |
|---|---|---|---|---|
| Hα | **−0.70 %** | +41.6 % | 7 % | 8 / 10 |
| Hβ | **−14.4 %** | −94.7 % | 20 % | 7 / 10 |
| [O III] 5007 | **+0.31 %** | −65.5 % | 6 % | 7 / 10 |
| [N II] 6583 | **−1.54 %** | −7.7 % | 10 % | 10 / 10 |

Redshift: median |Δv| **1.7 km/s**, max **18.5 km/s**.
Differential recovery on 4 real spectra with a line of exactly known flux
injected: 12 / 12 comparisons pass, Hα to 1–3 %.
Uncertainties: pull RMS **1.23** over the 40 synthetic line fluxes, where both
the value and its error are known independently (the separate implementation
gives 1.31 on the same rows). The twelve injection *increments* are excluded
from that number: their errors are a conservative quadrature sum of two
correlated fits, so pooling them would flatter the calibration.

**Independently validated.** A separately written program
(`validation/independent/`, RBF continuum, bounded optimiser, bootstrap, no shared
estimator) reproduces the redshift on **20/20** spectra to a median of 6.5 km/s,
and Hα and [N II] to **0.8 % and 0.2 %** median difference. It does **not**
reproduce Hβ: the two methods differ by a median of 17.5 %, and the truth lies
between them on 10 of 10 synthetic cases. See `validation/report.md` §5.

The worst-case column is not a failure — on those rows the reported σ is itself
large and carries the tolerance — but "passes" is weaker than "accurate", which
is why `validation/report.md` pins the medians instead.

## Figures

Five figures, each with the script that draws it, a vector PDF, and the measured
numbers underneath. Regenerate all of them with
`.venv/bin/python validation/figures/make_all.py` (~6 min); see
**`docs/figures.md`** for the astrophysics, the equations and the caveats.

| figure | what it shows | what it establishes |
|---|---|---|
| **1** S01 fit | spectrum, continuum, each Gaussian, total model, residuals | this spectrum is well described: χ²/ν = 1.051, all four lines within 1 % of truth |
| **2** Hβ configurations | signed Hβ error for four fitting configurations | the *sign* of the Hβ error depends on the configuration: −14.4 %, +19.6 %, +13.7 %, +18.1 % |
| **3** pull calibration | standardised residuals, both implementations | calibrated in aggregate (RMS 1.232 / 1.306 over n = 40), no pull beyond 3σ |
| **4** injection response | measured vs known line-flux increment | the pipeline transfers added flux: Hα 98 %, [O III] 98 %, Hβ 92 % |
| **5** diagnostics | BPT + Balmer decrement for the real spectra | 2 of 6 objects are measurable; R01's decrement exceeds case B by 3.3σ |

The single most important number in that table is Figure 2's: switching
**model selection** off moves Hβ by 34 percentage points (−14.4 % → +19.6 %),
which is larger than the tolerance the pipeline is graded against. A pipeline that
reports "Hβ to 14 %" without saying which configuration produced it is reporting a
choice, not a measurement.

## What this does not prove

* **Hβ is the weak link.** Median −14.4 %, and the independent implementation
  written for this project gets +18.1 % on the same data — the two bracket the
  truth on every case and agree with each other only to 17.5 % median. The 20 %
  band used in the original validation was an *admission threshold*, not an
  accuracy claim, and this repository does not claim Hβ to 20 %.
* **The stellar-absorption model is not shown to help.** The BIC selection
  rejects it on all ten synthetic spectra, and forcing it on makes Hβ *worse*
  (−11.1 % → +19.6 %). See `validation/experiments/FINDINGS.md`.
* **σ is not reproducible across environments** where an amplitude sits on its
  bound: one row moves 19.70 → 51.98 between a container and this host.
* **n = 20.** Nothing here is a statement about galaxies in general.

Full list, with evidence: **`docs/limitations.md`**.

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

Optional batch interface:

```bash
gspec measure --cases data/cases --out out/line_measurements.csv
```

## The fit, in one block

```
f(λ) = C(λ) · Π_i [1 − d_i·G_i(λ)]  +  Σ_j A_j·G_j(λ)

C(λ) = Σ_k c_k x^k                        one polynomial per line complex
       × multiplicative Balmer absorption, width tied to the gas width
                                          ([N II] 6583/6548 and [O III] 5007/4959
                                          ratios fixed at 2.94 / 2.98)
F = A · sqrt(σ_int² + σ_LSF²) · sqrt(2π)   integral of the observed profile
σ_F² = ∇Fᵀ · Cov(p) · ∇F                  with Cov = pinv(JᵀJ)·min(χ²/dof, 1.5)
```

Redshift: a matched filter over 11 optical lines proposes candidates; a
non-negative-least-squares χ² scan polishes each one. The two complexes are then
fitted jointly with a shared `z`, which is what breaks the redshift–width
degeneracy a single window suffers from.

Derivation, bounds and failure modes: **`docs/methods.md`**.

## Validate

```bash
pytest                     # 105 tests, ~170 s
```

The tests are written against analytic truth and exact increments, not against
the pipeline's own output — see `validation/report.md` for what each one detects,
and `docs/validation-history.md` for the adversarial audit that shaped them.

## Data

20 spectra, 2.0 MB: 10 synthetic with analytic truth, 6 real SDSS DR18 spectra
with the catalogue answers removed, 4 injection copies of real spectra with an
exactly known increment.

Provenance, redaction record and the required SDSS acknowledgment: **`DATA.md`**.
Code is MIT; text and figures are CC BY 4.0; the spectra are third-party material
and are neither.

## Layout

```
src/gspec/       the package: constants → io → redshift → fit → fluxes → flags
data/            20 case FITS + truth manifests + provenance
tools/           generator, injection builder, SDSS downloader
tests/           105 tests, one per scientific failure mode
validation/      baseline, independent implementation, report, experiments
validation/figures/  five figure scripts + PNG/PDF + FIGURE_MANIFEST.json
docs/            methods · figures · limitations · validation history · decision log
```

## Citing

`CITATION.cff`. If you use this work, acknowledge SDSS as described in `DATA.md`.
