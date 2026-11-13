# galaxy-spectroscopy

**Redshift estimation and emission-line analysis of optical galaxy spectra.**

`gspec` measures, for a flux-calibrated optical galaxy spectrum of unknown
redshift: the redshift, the continuum-subtracted integrated fluxes of
Hα λ6565, Hβ λ4863, [O III] λ5008 and [N II] λ6585, a 1σ uncertainty on each,
and a quality flag — by a joint fit of the Hα and Hβ complexes at one common
redshift.

> **Status: under construction (Stage 2 of 7).** The scientific core is
> migrated and reproduces the validated baseline; documentation, figures and
> the validation report are being written. This repository is not yet
> published.

## What it does, in one fit

```
f(λ) = C(λ) · Π_i [1 − d_i·G_i(λ)]  +  Σ_j A_j·G_j(λ)

C(λ) = Σ_k c_k x^k                       continuum, one polynomial per line complex
       × multiplicative Balmer absorption, width tied to the gas width
A_j·G_j                                  Gaussian emission, [N II] and [O III]
                                         doublet ratios fixed at 2.94 / 2.98
F = A · sqrt(σ_int² + σ_LSF²) · sqrt(2π)  integrated flux of the observed profile
σ_F² = ∇Fᵀ · Cov(p) · ∇F                 uncertainty, from the fit covariance
```

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

The dependency pins in `pyproject.toml` are part of the science: every number
in the documentation was calibrated against them.

## Run

```bash
gspec measure --cases data/cases --out out/line_measurements.csv
```

## Data

20 spectra: 10 synthetic with analytic truth, 6 real SDSS DR18 spectra
(redacted), 4 injection copies of real spectra with an exactly known increment.
Provenance, redaction and licensing: **`DATA.md`**. Code is MIT, text and
figures are CC BY 4.0 — see `LICENSE` and `LICENSE-content`.

## Layout

```
src/gspec/     the package (constants → io → redshift → fit → fluxes → flags)
data/          20 case FITS + truth manifests + provenance
tools/         the generator, the injection builder, the SDSS downloader
tests/         the scientific tests
validation/    baseline results, the independent implementation, the report
docs/          methods, limitations, validation history, decision log
```
