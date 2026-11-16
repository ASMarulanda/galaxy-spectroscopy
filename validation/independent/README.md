# Independent validation evidence

A second implementation of the same measurement, written separately, used to test
whether the Hβ behaviour reported in `docs/limitations.md` L2 is a property of
*this model* or of *the problem*.

## What is here

| file | what it is |
|---|---|
| `measure_independent2.py` | the independent analysis: RBF continuum, bounded optimiser, bootstrap uncertainties |
| `measure_independent.py` | its v1 module — FITS I/O and cross-correlation only; the v2 program imports it, which is why it ships |
| `run_independent.py` | a thin driver over `measure_independent2.py` |

Its outputs and the tables derived from them live in `../results/`:
`out_independent.csv` (80 rows) and `tables/*.csv`.

## Independence, stated precisely

The claim is about **dependency structure**, and it is checkable:

* `measure_independent2.py` imports `glob`, `math`, `os`, `sys`, `numpy`,
  `scipy.stats`, `scipy.optimize` and its own v1 module. It imports **nothing**
  from `gspec` and nothing from `tests/`.
* No function it defines is defined in `gspec`, except `_sig_pix` — the two-line
  line-spread quadrature `sqrt((λ σ_v/c)² + (wdisp Δλ)²)`. Both implementations
  write that physical relation; it is not copied machinery. Every multi-line
  routine is distinct.
* The fitting strategy differs: `gspec` scans trial redshifts with NNLS amplitudes
  and BIC selection over seven nested models; the independent program fits each
  emission complex in its own window with a radial-basis-function continuum and a
  bounded optimiser.

**Shared assumptions, honestly:** both adopt the same *physical model* — Gaussian
emission, multiplicative Balmer/G-band absorption, the same rest wavelengths and
atomic ratios. They share the input data. Neither shares an estimator, a continuum
basis, a fitter or a covariance computation with the other.

## Reproducing it

Dependencies are those of the package (`numpy==1.26.4`, `scipy==1.11.4`,
`astropy==7.0.0`):

```bash
cd validation/independent
python3 run_independent.py            # ../../data/cases -> ../results/out_independent_reproduced.csv
```

The original run took about 45 minutes for the 20 cases on 2 CPUs; roughly 25
minutes on a modern desktop. Run against the shipped cases, it reproduced
`out_independent.csv` byte-for-byte on the machine where it was first checked.
That agreement is environment-dependent rather than guaranteed: an earlier run on
different hardware differed in the third significant figure (median Hβ error
18.6 % then, 18.1 % here), and both land inside every acceptance band. Any figure
quoted from this output carries that spread.

> **Correction (independent re-run, 2026-11-13).** The "byte-for-byte" claim above
> did not hold on a fresh run. Re-running this program on the 20 shipped spectra
> reproduced the **redshifts exactly (20/20)** but not the fluxes: median
> difference 0.03 % (Hα), 0.68 % (Hβ), 0.15 % ([O III]), 0.03 % ([N II]);
> **maximum 14.2 % on one row** (S04 Hβ, 245.96 → 211.07), with 8 of 20 Hβ rows
> moving by more than 1 %. The program has a stochastic component (it resamples),
> so it is reproducible in its redshifts and at the sub-percent level in its
> fluxes, but not bit for bit.
>
> This does not weaken the comparison: the gap against `gspec` is stable across
> both runs (Hβ median −28.4 % vs −29.3 %, and the analytic truth falls between
> the two methods in 10 of 10 cases either way). It does mean that any single row
> of this output should be quoted with a per-row tolerance of ~15 %, not as an
> exact value. See `../report.md` §5.6.

## What the shipped outputs support

`../results/out_independent.csv` is the independent measurement of all twenty
cases, in the same column format as the reference baseline.

`../results/tables/per_line_calibration.csv` — per line and per implementation:
how many synthetic rows have the relative floor (rather than 3σ) setting the
tolerance, the worst error on those rows, the ratio floor/error, the median
absolute error, and how many rows exceed the floor raw. The honest reading of the
Hβ 20 % floor:

| | reference | independent |
|---|---|---|
| rows where the floor binds | 4 / 10 | 1 / 10 |
| worst error on those rows | 14.0 % | 1.9 % |
| median absolute error | 14.4 % | 18.1 % |
| rows above 20 % raw | 3 / 10 | 4 / 10 |

So 20 % is an **admission threshold**, not a per-row guarantee and not an
independently derived accuracy limit. For the independent implementation it barely
binds at all; most of its rows are decided by 3σ. `docs/limitations.md` L2 and
`../report.md` §5 say the same thing.

`../results/tables/hbeta_rows.csv` — every synthetic case, both implementations:
true flux, measured flux, relative error, 3σ/true, and whether the floor binds.

`../results/tables/hbeta_disagreement.csv` — per synthetic case, both Hβ fluxes
and `|F_ind − F_ref| / F_ref`. The median of those ten values is **40.3 %**; the
mean is dominated by S09 and is quoted only as such.

`../results/tables/redshift_calibration.csv` — per case: true redshift, both
redshifts, both σ_z, and both pulls (|dz|/σ_z).

`../results/tables/floor_margins.csv` and `floor_margins.md` — the same
floor-binding measurement as a population table: reference only, independent only
and pooled, for every line, with the number of binding rows and the case carrying
the worst error.

`../results/tables/calibration_by_gate.csv` — the eleven acceptance checks of the
validation, evaluated for both implementations.

## What is not shipped

The script that assembled these tables from the two CSVs is not included; the
tables ship as artifacts. Each is a direct aggregation of `out_independent.csv`,
`golden_reference.csv` and `data/truth/truth_synthetic.json`, with the definitions
given in `../report.md` §5, so any row can be recomputed by hand or in a few lines
of code.

## Scope

These files are validation evidence, not measurement code. Nothing in `src/gspec/`
imports them, they are not part of the installed package, and they are not covered
by the package's test suite.
