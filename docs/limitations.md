# Limitations

Ranked by how much they matter. Each one states the evidence, not just the
caveat. Nothing in this list is a plan — it is what the numbers currently do not
support.

---

**L1 — The sample is 20 spectra, and nothing about galaxies follows from it.**
Ten synthetic, six real, four injected. Every statement in this repository is
about a *measurement method*, not about a population. No luminosity function, no
metallicity relation, no SFR distribution can be inferred.

**L2 — Hβ is the weak link, and the bias is a model-misspecification signature.**
Median error on synthetic Hβ: **−14.4 %** with the shipped pipeline, **−11.1 %**
with absorption removed, **+19.6 %** with absorption forced on
(`validation/experiments/FINDINGS.md`). The independent implementation gets
**+18.1 %**. Biases in opposite directions on the same data: the two methods
differ from each other by a median of **17.5 %** on Hβ across all 20 spectra, and
the analytic truth falls **between** them on **10 of 10** synthetic cases
(`validation/report.md` §5.2). That is a systematic, not scatter — neither
configuration is right, and averaging the two would produce a number with no
error model. The 20 % acceptance band used in the validation history was an
*admission threshold*, not an accuracy claim, and this repository does not claim
Hβ to 20 %.

For contrast, the same two methods agree on Hα to **0.8 %** median and on
[N II] to **0.2 %** — which is what makes the Hβ gap interpretable rather than
just noisy.

The spread across **every** configuration shipped or documented in this
repository:

| configuration | Hβ median error | Hβ worst |
|---|---|---|
| `m1a` (default, BIC-selected, no absorption on synthetic) | **−14.4 %** | −94.7 % |
| `m1b` (alternative: quadratic continuum, wider windows, absorption kept) | **+13.7 %** | +491.7 % |
| absorption forced on (`E-Abs`) | **+19.6 %** | +462.1 % |
| independent implementation | **+18.1 %** | +116.5 % |

A **33 % range straddling the truth.** All four pass the acceptance band on every
row (the band widens by 3σ), so the band cannot separate them. The defensible
statement is not "Hβ is accurate to X %" but "Hβ is uncertain at the ±15–20 %
level in a regime-dependent direction, and no configuration in this repository is
demonstrably correct".

**L3 — The stellar-absorption term is not demonstrably useful, and the claim that
it is has been withdrawn.** On all ten synthetic spectra the BIC sweep rejects it;
forcing it on degrades Hβ; the fitted depth is 0.01–0.77 of the truth and its
width ratio lands on the 1.15 prior bound in four cases. It is kept because it is
selected on five of six real spectra and because removing it entirely loses the
I04 injection — not because it is shown to improve accuracy.

**L4 — σ_F is not reproducible across environments where an amplitude sits on a
bound.** Measured: S08 Hβ σ = 19.70 on this host, 51.98 in a Debian-slim
container, with the flag flipping `ok` → `nondetection`; S04 Hα and S04 Hβ move
by 59–71 % between the same two environments. The mechanism is structural: on
S01, 81 % of σ_F(Hβ)² comes from the amplitude of the *broad* component — a
parameter the model selection has effectively zeroed — and the corresponding
eigen-direction is 8 orders of magnitude below λ_max.

The obvious fix was tested and **rejected**: conditioning the covariance on the
bound (dropping pinned directions) shrinks σ by up to 5.8×, flips six flags from
`nondetection` to `ok`, and degrades calibration (pull RMS 1.23 → 1.74, worst
pull 2.50 → 6.05). The fragility is the price of a correctly calibrated σ. Any
future replacement must be both environment-stable **and** calibrated; none of the
tested candidates is.

**L5 — Atomic ratios are applied exactly, and they are not exact.**
`I(6583)/I(6548) = 2.94` and `I(5007)/I(4959) = 2.98`. Real gas with dust,
shocks or a different density varies by tens of percent. Tying them helps
identifiability and is standard practice; treating them as exact is a modelling
choice, not a measurement.

**L6 — The detection flag is a convention on a knife-edge.**
`|F|/σ ≥ 3` is not a detection theory. 9 of 80 rows sit within a factor 1.5 of
that boundary; on S08 Hβ a 1.36× change in σ flips it. Any comparison of flag
counts between environments is therefore fragile, and the affected rows are
identifiable from the diagnostics.

Two further facts about the flag contract, both found by the Stage-3 test review
and neither previously documented:

* **`unresolved_blend` is unreachable.** `nondetection` is tested first, so
  reaching the blend branch requires `|F|/σ ≥ 3`, which implies `σ ≤ |F|/3`, which
  contradicts its own `σ > 0.9|F|`. The value survives in the output contract for
  compatibility with the original implementation. The ordering is **not** changed
  here, because reordering it would change published outputs; the contradiction is
  pinned by `tests/test_core_identities.py`.
* **The independent implementation flags 9 rows `ok` that `gspec` calls
  `nondetection`** — and never the reverse (8 of the 9 are Hβ). The flag is a
  statement about σ, and σ is the quantity the two methods disagree on most
  (`validation/report.md` §5.4).

**L7 — The continuum initialisation preserves a quirk of the original code.**
Only `c0` and `c1` are initialised (to the median flux of each window) and
bounded; `c1` is the *slope* of the first complex, not the constant of the second,
so the second complex's constant term starts at 0.0 and every other coefficient is
unbounded. This is preserved verbatim so the port reproduces the baseline
(`fit.build_parameters`, marked in place). It has no measured consequence — the
fits converge to χ²/dof ≈ 1 — but it is not what the code looks like it does.

**L8 — Absolute flux calibration is below the tightest floor and never modelled.**
SDSS spectrophotometry carries a few-percent absolute scale uncertainty, which is
smaller than the 7 % Hα floor, so it never binds. Synthetic truth is
calibration-free and injections are differential, so the validation is not
affected — but a real-world flux would inherit that error and it is not in σ.

**L9 — No template-based stellar continuum.** The polynomial is a local
approximation; a pPXF/starlight fit would be the next scientific step and is the
most likely route to a defensible Hβ (L2/L3). It is out of scope here.

**L10 — Real-data fluxes rest on four injection pairs and a factor-2 band.**
There is no exact flux truth for the six real spectra. The archival comparison is
a consistency check: the SDSS pipeline and MPA-JHU disagree by a median of 139 %
per line on this sample, so the band is wide by necessity.

**L11 — The broad-line AGN (R05) is not validated at all.**
No injection copy exists, because no valid method recovered a +144 Hα increment
reproducibly (+55 %, +99 %, +100 % under three configurations; χ²/dof 2.89). Its
Hβ σ is 690 under the shipped estimator and 149 conditioned (L4), which is a
measure of how little the model constrains it.

**L12 — Cross-version behaviour of scipy is untested.** The solver settings are
pinned (`ftol=xtol=1e-10`, `max_nfev=3000`, `x_scale="jac"`) and the environment
is pinned in `pyproject.toml`, but no experiment has compared `least_squares`
across scipy versions. Treat any change of scipy major or minor version as a
change of numerics.

**L13 — Storage precision caps every reproducibility claim.** The case files are
float32 (`format='E'`), i.e. ~1e-7 relative quantisation. A claim of agreement
tighter than ~1e-6 relative cannot be tested against this data
(`DATA.md` §5, `tests/test_io.py`).

**L14 — The acceptance band cannot tell a good model from a bad one.**
Every row of all four Hβ configurations in `fig2` passes the acceptance test,
including a **491.7 %** error on S06 and a **−94.7 %** error on S09; the large
errors pass because the reported σ is itself large and widens the tolerance to 3σ.
The band is a rejection test for broken code, not a certificate of accuracy. This
is why the medians are pinned by tests and why the README reports the spread
across configurations rather than a single accuracy figure.

**L15 — Figures are reproducible only on the pinned stack.**
The scripts strip timestamps, so a re-run on the same versions is byte-identical
(asserted in `tests/test_plotting.py`). A different matplotlib, freetype or font
backend changes antialiasing and text metrics, and the bytes differ even though the
science does not. `FIGURE_MANIFEST.json` records the stack used, so a mismatch is
diagnosable rather than mysterious.

**L16 — The diagnostics are two objects.**
2 of 6 real spectra can be placed on a BPT diagram and 2 have a Balmer decrement,
because [O III]/Hβ needs a measured Hβ and four objects have none. No population
statement, no reddening law and no classification follows from n = 2. R01's
A_V = 1.92 rests on an Hβ that this document (L2, L3) calls the least reliable
number the pipeline produces: a 20 % low Hβ moves it to 2.62. The *measured*
statement is the 3.3σ excess over case B, not the extinction.

**L17 — The numbers reproduce on this stack, not across BLAS builds.**
V0.3.0's CI ran the suite on GitHub's `ubuntu-latest` (numpy 1.26.4 wheel, bundled
OpenBLAS) and 6 of 195 tests failed. Reproduced exactly on a wheel-based virtualenv
here — the same six tests, the same numbers — so the cause is the linear-algebra
build, not the runner.

How far apart, measured across all 80 rows of the baseline (wheel stack against the
reference stack):

| field | worst relative | where | rows over 1e-4 | over 1e-3 |
|---|---|---|---|---|
| redshift | 2.0e-8 | I01 Hα | 0 | 0 |
| flux | **6.4e-4** | I01 Hβ | 2 | 0 |
| `redshift_err` | 5.0e-3 | S04 Hα | 8 | 4 |
| `flux_err` | **7.1e-1** | S08 Hα | 13 | 11 |

And two discrete outputs move:

* **one flag flips** — S08 Hβ is `ok` on the reference stack and `nondetection` on the
  wheel stack, because its σ moves by 2.6× and the classification follows σ;
* **the model-selection sweep lands on a different configuration for 2 of 20 cases**
  (R01 and I01: `with-abs, ndeg=3` versus `with-abs, ndeg=3, no broad Ha`). The χ²
  difference between those two candidates is ~1e-4 — a margin smaller than the
  environmental noise in χ², so the parsimony decision is a coin toss. The flux
  consequence is the 6.4e-4 above: I01 Hβ 185.4700 → 185.3509.

Three behaviours, then, and they need different tests.

*Redshifts are stable* to 2e-8, and *χ² to 1e-11*. A different summation order inside
the least-squares kernel moves the solution in its last digits, and on 78 of 80 rows
that is the whole of it.

*Fluxes are stable on 78 of 80 rows* (to ~1e-6) and move by **6.4e-4** on the two
where the model-selection sweep chose differently (above). That is 0.064 % — 100×
below the smallest accuracy floor in this project (6 % for [O III]) — but it is not
rounding, and the portable tolerance has to accommodate it.

*σ on bound-pinned rows is not stable at all*, by up to a factor of 3.5, and it drags
a flag with it. The mechanism, measured on S08 rather than asserted:

1. The model carries a broad component in each line complex whose amplitude the fit
   drives to zero — `amp_Ha_b = 1.6e-12`, `amp_Hb_b = 3.5e-21`, both at their bound,
   both reported in `at_bounds`.
2. Its **width** (`broad`) is therefore unconstrained by the data. At `broad = 5.53`
   and `broad = 7.50` — a 26 % difference — χ² is identical to **1.6e-11**. The
   parameter is unidentifiable, and the optimiser stops wherever it enters that flat
   valley.
3. The width still sets the *coupling* between the broad and narrow amplitudes.
   Along the valley the second-smallest eigenvalue of JᵀJ moves from **7.7e-3 to
   1.3e-1** (17×), in the direction of the broad/narrow Hβ split.
4. σ = sqrt(gᵀCg) inherits that change, and moves by 3.5× — while every reported
   flux, the redshift and χ² are identical.
5. The flag follows σ. `classify()` compares the flux with its own uncertainty, so a
   2.6× larger σ turns S08's Hβ from `ok` into `nondetection`. A discrete scientific
   output moves because a continuous one did.

The pseudo-inverse truncation is *not* the culprit: both stacks drop the same
direction (the one at eigenvalue ~1e-15 dominated by `broad`). What differs is where
in the valley the optimiser lands, and that is decided by the last digits of a BLAS
kernel.

**What the golden baseline guarantees, and what it does not.** It guarantees that on
the reference stack the pipeline produces these exact numbers, so any change in the
code that moves them is visible (`characterisation` tests). It does **not** guarantee
that another environment produces them. Until V0.3.1 the suite asserted the stronger
claim — including a CI step that required the baseline hash to match — which was a
defect in the tests, not a property of the pipeline.

**Which tests are portable, and which are not.**

| suite | runs | asserts |
|---|---|---|
| `pytest -m "not characterisation"` | either stack measured, and in CI | redshift to 1e-6; flux to 1e-3; σ, `redshift_err` and `flux_err` to their baseline values on the 44 rows whose line complex has no active bound (of 80); flags on those same rows; the covariance *rule* (analytic), finiteness, positive semi-definiteness, bound reporting; that each flag is consistent with the flux and σ the same result reports, everywhere; accuracy against analytic truth within the acceptance floors; failure handling |
| `pytest -m characterisation` | the reference stack only | the golden baseline to its printed digits, and σ on every row |

Tolerances and their justification are in `tests/test_regression.py`. **None of them is
a portability guarantee.** The flux bound (1e-3) is 1.6× the largest cross-environment
flux drift measured (6.4e-4); the redshift bound (1e-6) is 50× above its drift (2e-8),
which leaves real room. The error columns are compared only on rows whose complex has
no active bound, because on bound rows their drift reaches 5e-3 (`redshift_err`) and
7.1e-1 (`flux_err`). The 202 tests run on both stacks that
have been examined because they assert what was measured to hold across those two — not
because the numbers are shown to reproduce everywhere, which would need stacks nobody
has run.

**What is deliberately not fixed.** A component pinned at zero amplitude keeps a
free width, which is where the flat direction comes from. Removing it would change
the model and therefore the numbers, so it is recorded as open question O6 rather
than changed in a maintenance release. Nothing here claims that σ is reproducible
across BLAS builds, or that it is calibrated — calibration is a separate question,
tested by pull distributions within a stack (`tests/test_uncertainty.py`).
