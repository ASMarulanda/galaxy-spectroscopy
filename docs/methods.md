# Methods

Every equation in this document is implemented in one named function; the
reference is given beside it. Nothing here is aspirational — if a number is
quoted, it was measured on the 20 spectra in `data/`.

```
src/gspec/constants.py     the physical vocabulary
src/gspec/io.py            the FITS contract and the mask policy
src/gspec/lsf.py           the width relation
src/gspec/lines.py         which lines exist, in which window, tied how
src/gspec/continuum.py     the continuum and the absorption product
src/gspec/redshift.py      candidate generation and the linear scan
src/gspec/fit.py           the bounded solve
src/gspec/covariance.py    the uncertainty
src/gspec/fluxes.py        the integrated flux
src/gspec/select.py        the parsimony sweep
src/gspec/diagnostics.py   line ratios, CCM89 reddening, BPT curves
src/gspec/plotting.py      the figures, re-derived from the fitted parameters
```

## 1. The measurement problem

A single-fibre optical spectrum is a flux **density** `f_λ(λ)` in
`1e-17 erg s⁻¹ cm⁻² Å⁻¹` on a log-λ grid. It is the sum of four things that
overlap in wavelength:

```
f_λ(λ) = stellar continuum          (smooth, 10-100 Å scales)
       × stellar absorption          (multiplicative dips: Balmer, Ca II, G band)
       + nebular emission            (Gaussian lines, 50-500 km/s)
       + noise                       (photon, read, sky subtraction)
```

The task is to recover, for each object, one redshift and four integrated line
fluxes with uncertainties. The four lines are chosen because of what they
diagnose *together*:

| line | rest (vacuum, Å) | diagnoses |
|---|---|---|
| Hα | 6564.61 | star formation, ionised gas mass |
| Hβ | 4862.68 | Balmer decrement with Hα → extinction |
| [O III] 5007 | 5008.24 | ionisation parameter, AGN/shock |
| [N II] 6583 | 6585.27 | metallicity proxy, BPT classification |

Hα/Hβ gives the Balmer decrement; [N II]/Hα and [O III]/Hβ give the BPT axes.
None of them works without deblending and continuum subtraction, which is why the
project measures all four and not just the brightest.

## 2. Redshift

```
λ_obs = λ_rest (1 + z)
```

Two-stage search, because each stage alone fails:

**2a. Candidate generation** — `redshift.matched_filter_scan`

A weighted matched filter over 11 optical lines. For each trial redshift and each
template line, the local continuum is estimated from sidebands 2–5 pixels away:

```
S(z) = Σ_l w_l · Σ_i (f_i − median(sideband)) · K_i · ivar_i / sqrt(Σ_i K_i² ivar_i)
K_i  = exp(−½((λ_i − λ_c)/σ_a)²)          σ_a = 2.2 px
```

with the relative weights `ZFILTER` (Hα 1.00, [N II] 0.35/0.12, Hβ 0.35,
[O III] 0.32/0.11, Hγ 0.16, [S II] 0.13/0.11, [O II] 0.45/0.45).
1400 trial redshifts over 0.008–0.102.

*Why a filter and not a χ² scan:* a global χ² is sensitive to the continuum
shape and, at low S/N, aliases onto the wrong member of the Hα/[N II] doublet —
a 943 km/s trap. The filter is insensitive to the slowly varying continuum, but
locates a peak coarsely.

**2b. Polish** — `redshift.redshift_peaks` → `linear_redshift_scan`

Within ±0.004 of each of the top candidates (kept ≥0.0020 apart in z), the
redshift is refined by a 41-point χ² scan in which **all amplitudes are solved
exactly** by non-negative least squares at each trial z. Two fixed widths
(110 and 250 km/s) bracket the intrinsic width.

*Why NNLS:* at fixed `z` the model is linear in every amplitude, so the
amplitudes have a closed-form non-negative solution and only `z` and the widths
need searching. This is what makes a 41-point scan cheap.

**2c. Final refinement** — inside `fit.fit_case`

`z` becomes a free parameter bounded to `z₀ ± 0.0035` (about 1000 km/s). That
window is the pipeline's main external assumption: the seed must already be close.

Measured on the ten synthetic cases: median |Δv| **1.7 km/s**, max **18.5 km/s**.

## 3. The model

`continuum.continuum_at`, `model.model_segments`, `model.residuals`.

```
f(λ) = C(λ) · Π_i [1 − d_i G_i(λ)]  +  Σ_j A_j G_j(λ)

C(λ) = Σ_{k=0..n} c_k x^k ,   x = (λ/(1+z) − λ̄_c) / Δλ_c
```

One polynomial per complex (independent coefficients), `n = 3` for the reference.

**Two windows, one redshift.** The Hα complex (6420–6790 Å rest) and the Hβ
complex (4250–5090 Å rest) are fitted simultaneously with a shared `z`. A single
complex cannot separate redshift from line width; two complexes 1700 Å apart can,
because `z` stretches both by the same fraction while the widths are independent.
The continuum is allowed to differ between them, which absorbs the real shape
difference without forcing a global model.

| complex | components |
|---|---|
| Hα | [N II] 6548 (tied), Hα narrow, Hα broad, [N II] 6583, [S II] 6716, [S II] 6731 |
| Hβ | Hγ, Hβ narrow, Hβ broad, [O III] 4959 (tied), [O III] 5007 |

Sibling lines ([S II], [O I], Hγ, higher Balmer, [O II]) are present in the data
and must be modelled, or they contaminate the graded lines.

**Ties.** `I(6583)/I(6548) = 2.94`, `I(5007)/I(4959) = 2.98`
(`constants.NII_RATIO_*`). These are set by transition probabilities, not by the
galaxy, so they are applied exactly. Treating them as exact is a modelling
choice — real gas with dust and shocks varies by tens of percent (L5).

**Bounds** (`fit.build_parameters`), 25 parameters for the full reference model:

| parameter | count | bounds |
|---|---|---|
| `z` | 1 | z₀ ± 0.0035 |
| `narrowA`, `narrowB`, `narrowO` (log σ_v) | 3 | ln[25, 600] km/s |
| `broad` (log σ_v) | 1 | ln[250, 1800] km/s |
| `abs_scale`, `abs_ratio`, `gband` | 3 | [0,1.5], ln[1.15,6], [0,1] |
| amplitudes ≥ 0 | 9 | [0, ∞) |
| continuum coefficients | 8 | −50…500 (see L7) |

## 4. Stellar absorption, and the identifiability prior

The absorption width is **not free**. It is tied to the narrow emission width:

```
σ_abs = σ_narrow,B · exp(abs_ratio),      abs_ratio ∈ [1.15, 6.0]
```

That is a physical prior, not a fit trick: stars fill the whole galaxy while H II
regions are kinematically cold, so stellar absorption must be broader than the
nebular lines. Without the tie, "deep narrow absorption" and "shallow broad
absorption" fit the same data equally well and the fitted Hβ amplitude absorbs
the difference — measured at +10…+80 % bias.

### 4.2 What the term actually does — a negative result

The fitter runs a BIC sweep over nested models (`select.fit_with_model_selection`).
On this dataset the outcome is:

| family | cases with absorption selected |
|---|---|
| synthetic (S01–S10) | **0 of 10** |
| real (R01–R06) | **5 of 6** (all but R03) |

And `validation/experiments/e_abs_absorption.py` shows that forcing the term on
makes Hβ *worse* (median +19.6 % instead of −11.1 %), that the basis is not the
reason (matching the generator's own weights changes little), and that the fitted
depth is 0.01–0.77 of the true depth with the width ratio pinned at its 1.15
bound in four cases.

**Consequence for what this project may claim:**

* permitted — "the fitter contains a multiplicative Balmer absorption model, and
  the BIC selection keeps it on 5 of 6 real spectra";
* **not** permitted — "the absorption model improves Hβ", "absorption modelling
  is required for Hβ", "the absorption tie is why Hβ is accurate".

Neither configuration is correct: without the term Hβ is biased low by ~11–14 %,
with it biased high by ~17–20 %. See `docs/limitations.md` L2 and
`validation/experiments/FINDINGS.md`.

## 5. Integrated flux

`fluxes.integrated_flux`

```
σ_int(Å)  = λ_obs · σ_v / c
σ_LSF(Å)  = wdisp[pixels] · Δλ(pixel)          Δλ local, never global
σ_obs(Å)  = sqrt(σ_int² + σ_LSF²)
F_line    = Σ_components A · σ_obs · sqrt(2π)
```

Three things must be right, and each has a measured failure mode:

1. **Quadrature.** Both the intrinsic profile and the instrumental response are
   Gaussian, so their convolution adds *variances*. Adding widths linearly
   inflates F by **+39.6 %** on S01 Hα; ignoring the LSF deflates it by 18.9 %.
2. **The integral is of the observed profile.** Convolution preserves area, so
   no deconvolution before integrating.
3. **Sum over components.** Hα means narrow + broad (a definition, not a
   measurement). [N II] 6583 means the 6583 component only.

**Units.** The result is the FITS density unit times Å = `1e-17 erg s⁻¹ cm⁻²`.
No `(1+z)` correction is applied: this is the observed flux.

**Density vs integral.** `F/A = sqrt(2π)·σ_obs` — a factor of **5.98** on S01 Hα,
not the "two orders of magnitude" that an older note claimed. Reporting the peak
density instead of the integral is wrong by that factor, not by 100.

## 6. Uncertainty

`covariance.covariance`, `covariance.sigma_of_flux`

```
Cov(p)    = pinv(JᵀJ) · min(χ²/dof, 1.5)
σ_F²      = ∇Fᵀ · Cov(p) · ∇F
∇F        a forward finite difference of the same flux function that is reported
```

**The cap.** For a least-squares problem whose residuals already carry `√ivar`
weights, `JᵀJ` is the inverse covariance; scaling by `χ²/dof` is the standard
correction for underestimated errors. On real spectra `χ²/dof` measures *model
mismatch*, not photon noise, so the full correction over-corrects. Measured on
this dataset, the cap binds on 2 of 8 cases examined: R01 (χ²/dof 1.5217 →
×1.01) and R05 (2.8914 → ×1.39). The widely quoted "35 % on 40σ lines" came from
an earlier configuration and is history, not current behaviour.

**What σ is not.** It is the local linearised statistical error. It excludes
absolute flux calibration (~few %), the choice of continuum basis, and the
Balmer-ratio and absorption assumptions.

**Calibration, measured:** pull RMS **1.232** over the 40 synthetic line fluxes
(`median σ/F = 0.052`, no pull beyond 5σ; the independent implementation gives
1.306 on the same rows). The historical "1.08 over 52 points" pooled these with
the twelve injection increments, whose errors are a conservative quadrature sum —
it is retained in the report as a regression value, not as a calibration
statement. `tests/test_uncertainty.py` enforces both.

## 7. Model selection

`select.bic`, `select.fit_with_model_selection`

```
BIC = χ² + k · ln(N)
```

Two redshift seeds are fitted; the lower χ² wins. Then a parsimony sweep over
seven nested models (broad-Hβ on/off, absorption on/off, `ndeg` 1–3, broad-Hα
on/off) accepts a simpler model whenever its BIC is lower. A broad component or
a full cubic can otherwise absorb a weak line's flux at low S/N.

The sweep runs **only** for `m1a`: it is a statement about the reference model
family, not about a variant that has already fixed its shape.

The selected model travels with the result (`SpectrumResult.model_selected`),
because §4.2 shows that *which* model won is itself a scientific result.

## 8. Flags

`flags.classify`

```
failed             no finite estimate, or σ ≤ 0
nondetection       |F| / σ < 3
unresolved_blend   broad component present AND σ > 0.9 |F|      <- unreachable
ok                 otherwise
```

The flag describes the measurement, using the numbers in the same row. There is
no hidden significance threshold. The 3σ boundary is a convention: 9 of the 80
rows sit within a factor 1.5 of it, so for ~11 % of rows the flag is a
knife-edge (L6).

**`unresolved_blend` cannot be emitted.** Because `nondetection` is tested first,
reaching the blend branch requires `|F|/σ ≥ 3`, i.e. `σ ≤ |F|/3`, which
contradicts its own `σ > 0.9|F|`. The value survives in the output contract for
compatibility with the original implementation, and the ordering is deliberately
not changed here — reordering it would change published outputs. The
contradiction is pinned by a test so that a future change cannot pass unnoticed
(`tests/test_core_identities.py`).

## 9. Variants

`variants.py` — configuration objects, never forked code.

| | continuum | windows | Hγ weight | selection |
|---|---|---|---|---|
| `m1a` | cubic | 6420–6790 / 4250–5090 | 0.95 | BIC sweep |
| `m1b` | quadratic | widened by 8 Å | 1.05 | none |
| `m2` | cubic | as m1a | — | no absorption term |

`m1b` is an **alternative configuration, not an independent implementation**: it
shares this fitter and differs in three constants. Genuine independence is the
separately written program in `validation/independent/`.

### 9.1 What `m1b` actually does, measured

It is not a subtle variation. On the ten synthetic spectra:

| line | `m1a` median | `m1a` worst | `m1b` median | `m1b` worst |
|---|---|---|---|---|
| Hα | −0.7 % | +41.6 % | −0.6 % | +49.4 % |
| Hβ | **−14.4 %** | −94.7 % | **+13.7 %** | **+491.7 %** |
| [O III] | +0.3 % | −65.5 % | −0.1 % | −63.4 % |
| [N II] | −1.5 % | −7.7 % | −1.9 % | −7.7 % |

Three observations that matter:

1. **Both configurations pass the acceptance band on every row** (10/10 for each,
   every line) because the band widens by 3σ. So the band does not separate them,
   which is another reminder that "passes" ≠ "accurate".
2. **`m1b` moves Hβ to the other side of the truth** (median +13.7 % against
   −14.4 %) with a similar median magnitude, and it is catastrophic on the worst
   row (S06: 295.9 against a true 50.0). On the same faint cases the reference's
   BIC sweep *prevents* that overshoot — a benefit of model selection the earlier
   documentation did not state.
3. **`m1b` is not systematically better or worse**: it wins the Hβ comparison on
   4 of 10 cases, loses on 6.

The honest summary across every configuration in this repository:

```
m1a  (shipped, pre-selection)      Hbeta median  -14.4 %
abs_on  (E-Abs, forced)            Hbeta median  +19.6 %
m1b  (alternative configuration)   Hbeta median  +13.7 %
independent implementation         Hbeta median  +18.1 %
```

Four defensible configurations spread Hβ across a **33 % range**, straddling the
truth. That spread, not any single number, is the measurement's real Hβ
uncertainty — see `docs/limitations.md` L2.

## 10. Diagnostics beyond the flux

`diagnostics.py`, documented in **`docs/figures.md` §5**.

Three quantities are derived from the reported fluxes and nothing else:

```
R  = F1/F2                    σ_R² = (σ1/F2)² + (F1·σ2/F2²)²
D  = Hα/Hβ
E(B−V) = 2.5/(k(Hβ)−k(Hα)) · log10(D/2.86)
```

with `k(λ)` the Cardelli, Clayton & Mathis (1989) curve at `R_V = 3.1`, validated
against published `R_λ` values in `tests/test_diagnostics.py`.

The distinction the code enforces is between a **measured** ratio and an
**inferred** property: `Ratio` carries `detected` and refuses to hold a value for
a line below 3σ, and `Extinction` carries `valid` plus an `av_upper` limit when the
decrement does not exceed case B. An inference is never returned without the four
assumptions that produce it (case B, foreground screen, negligible stellar Balmer
absorption, CCM89 R_V = 3.1), because those are what a reader needs in order to
discount the number.
