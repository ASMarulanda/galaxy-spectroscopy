# Figures and astrophysical diagnostics

Every figure in this project is produced by a script in `validation/figures/`, from
committed data, with no number typed in by hand. Each section below states the
astrophysical question, the equations, the code that draws it, what the figure
supports, and — the part that matters most — what it does not establish.

```bash
.venv/bin/python validation/figures/make_all.py      # all five, ~6 min
```

(matplotlib 3.6.3 is the figure stack and is installed by the `[dev]` extra; the
core package deliberately does not depend on a plotting library.)

Individual figures:

```bash
.venv/bin/python validation/figures/fig1_s01_fit.py            # ~5 s
.venv/bin/python validation/figures/fig2_hbeta_configurations.py   # ~4 min
.venv/bin/python validation/figures/fig3_pull_calibration.py   # ~2 s
.venv/bin/python validation/figures/fig4_injection_recovery.py # ~1 min
.venv/bin/python validation/figures/fig5_diagnostics.py        # ~40 s
```

`make_all.py` writes `FIGURE_MANIFEST.json`: software versions, the SHA-256 of
every input, and the SHA-256 of every output. The figures carry no timestamps, so
a re-run on the same stack reproduces them **byte for byte** — which is asserted
in `tests/test_plotting.py::test_figures_are_byte_reproducible`.

---

## 1. The fit itself — `fig1_s01_fit.{png,pdf}`

**Question.** Does the model actually describe the spectrum, or is it a smooth
function threaded through the line peaks?

**Data.** `data/cases/S01.fits`, a synthetic spectrum generated from known
analytic parameters at z = 0.040000, so every curve in the figure can be checked
against an exact answer.

**Method.** The figure re-evaluates the fitted model from `result.parameters`
alone — it does not reuse any array the optimiser kept:

```
f(λ) = C(λ) + Σ_j A_j · G(λ; μ_j, σ_obs,j)          μ_j = λ_j(1+z)
C(λ) = Σ_k c_k x^k  ×  Π_i [1 − d_i G_i(λ)]
σ_obs  = sqrt(σ_int² + σ_LSF²)                      (quadrature, not sum)
```

**Code.** `gspec.plotting.model_curves` (decomposition) and
`gspec.plotting.fig_spectrum_fit` (layout). The decomposition is verified twice in
the test suite: `total == continuum + Σ components` to 1e-12, and
`Σ ((data − model)/σ)² == result.chi2` to 1e-8 — i.e. the replotted model *is*
the fitted model.

**How to read it.** Panel (a) shows the whole spectrum with the two fitted windows
shaded: the fit sees 20 % of the pixels, and the residual panel (d) is drawn over
every pixel it used, so structure outside the zooms cannot hide. In (b) and (c)
the dashed curve is the continuum and each filled shape is one Gaussian added on
top of it; the red curve is their sum. The narrow and broad components of a line
are separate shapes with separate widths, and only the sum is reported as the
line flux.

**What it supports.** On S01 the four lines land at Hα +0.15 %, Hβ −0.99 %,
[O III] +0.28 %, [N II] +0.04 % of truth with χ²/ν = 1.051 for 20 parameters over
1009 degrees of freedom. The residuals in (d) are structureless at the ±1σ level:
this particular spectrum is well described.

**What it does not establish.** One well-fitted case is not a validated pipeline.
On S06 and S09 the same model misses Hβ by −63 % and −94.7 %, and those residuals
are *also* unremarkable — the fit is not failing loudly, it is answering a
slightly different question. The figure also cannot show that the decomposition is
unique: a narrow Gaussian plus a broad one can mimic a single slightly broader
line, which is exactly the degeneracy §2 measures.

---

## 2. Hβ depends on the configuration — `fig2_hbeta_configurations.{png,pdf}`

**Question.** Hβ is the weakest number the pipeline produces. Is that a property
of the *data*, or of the *model*?

**Data.** `data/cases/S01..S10.fits` with `data/truth/truth_synthetic.json`
(analytic truth), plus `validation/results/out_independent.csv` for the fourth
configuration.

**Method.** Four configurations, each fit end-to-end on all ten spectra:

| label | what it is | what it isolates |
|---|---|---|
| `BIC selection (default)` | cubic continuum, absorption available, 7-candidate BIC sweep | the shipped pipeline |
| `no selection` | identical constants, sweep off | the effect of **model selection alone** |
| `m1b` | quadratic continuum, windows 8 Å wider, Hγ weight 1.05, no sweep | the effect of the other constants |
| `independent` | separately written program (RBF continuum, bounded optimiser, bootstrap) | a genuinely different method |

Signed error is `100·(F_meas − F_true)/F_true`. Panel (a) plots every case so the
outliers stay visible; panel (b) shows medians with the cases overlaid, because a
mean would hide S06.

**Code.** `gspec.plotting.fig_hbeta_configurations`; the driver is
`validation/figures/fig2_hbeta_configurations.py`.

**Measured.**

| configuration | median signed error | worst case |
|---|---|---|
| BIC selection (default) | **−14.4 %** | −94.7 % (S09) |
| no selection | **+19.6 %** | +462.1 % (S06) |
| m1b | **+13.7 %** | +491.7 % (S06) |
| independent | **+18.1 %** | +116.5 % (S06) |

**What it supports.** The sign of the Hβ error flips with the configuration. The
default under-predicts by 14 %; keeping the absorption term but dropping the
sweep over-predicts by 20 %; m1b over-predicts by 14 %. Switching model selection
off moves Hβ by 34 percentage points — **more than the total spread the acceptance
band allows**. It also reproduces, from an independent route, the +19.6 % that
experiment E-Abs measured when absorption was forced on.

**What it does not establish.** That any of the four is right. All four pass the
acceptance test on every row, which says the test is loose, not that the answers
agree. The configurations are also not a factorial design — two constants change at
once between `default` and `m1b` — so this is a sensitivity map, not a model
comparison. And nothing here transfers to real spectra, where the continuum and
the absorption are genuinely present rather than synthetic.

---

## 3. Are the uncertainties honest? — `fig3_pull_calibration.{png,pdf}`

**Question.** Does the reported σ describe the error actually made?

**Method.** For a measurement `F`, reported error `σ_F`, and known truth `F_true`:

```
pull  p = (F − F_true) / σ_F
```

If σ_F is the standard deviation of the measurement error, the pulls are drawn from
N(0,1): RMS ≈ 1, 68 % within ±1. RMS ≫ 1 means the errors are underestimated; RMS
≪ 1 means the uncertainty is padded.

**Data.** `validation/results/golden_reference.csv` (this pipeline) and
`out_independent.csv` (the separate implementation), against the analytic truth of
the ten synthetic spectra: **n = 40** (4 lines × 10 spectra).

**Code.** `gspec.plotting.fig_pulls`;
`validation/figures/fig3_pull_calibration.py::load_pulls`.

**Measured.**

| estimator | n | pull RMS | median \|pull\| | within 1σ | max \|pull\| |
|---|---|---|---|---|---|
| gspec | 40 | **1.232** | 0.958 | 57 % | 2.50 |
| independent | 40 | **1.306** | 0.712 | 65 % | 2.99 |

**What it supports.** Both estimators are calibrated in the weak sense that
matters operationally: no pull exceeds 3σ, none exceeds the 5σ failure threshold,
and the two independent estimators land within 6 % of each other on this statistic.

**What it does not establish.** Calibration. Aggregate RMS near 1 is *necessary,
not sufficient*: n = 40, the four lines of a case come from one fit and are
therefore correlated, and the synthetic spectra are smoother than real data — a
pull distribution computed on idealised noise can be correct while the pipeline is
wrong on galaxies. It also says nothing about correlated errors between lines,
which is what a line *ratio* needs.

Two further caveats, both measured. First, `docs/limitations.md` L4: on rows where
an amplitude sits on its bound, σ moves by up to 2.6× between a container and this
host. Second, the number quoted in earlier versions of these documents —
1.084 over 52 points — **pooled these 40 pulls with the twelve injection
increments**, whose errors are a conservative quadrature sum and whose pulls
therefore sit at RMS 0.195. That pooled figure flatters the calibration and is no
longer used as a calibration statement; it survives only as a regression value in
`tests/test_uncertainty.py`.

---

## 4. Does the pipeline transfer added flux? — `fig4_injection_recovery.{png,pdf}`

**Question.** If a line's flux changes by a known amount, does the measurement
change by that amount?

**Data.** `data/cases/I01..I04.fits`, each its base real spectrum plus one Gaussian
line of exactly known flux (40 % of the base pipeline flux), and
`data/truth/truth_injections.json`.

**Method.** A *differential* measurement:

```
ΔF_measured = F(I_case) − F(base_case)         vs      ΔF_known
```

Both spectra share the same continuum, the same absorption strength and the same
neighbouring lines, so a common error in the continuum cancels. Error bars are
`sqrt(σ_I² + σ_base²)` from the two fits — deliberately conservative, since the two
fits are strongly correlated and the true uncertainty on the difference is smaller.

**Code.** `gspec.plotting.fig_injection`;
`validation/figures/fig4_injection_recovery.py`.

**Measured.** Median recovery fraction: Hα **98.1 %**, [O III] **98.2 %**, Hβ
**92.0 %** (n = 4 each). Individual cases: Hα 97.7–99.9 %, [O III] 88.5–100.3 %,
Hβ 80.6–103.5 %.

**What it supports.** The pipeline is close to linear in added line flux and does
not systematically swallow a real increment. Hβ's 92 % median with two cases near
81 % is the same weakness seen in §2, now visible differentially.

**What it does not establish.** The absolute flux scale. A pipeline that
underestimated every Hβ flux by 40 % would still recover a 40 % increment
correctly — this test measures the *slope*, not the intercept. Nor does it test
redshift robustness, the noise model, or behaviour at flux levels below the
injections used here (8–570 units).

---

## 5. Line-ratio diagnostics — `fig5_diagnostics.{png,pdf}`

**Question.** Do the six real spectra support the standard emission-line
diagnostics — the BPT classification and a Balmer-decrement extinction?

**Answer, measured first: mostly not.** Detection at the pipeline's own 3σ
criterion:

| case | Hα | Hβ | [O III] 5007 | [N II] 6583 | BPT position | Hα/Hβ |
|---|---|---|---|---|---|---|
| R01 | 61.7σ | 7.3σ | 4.2σ | 43.1σ | yes | **5.270 ± 0.731** |
| R02 | 99.3σ | 38.2σ | 117.9σ | 32.0σ | yes | **3.264 ± 0.092** |
| R03 | 6.9σ | 1.0σ | 2.6σ | 31.3σ | no | — |
| R04 | 4.0σ | 0.5σ | 6.0σ | 23.0σ | no | — |
| R05 | 7.0σ | 0.3σ | 20.3σ | 19.2σ | no | — |
| R06 | 2.0σ | 0.1σ | 0.6σ | 4.9σ | no | — |

**Two of six objects can be placed on a BPT diagram, because [O III]/Hβ needs a
measured Hβ.** That is the finding, and the figure states it rather than
compensating for it.

### 5.1 BPT

```
x = log10([N II] 6583 / Hα)          y = log10([O III] 5007 / Hβ)
```

Reference curves, drawn as references and not as a classification:

* **Kauffmann et al. (2003)**, the empirical star-formation ceiling:
  `y = 0.61/(x − 0.05) + 1.3`, valid for `x < 0.05`.
* **Kewley et al. (2001)**, the theoretical maximum-starburst line:
  `y = 0.61/(x − 0.47) + 1.19`, valid for `x < 0.47`.

Both are plotted only over their stated domains; the region below Kauffmann is
star-forming, between the curves is composite, above Kewley is AGN.

| case | log([N II]/Hα) | log([O III]/Hβ) | position |
|---|---|---|---|
| R01 | −0.471 ± 0.012 | −0.787 ± 0.120 | below both curves — star-forming |
| R02 | −1.268 ± 0.014 | +0.533 ± 0.012 | below Kauffmann — star-forming |

**What it supports.** Both measurable objects sit in the star-forming region, with
error bars small enough that the placement is not driven by uncertainty (R02 sits
25σ below the Kauffmann line, R01 8σ below it).

**What it does not establish.** Any classification of a *sample*, and any statement
about the four objects that could not be placed — they are not "AGN", they are
unmeasured. The BPT diagram also assumes the lines come from one ionisation
mechanism and that the ratios are extinction-corrected; R01's [N II]/Hα is measured
on an uncorrected spectrum, and the correction itself is not available (§5.2).

### 5.2 Balmer decrement and extinction

```
D = Hα/Hβ (observed)
E(B−V) = 2.5 / (k(Hβ) − k(Hα)) · log10(D / 2.86)
A_V = R_V · E(B−V),  R_V = 3.1
```

`k(λ)` is the Cardelli, Clayton & Mathis (1989) curve. It is implemented in
`gspec.diagnostics.ccm89_k` and anchored in `tests/test_diagnostics.py` to the
published `R_λ = A_λ/A_V` values — V = 1.000, B = 1.324, J = 0.282, K = 0.114 — so
a wrong polynomial branch fails a test rather than a reader. For Hα/Hβ the
coefficient `2.5/(k(Hβ) − k(Hα))` evaluates to **2.33**.

| case | Hα/Hβ | E(B−V) | A_V |
|---|---|---|---|
| R01 | 5.270 ± 0.731 | 0.618 ± 0.108 | **1.92 ± 0.33** |
| R02 | 3.264 ± 0.092 | 0.134 ± 0.013 | **0.41 ± 0.04** |

**The assumptions, which the figure's secondary axis only means if all four hold:**

1. **Case-B recombination**, T ≈ 10⁴ K, n_e ≈ 10²–10⁴ cm⁻³ (Osterbrock & Ferland),
   so the intrinsic ratio is 2.86.
2. **A foreground dust screen**, not a mixed geometry.
3. **The Balmer lines are not absorbed by stellar photospheres.** This is the weak
   link for any post-starburst spectrum: Balmer absorption depresses Hβ
   disproportionately and inflates the apparent decrement exactly as dust does.
4. **CCM89 with R_V = 3.1.**

**What it supports.** Two measured line ratios and, under those assumptions, two
extinctions: R02 at A_V ≈ 0.4 mag, ordinary for a star-forming galaxy; R01 at
A_V ≈ 1.9 mag, which is high.

**What it does not establish.** That R01 is reddened by 1.9 mag. The same
decrement is produced by an underestimated Hβ, and Hβ is the quantity §2 shows to
be model-dependent at the 20–35 % level — a 20 % low Hβ moves A_V from 1.92 to
2.4. **The honest reading of R01 is "the decrement is inconsistent with case B by
3.3σ", which is a statement about the measurement; the extinction is the
interpretation.** A 20 % low Hβ — well inside the ±14–20 % spread of §2 — moves
A_V from 1.92 to **2.62**. Only two objects have a decrement at all, so no extinction law,
no reddening distribution and no population statement can be derived.

---

## 6. What was deliberately not plotted

* **A BPT diagram of all six objects.** Four have no measurable Hβ; drawing them
  at a nominal position would invent a measurement.
* **A colour–magnitude or mass–metallicity relation.** No photometry, no stellar
  masses, n = 6.
* **An extinction-corrected spectrum.** The correction needs a measured Hβ, which
  exists for two objects, and the correction factor at that A_V is 1.2–3× — the
  corrected spectrum would be dominated by the uncertainty of the correction.
* **A fit-quality "score".** A single number would hide the structure §1–§4 show.

## 7. How the figures are kept honest

| rule | enforced by |
|---|---|
| every plotted value comes from code, not a caption | the scripts compute; nothing is hard-coded |
| the decomposition is the fitted model | `test_curves_reproduce_the_fitted_chi2` |
| a non-detection is never a point | `test_nondetections_are_not_drawn_as_measurements` |
| the extinction curve is the published one | `test_ccm89_matches_published_r_lambda` |
| figures are reproducible byte for byte | `test_figures_are_byte_reproducible` |
| inputs and versions are recorded | `FIGURE_MANIFEST.json` from `make_all.py` |
