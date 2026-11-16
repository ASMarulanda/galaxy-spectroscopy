# Experimental findings

Two controlled experiments, run once, recorded here. **Neither is part of the
default pipeline**: the package's behaviour is unchanged by both, and
`validation/results/golden_reference.csv` still reproduces bit for bit.

Each experiment exists because a claim in the documentation needed testing
rather than assertion. Both produced results that contradicted the expectation
going in — which is the point of running them.

---

## E-Abs — does the stellar-absorption term do anything?

Script: `e_abs_absorption.py` · Data: `results/e_abs_absorption.json` · 10 synthetic cases

### Design

Four configurations, all fitted at the **same redshift seed** (the value the
default pipeline reports), so the only thing that varies is the model:

| config | model |
|---|---|
| `bic` | the shipped pipeline: 3-line absorption on the table, BIC decides |
| `abs_on` | 3-line absorption forced on, `ndeg=3`, no selection |
| `abs_off` | no absorption term, `ndeg=3`, no selection |
| `abs_gen` | absorption forced on with the **generator's** relative weights (Hγ 1.00, G 0.45, Hβ 1.00 instead of 0.95 / 0.60 / 1.00) |

### Result 1 — forcing absorption on makes Hβ substantially worse

Median error per line, over the ten cases (%):

| config | Hα | Hβ | [O III] | [N II] | mean abs. error |
|---|---|---|---|---|---|
| `bic` (shipped) | −0.70 | **−14.42** | +0.31 | −1.54 | 11.52 |
| `abs_on` | −0.75 | **+19.63** | +0.02 | −1.01 | **31.27** |
| `abs_off` | −0.75 | **−11.13** | +0.35 | −1.01 | **14.40** |
| `abs_gen` | −0.75 | +16.96 | +0.02 | −1.01 | 30.66 |

Worst-case |error| on Hβ: 94.7 % (shipped), **462.1 %** (`abs_on`), 154.6 % (`abs_off`).

The sign flips. Without absorption the fitter under-estimates Hβ by ~11–14 %;
with absorption forced on it over-estimates by ~17–20 %. The true flux is
bracketed by the two configurations, and neither is right.

### Result 2 — the basis is not the problem

`abs_gen` (the generator's own weights) differs from `abs_on` by 2.7 pp on the
Hβ median (19.63 % vs 16.96 %). Weight mismatch is therefore *not* what makes
the term harmful. A richer basis would not rescue it.

### Result 3 — BIC prefers absorption on 1 of 10 cases

ΔBIC = `abs_on` − `abs_off` (negative means absorption wins):

```
S01  +1.8   S02  -2.4   S03  +9.4   S04  +1.5   S05 +16.6
S06 +13.7   S07 +13.1   S08 +19.8   S09  +8.3   S10 +11.5
```

χ² always improves (−0.9 to −23.2) but never enough to pay for three
parameters (the price is 3·ln(1027) ≈ 20.8). Only S02 clears it, by 2.4 units —
a marginal preference, not evidence.

### Result 4 — the fit does not recover the absorption that is there

With absorption forced on, the fitted depth against the generator's true depth:

```
case  true abs_scale  fitted  ratio      true σ*/σ_gas  fitted σ*/σ_gas
S01        0.670       0.054   0.08            2.93          3.64
S02        0.927       0.146   0.16            2.04          2.05
S03        0.860       0.158   0.18            2.47          1.90
S04        1.322       0.121   0.09            1.35          1.54
S05        0.657       0.273   0.42            2.47          1.50
S06        1.087       0.838   0.77            1.88          1.15
S07        0.944       0.008   0.01            1.35          4.52
S08        0.785       0.286   0.36            1.98          1.15
S09        1.303       0.425   0.33            1.61          1.15
S10        1.165       0.540   0.46            1.35          1.15
```

The depth is systematically under-recovered (ratio 0.01–0.77, median ≈ 0.22), and
the width ratio lands **exactly on its lower bound, 1.15, in four cases** — the
prior, not the data, is setting it. The term is not identifiable on these spectra
with this basis.

### Verdict

**The absorption model is retained, and the claim that it improves Hβ is
withdrawn.** The permitted statements are:

* the fitter contains an explicit multiplicative Balmer absorption model, and its
  BIC selection keeps that term on 5 of 6 real spectra;
* on synthetic spectra with analytic truth the term is *not* supported: it is
  rejected by BIC on 9 of 10 cases and degrades Hβ when forced on.

The Hβ bias is therefore best described as a **model-misspecification signature**:
the fitter has a 3-line approximation of a 12-line effect, and the polynomial
continuum absorbs the smooth part. Neither including it nor excluding it gives a
correct Hβ. This is limitations L2/L3 and it is the most interesting negative
result in the repository.

---

## E-Cov — what happens to σ_F when a parameter sits on its bound?

Script: `e_cov_bound_pinned.py` · Data: `results/e_cov_bound_pinned.json` · all 20 cases

### Design

`cov = pinv(JᵀJ) · min(χ²/dof, 1.5)` is the shipped estimator and is **not**
modified. Two estimators are compared:

* **shipped** — pseudoinverse over all parameters, exactly as the package reports;
* **conditioned** — the parameters pinned at a bound are removed from the problem
  before the covariance is formed, i.e. the uncertainty is computed *given* the
  model holds them fixed.

### Result 1 — the pinned amplitude dominates σ, and removing it shrinks σ a lot

27 of 80 rows move by more than 10 %. The largest:

| case | line | shipped σ | conditioned σ | ratio | dominant direction | share of σ² | λ/cutoff |
|---|---|---|---|---|---|---|---|
| R03 | Hβ | 14.81 | 2.57 | **5.77** | `amp_Hb_b` | 77.3 % | 1.3e6 |
| R05 | Hβ | 690.05 | 149.42 | **4.62** | `amp_Hb_n` | 99.7 % | 1.1e2 |
| I03 | Hβ | 14.84 | 3.27 | 4.54 | `amp_Hb_b` | 73.4 % | 8.6e5 |
| S09 | Hβ | 19.20 | 4.51 | 4.25 | `amp_Hb_b` | 90.9 % | 9.9e6 |
| S06 | Hβ | 58.97 | 15.09 | 3.91 | `amp_Hb_b` | 77.9 % | 9.0e5 |
| S01 | Hβ | 1.384 | 0.518 | 2.67 | `amp_Hb_b` | 80.7 % | 8.4e3 |
| S08 | Hβ | 19.70 | 14.98 | 1.32 | `narrowB` | 48.2 % | 3.9e6 |

On S01, 81 % of σ_F(Hβ)² comes from the amplitude of the **broad** component — a
parameter the BIC selection has effectively zeroed. That is the structural fact
behind the whole question.

### Result 2 — conditioning would flip six flags from `nondetection` to `ok`

```
S06 Hbeta, R03 Hbeta, R04 Hbeta, R06 Halpha, I03 Hbeta, I04 Hbeta
```

All six are lines whose significance sits just under 3σ and crosses it when σ
shrinks.

### Result 3 — but conditioning degrades the calibration

Against analytic truth on the 40 synthetic rows:

| estimator | pull RMS | median \|pull\| | max \|pull\| |
|---|---|---|---|
| shipped | **1.231** | 0.928 | **2.50** |
| conditioned | 1.739 | 0.968 | **6.05** |

The conditioned estimator is over-confident: its worst pull goes past the 5σ hard
limit that the validated calibration sets, and its RMS moves from "good" to
"borderline" inside the accepted band [0.5, 2.5].

### Verdict

**Keep the shipped estimator.** The pinned-amplitude direction carries variance
the data genuinely do not constrain, so removing it manufactures confidence.
Dropping bound-pinned directions was the hypothesised improvement
(`HANDOVER/58` §F.4, option b); the measurement rejects it.

The fragility recorded in `docs/limitations.md` L4 is therefore the *price* of a
correctly calibrated σ, not a bug to be optimised away. What remains open is the
cross-environment instability of that one direction (S08 Hβ, 19.70 vs 51.98): a
future estimator would have to be **both** environment-stable and calibrated,
and conditioning alone achieves neither.

---

## Reproducing

```bash
.venv/bin/python validation/experiments/e_abs_absorption.py     # ~2 min
.venv/bin/python validation/experiments/e_cov_bound_pinned.py   # ~2 min
```

Both are importable, print their tables, and write JSON into `results/`. They are
marked `experiment` in the test configuration and are never run by the test suite.
