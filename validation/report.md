# Validation report

What was measured, on what, and what it does and does not establish.
Every number here was produced on the pinned environment
(python 3.12.3, numpy 1.26.4, scipy 1.11.4, astropy 7.0.0) and is reproducible
with the commands in §7.

---

## 1. The baseline

`validation/results/golden_reference.csv` — 80 rows, produced by the **frozen
pre-migration implementation** and, independently, by `gspec`. sha256:

```
460913a261a8ebecd88b89d63266c87d0f1ad906c71ff2fa8907d4fb76a2211d
```

It is byte-identical to two independent runs of the frozen code (two runs, same
host, one hash), and `gspec` reproduces every printed value. Machine-readable
metadata — environment, per-case timings, accuracy, redshift residuals, flag
counts — is in `golden_reference.meta.json`.

## 2. Accuracy against analytic truth (10 synthetic spectra)

| line | median error | worst | acceptance band | rows inside it alone |
|---|---|---|---|---|
| Hα | **−0.70 %** | +41.6 % | 7 % | 8 / 10 |
| Hβ | **−14.42 %** | −94.7 % | 20 % | 7 / 10 |
| [O III] 5007 | **+0.31 %** | −65.5 % | 6 % | 7 / 10 |
| [N II] 6583 | **−1.54 %** | −7.7 % | 10 % | 10 / 10 |

**Read the worst column carefully.** It is not a failure: on those rows the
reported σ is itself large, and the acceptance rule `max(band, 3σ)` is carried by
the uncertainty term. S09 Hβ is −94.7 % with σ = 96 % of the truth. That is the
honest-uncertainty design working — but it also means "passes" is weaker than
"accurate", which is why the medians are pinned by a test.

Redshift: median |Δv| **1.7 km/s**, max **18.5 km/s**.

## 3. Differential recovery on real spectra (4 injection pairs)

The strongest test in the repository, because the truth cannot come from any
measurement of the data: the copy and its base share one noise realisation, so
`F(copy) − F(base) = F_injected` exactly, whatever the continuum does.

| pair | Hα | Hβ | [O III] |
|---|---|---|---|
| I01 ← R01 | −1.4 % | +2.7 % | +0.3 % |
| I02 ← R02 | −0.1 % | +3.5 % | −1.2 % |
| I03 ← R03 | −2.3 % | −19.4 % | −2.4 % |
| I04 ← R04 | −3.4 % | −18.7 % | −11.5 % |

12 of 12 comparisons pass the acceptance rule. The Hβ misses on I03/I04 are at
the level of the per-measurement error on increments of 8 and 12 units.

## 4. Uncertainty calibration

| quantity | measured | accepted range |
|---|---|---|
| pull RMS, 40 synthetic line fluxes | **1.232** | 0.5 – 2.5 |
| pull RMS, pooled 40 + 12 increments | 1.084 | 0.5 – 2.5 (uninformative) |
| pulls beyond 5σ | 0 % | ≤ 20 % |
| median σ/F | **0.052** | — (reported) |
| max σ / true flux | ≤ 1.18 | ≤ 1.5 |

The pooled figure is quoted for continuity with the historical record and is
**not** a calibration measure: the twelve injection increments carry an error
that is deliberately conservative (see §3), their pulls have RMS 0.195, and
mixing them with the flux pulls lowers the aggregate to 1.084. The quantity the
calibration claim rests on is the 40-row one.

## 5. Independent comparison — `gspec` against a separately written program

Sections 1–4 compare the pipeline with *itself* (a frozen copy) and with truth.
This section compares it with `validation/independent/measure_independent2.py`:
an RBF continuum, a bounded optimiser and a bootstrap, sharing no estimator with
`gspec`. Agreement here is evidence rather than a regression check. The default
estimator was **not** changed to improve agreement: agreement is evidence, not an
objective.

### 5.1 Redshift — independently confirmed

| statistic | value |
|---|---|
| median \|Δv\| between the two methods | **6.5 km/s** |
| worst | **34.8 km/s** (R05, the broad-line AGN) |
| cases agreeing within 60 km/s (one instrumental σ) | **20 / 20** |
| cases agreeing within 150 km/s | **20 / 20** |

Two entirely different estimators — a matched filter with an NNLS polish, and a
cross-correlation with a per-window bounded fit — land on the same velocity to a
tenth of the resolution element, on every case including the six real spectra.
**This is the strongest independent result in the repository.**

### 5.2 Flux agreement — good for Hα and [N II], poor for Hβ

| line | median difference | mean | worst | worst case |
|---|---|---|---|---|
| Hα | **+0.8 %** | +5.3 % | +29.4 % | S06 |
| Hβ | **−17.5 %** | +6.1 % | +216 % | R06 |
| [O III] 5007 | +6.3 % | +81.5 % | +1282 % | R06 |
| [N II] 6583 | **+0.2 %** | −1.2 % | −17.6 % | S09 |

Hα and [N II] agree to under a percent at the median, across 20 spectra and two
independent estimators. [O III] agrees on bright rows; its mean is dominated by
two non-detections where both values are noise (R06: 1.63 vs 0.12).

**Hβ disagrees systematically, and the two methods bracket the truth:**

| case | gspec | independent | gap | true | gspec/true | independent/true |
|---|---|---|---|---|---|---|
| S01 | 312.65 | 321.86 | −2.9 % | 315.79 | −1.0 % | +1.9 % |
| S02 | 123.44 | 162.26 | −23.9 % | 135.48 | −8.9 % | +19.8 % |
| S03 | 52.74 | 92.12 | −42.7 % | 76.47 | −31.0 % | +20.5 % |
| S04 | 128.64 | 245.96 | −47.7 % | 157.89 | −18.5 % | +55.8 % |
| S05 | 159.69 | 194.99 | −18.1 % | 185.71 | −14.0 % | +5.0 % |
| S06 | 47.24 | 108.25 | −56.4 % | 50.00 | −5.5 % | +116.5 % |
| S07 | 147.11 | 178.65 | −17.7 % | 172.73 | −14.8 % | +3.4 % |
| S08 | 80.40 | 119.93 | −33.0 % | 103.03 | −22.0 % | +16.4 % |
| S09 | 1.05 | 42.80 | −97.5 % | 20.00 | −94.7 % | +114.0 % |
| S10 | 132.18 | 148.46 | −11.0 % | 135.00 | −2.1 % | +10.0 % |

**Median gap −28.4 %, and the truth lies between the two methods in 10 of 10
cases.** That is the single most informative number here: the Hβ disagreement is
not random scatter between two noisy estimators, it is a systematic that puts one
method consistently low and the other consistently high. Neither is right, and
averaging them would not produce a measured value — it would produce a number with
no error model (`docs/limitations.md` L2).

The mechanism is identified in `validation/experiments/FINDINGS.md`: `gspec` fits
no stellar absorption on these spectra (BIC rejects the term), so Hβ emission is
partly eaten by unmodelled absorption and recovered low; the independent
implementation fits the continuum jointly with an RBF basis and recovers it high.
Forcing absorption on in `gspec` overshoots in the other direction (+19.6 %).

### 5.3 Uncertainties — calibrated in aggregate, different per row

| estimator | n | pull RMS | median \|pull\| | max \|pull\| | beyond 5σ |
|---|---|---|---|---|---|
| gspec | 40 | **1.232** | 0.958 | 2.50 | 0.0 % |
| independent | 40 | **1.306** | 0.712 | 2.99 | 0.0 % |
| gspec, pooled with increments | 52 | 1.084 | 0.655 | 2.50 | 0.0 % |
| independent, pooled with increments | 52 | 1.151 | 0.440 | 2.99 | 0.0 % |

`n = 40` is the ten synthetic spectra × four lines, where the truth is analytic
and the error is the estimator's own. The pooled rows add the twelve injection
increments, whose uncertainties are conservative by construction; they are shown
only because earlier versions of these documents quoted them.

Both are calibrated; neither has a pathological tail. But their *magnitudes* differ
per row — median σ ratio (gspec / independent) is 1.48 for Hα, 1.64 for Hβ, 1.10
for [O III], 0.72 for [N II], with individual rows ranging from 0.12 to 119.

**The honest statement about the uncertainties:** aggregate calibration is
independently reproduced; the per-row magnitude is an estimator-dependent
quantity, and no claim should be made about a single row's σ to better than a
factor of a few.

### 5.4 Detection flags — 71 / 80 agree, and the 9 conflicts are one-directional

| gspec | independent | rows |
|---|---|---|
| ok | ok | 66 |
| nondetection | ok | **9** |
| nondetection | nondetection | 5 |

Every conflict is `gspec` calling a line a non-detection where the independent
method calls it `ok` — never the reverse. Eight of the nine are **Hβ**, the ninth
is R06 Hα. The cause is visible in the significances:

```
I03 Hbeta   |F|/sigma  gspec  1.41   independent  4.78
I04 Hbeta              gspec  0.52   independent 11.99
R04 Hbeta              gspec  0.48   independent 12.91
S06 Hbeta              gspec  0.80   independent  5.55
S09 Hbeta              gspec  0.05   independent  3.95
```

`gspec` reports a larger σ on these rows (§5.3), so its significance falls below
the 3σ line. Because the flag is *defined* as consistency with the reported
numbers, both flags are internally correct — the disagreement is about σ, not
about the flux. And σ on exactly this shape of row is the documented fragility
(`docs/limitations.md` L4): a pinned broad-Hβ amplitude carries most of the
variance.

### 5.5 Verdict

| claim | status after the independent comparison |
|---|---|
| redshift | **independently confirmed** (20/20 within one instrumental σ) |
| Hα flux | **independently confirmed** (median 0.8 %) |
| [N II] flux | **independently confirmed** (median 0.2 %) |
| [O III] flux | confirmed on bright rows; uninformative on non-detections |
| Hβ flux | **NOT confirmed** — 17.5 % systematic gap, truth bracketed 10/10 |
| σ calibration | independently confirmed in aggregate; per-row magnitude is estimator-dependent |
| flags | 71/80; conflicts traceable to σ, and all in one direction |

### 5.6 How reproducible is the independent implementation itself?

It is **not** bit-reproducible, and this repository previously claimed it was.
Re-running `validation/independent/run_independent.py` on the 20 shipped spectra
(27 min on this host) against the shipped `out_independent.csv`:

| quantity | result |
|---|---|
| redshifts | **identical, 20 / 20** |
| flux, median difference | 0.03 % Hα · 0.68 % Hβ · 0.15 % [O III] · 0.03 % [N II] |
| flux, worst row | **14.2 %** (S04 Hβ: 245.96 → 211.07) |
| rows moving more than 1 % | 4/20 Hα, **8/20 Hβ**, 6/20 [O III], 3/20 [N II] |

The program resamples, so it carries a stochastic component. What matters is that
**the comparison in §5.2 is stable across both runs**: the synthetic Hβ median gap
is −28.4 % against the shipped output and −29.3 % against the re-run, and the
analytic truth falls between the two methods in **10 of 10** cases either way. The
systematic disagreement is therefore an order of magnitude larger than the
method's own run-to-run noise (median 0.68 %).

The practical consequence, now stated in `validation/independent/README.md`: a
single row of that output should be quoted with a per-row tolerance of ~15 %, not
as an exact value.

## 6. The negative control

`m2` — the same fitter with no stellar-absorption term — must be rejected.

| check | m1a | m2 |
|---|---|---|
| injection differential | 12/12 pass | **fails**: I04 loses all three lines |
| I04 redshift | 0.08258 | **0.01950** (18 910 km/s off) |
| archival redshift, R04/R06 | within 300 km/s | 18 910 / 14 885 km/s off |

Note what the control does **not** fail: its synthetic Hβ median is −11.1 %,
better than the reference's −14.4 %. The control is therefore tested against the
property it demonstrably breaks (the differential), not against an accuracy claim
it happens to survive. See `docs/limitations.md` L3.

## 7. Reproducing

```bash
python -m venv --system-site-packages .venv && source .venv/bin/activate
pip install -e ".[dev]"

pytest                                                    # 105 tests, ~240 s
gspec measure --cases data/cases --out out/line_measurements.csv
diff <(sort out/line_measurements.csv) \
     <(sort validation/results/golden_reference.csv)       # identical

python validation/experiments/e_abs_absorption.py          # ~2 min
python validation/experiments/e_cov_bound_pinned.py        # ~2 min

# the independent implementation used in §5: ~15 min on this host (46 s/case)
python validation/independent/run_independent.py \
       --cases data/cases --out /tmp/out_independent_reproduced.csv
```

## 8. Reproducibility protocol, and its limits

Four tiers, each answering a different question:

| tier | comparison | tolerance | result |
|---|---|---|---|
| 1 | same environment, two runs | bit-identical | ✓ two frozen runs, one sha256 |
| 2 | ported code vs frozen code, same host | rtol 1e-12 | ✓ **0.0** on all 80 rows × 4 fields |
| 3 | ported code vs the golden CSV | rtol 1e-6 | ✓ identical to the printed precision |
| 4 | across environments | no equality claim | differs — see below |

**Tier 4, measured.** The same code in a `python:3.12-slim-bookworm` container
produces the same fluxes to ≤ 6.4e-4 relative but a *different* σ on four rows
(S04 Hα 14.62 vs 50.68, S04 Hβ 16.82 vs 41.01, S08 Hα 11.50 vs 40.15, S08 Hβ
19.70 vs 51.98) and one flag flips. So:

> This repository is reproducible to the precision of its input (≈1e-6 relative)
> under the pinned environment, **with two documented exceptions**: σ on
> bound-pinned amplitudes (L4) and discrete flag/model selections.

It is **not** bit-reproducible across environments, and no claim in this
repository depends on that.

## 9. What this validation does not establish

| claim | status |
|---|---|
| "the method is competitive with template fitting" | **NOT TESTED** — no pPXF comparison exists |
| "Hβ is accurate to 20 %" | **NOT DEMONSTRATED** — median −14.4 % for the reference and +18.1 % for the independent method, a 28 % disagreement between them |
| "the absorption model improves Hβ" | **CONTRADICTED** — §6, `experiments/FINDINGS.md` |
| "real-data fluxes are scientifically validated" | **PARTLY** — 4 injections and a factor-2 band |
| "the results generalise to other surveys" | **NOT TESTED** — one instrument, one reduction |
| "the sample supports a population statement" | **FALSE** — n = 20 |

## 10. Provenance of this validation

The scientific core, the datasets and the ground truths come from an earlier
project of the same author, where the pipeline was developed and where its
acceptance criteria were attacked by an independent adversarial review that found
and closed five real vulnerabilities (docs/validation-history.md §2). The
scientific findings were carried forward; the audit scaffolding was not, and
nothing in this repository scores a third-party measurement. What survives is the
set of requirements the audit produced, now expressed as tests in `tests/`.
