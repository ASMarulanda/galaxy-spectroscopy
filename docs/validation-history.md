# Validation history

How these measurements came to be trusted — and what it cost.

This is a condensed record of an adversarial audit of the *validation procedure
itself*, run before the code was restructured into a public package. It is kept
because the process, not just the result, is part of the claim: the acceptance
criteria were deliberately attacked, the attacks succeeded, and the defences are
what you see in `tests/` today.

---

## 1. The four independent ground truths

Before any attack, the design fixed how truth is established, and deliberately did
not use one mechanism:

| mechanism | n | why it cannot leak from the method |
|---|---|---|
| analytic synthetic | 10 | the generator writes each Gaussian so its **integral** is the declared flux; a self-test re-integrates on the pixel grid to 8.8e-9 |
| exact injection | 4 | copy and base share **one noise realisation**, so the difference of two measurements *is* the increment — the continuum cancels |
| archival | 6 | external catalogue values, used only as a broad consistency band (±300 km/s, factor 2) |
| independent implementation | 1 program | written separately: RBF continuum, bounded optimiser, bootstrap; shares no estimator with this package |

## 2. The adversarial audit — five findings, all real

An independent reviewer was asked to break the acceptance criteria with no access
to the fitter. It broke them five times. In the worst demonstration, 30 of 40
synthetic fluxes were reported at **20× the true value**, injected increments were
delivered at 50 %, and the pass/fail summary still read "0 failures" on both
accuracy checks.

The root cause was the same in every case: the tolerance was
`max(floor·truth, 3·reported_σ)`, and **nothing bounded `reported_σ`**. An
implementation can always place its own residual inside 3σ by reporting
`σ = residual / 3`. The criterion reads as a scientific statement about accuracy
and behaves as a statement about the reporter's honesty.

Five findings, and what closed each:

**F1 (critical) — the flux accuracy check was purchasable with a fabricated σ.**
Closed by a per-row ceiling: no reported σ may exceed 1.5× the true line flux.
Both valid implementations stay below 1.18, and only on faint non-detections.

**F2 (high) — the aggregate pull statistic could be diluted.** A subset of
well-chosen pulls could sit inside the accepted RMS band while a third of the
rows were wrong. Closed by making the aggregate two-sided *and* independent of σ:
at most 6 of 40 rows may sit outside 3× their relative floor, with failed rows
counted in that budget because they carry no measurement either.

**F3 (critical) — the injection differential was purchasable.** Both rows of a
comparison came from the same implementation, so it set the tolerance for its own
differential. Closed by capping the tolerance at 2× the relative floor: the
increment is known exactly, so a 50 % miss cannot be excused by a large claimed
uncertainty.

**F4 — a computed statistic was never tested.** `median σ/F` cleanly separated
every counterexample (6.55) from both valid implementations (0.07) and had been
discarded. It is now part of the uncertainty contract.

**F5 — an archival escape, found later.** Two faint rows could be pushed to 9000
with the σ escape on the real-data comparison. Closed by a flag-independent
magnitude cap: no row may differ from the published catalogue by more than 10×,
whatever its σ or its flag.

## 3. What the defences cost, and why some were refused

Three proposed defences were **measured and rejected** rather than adopted,
because they could not separate a wrong implementation from a right one:

* a σ-inflation bound of 10× — the reference's own worst row is 19.4×;
* a hidden photon-significance threshold for "detected" — the achievable
  precision here is model-limited, so the threshold punished honest
  non-detections;
* a lenient allowance for non-detections — it made the flag change which rows
  were graded, which is exactly the coupling the contract forbids.

Two deliberate, *stated* allowances stood in the original procedure and are
reproduced here as facts rather than hidden: six of forty synthetic rows may sit
up to ~4.35× from the truth, and two of twenty-two archival rows may sit outside
the band. They existed because the reference implementation itself needed them.
This repository does not enforce any of it — it reports the same weakness as
limitation L2 in `docs/limitations.md`.

## 4. The negative control

`m2` (in `src/gspec/variants.py`) removes the stellar-absorption term and changes
nothing else. It was constructed to be **rejected**, and it is: it loses the I04
injection entirely and collapses that object's redshift to z = 0.0195 instead of
0.0826 — an 18 910 km s⁻¹ error.

What it does *not* lose is synthetic accuracy — and that observation is what led,
two stages later, to the finding that the absorption term is not supported by the
data at all (`docs/limitations.md` L3, `validation/experiments/FINDINGS.md`). A
control built to fail ended up falsifying a claim in the documentation.

## 5. What the audit taught, in one paragraph

Every vulnerability had the same structure: a rule that looked quantitative was
actually negotiable by the party being measured. The remedy in each case was a
quantity the implementation cannot choose — an exact increment, a truth-derived
floor, a catalogue magnitude. That principle survives the move to a public
library: the tests in `tests/` are written against analytic truth and exact
increments, not against the pipeline's own output, and the one test that could
have become circular (`test_regression`) compares against a baseline produced by
*frozen* code.

## 6. Chronology of the science

| step | what it established |
|---|---|
| units and masks verified against the archive | `LINEAREA/(LINECONTLEVEL×LINEEW) = 1.0000` over 32 measurements, so the catalogue flux scale is integrated flux, not density |
| three catalogue conventions compared | pipeline / MPA-JHU / Portsmouth disagree by a median of 139 % per line — hence a band, not an accuracy test |
| the injection method validated | recovery −0.1 % to +0.1 % on well-modelled objects |
| the absorption width tied to the gas width | removed a measured +10…+80 % Hβ bias from an unidentifiable width |
| the broad component bounded | an unbounded 2500 km/s component acted as a pseudo-continuum (44 % amplitude swing for a 16 % perturbation) |
| the covariance cap introduced | uncapped `χ²/dof` scaling over-corrected on real spectra |
| the independent implementation built | a second, estimator-disjoint method reaching the same acceptance |

The order matters: each fix was made because an experiment showed the previous
choice failing, and the later absorption result (§4) shows that at least one of
them — the absorption model itself — fixed the symptom and not the cause.
