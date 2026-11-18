# Decision log

Who decided what, on what evidence. Kept because a scientific codebase is a
sequence of choices, and a reader deserves to know which ones were measured,
which were assumed, and which were made by whom.

**Attribution keys**

| key | meaning |
|---|---|
| **[ASMD]** | the author decided this, explicitly |
| **[HERMES]** | the assistant proposed it and implemented it; approved by the author |
| **[MEASURED]** | the evidence is a number produced by a run recorded in this repository |
| **[ASSUMED]** | a modelling choice with a physical justification but no separating experiment |

---

## Phase 1 — data selection and units

| # | decision | who | evidence |
|---|---|---|---|
| 1.1 | use SDSS DR18 spectra of the SDSS-I/II Legacy survey, `run2d=26` | [HERMES] [MEASURED] | 10 files downloaded twice, byte-identical |
| 1.2 | select ten galaxies spanning S/N 7–39, z 0.021–0.087 and emission-line strength | [HERMES] [MEASURED] | `data/provenance/sdss_source_picks.csv` |
| 1.3 | treat `SPZLINE.LINEAREA` as integrated flux in 1e-17 erg s⁻¹ cm⁻² | [HERMES] [MEASURED] | `LINEAREA/(LINECONTLEVEL×LINEEW) = 1.0000` over 32 measurements |
| 1.4 | compare three catalogue conventions rather than trusting one | [HERMES] [MEASURED] | median disagreement 139 % per line |
| 1.5 | do **not** investigate the Portsmouth flux-scale anomaly further; exclude those fluxes | **[ASMD]** | an unresolved scale discrepancy (Hα ≈ 2.4× high) is recorded, not explained |

## Phase 2 — the contract

| # | decision | who | evidence |
|---|---|---|---|
| 2.1 | analytic synthetic truth + exact injection are the primary verification; archival data only for consistency | **[ASMD]** | reference-solution outputs are kept out of the validation path |
| 2.2 | Hα means **total** narrow + broad integrated flux | **[ASMD]** | stated as a definition, not a measurement |
| 2.3 | synthetic continuum = three power laws × absorption, deliberately *not* the fitter's polynomial basis | [HERMES] [ASSUMED] | prevents passing by assuming the generator's model |
| 2.4 | absorption width tied to the gas width, ratio ≥ 1.15 | [HERMES] [MEASURED] | removed a +10…+80 % Hβ bias from an unidentifiable width |
| 2.5 | absorption depth in the generator fixed at 0.11 | [HERMES] [ASSUMED] | deeper values re-created the degeneracy; a fairness choice that also sets the difficulty of the hardest line |
| 2.6 | broad Hα/Hβ widths tied and bounded to [250, 1800] km/s | [HERMES] [MEASURED] | an unbounded 2500 km/s component acted as a pseudo-continuum |
| 2.7 | covariance scaled by `min(χ²/dof, 1.5)` | [HERMES] [MEASURED] | uncapped scaling over-corrected on real spectra |
| 2.8 | BIC model selection over seven nested models | [HERMES] [MEASURED] | a broad component or a cubic absorbs a weak line at low S/N |
| 2.9 | keep S06 at S/N 7 as an honest marginal case | [HERMES] [MEASURED] | it passes; removing it would have made the set easier |
| 2.10 | R05 has no injection copy | [HERMES] [MEASURED] | +55 %, +99 %, +100 % recovery under three configurations — not reproducible |

## Phase 3 — validation design

| # | decision | who | evidence |
|---|---|---|---|
| 3.1 | an independent, code-disjoint implementation is required | **[ASMD]** (via 2.1) | `validation/independent/` shares no estimator with the package |
| 3.2 | tolerances from measurement, never from the reference's output | [HERMES] [MEASURED] | `validation/results/tables/floor_margins.csv` |
| 3.3 | Hβ floor 20 %, and its justification re-derived on independent evidence | [HERMES], approved **[ASMD]** | the earlier "inter-method spread" argument was withdrawn because both methods shared a fitter |
| 3.4 | treat flags as consistency statements about the reported numbers, never as quality judgements | [HERMES] [MEASURED] | a hidden significance threshold punished honest non-detections |
| 3.5 | a negative control (`m2`) is a required part of the suite | [HERMES] [MEASURED] | it is rejected — and it later falsified the absorption claim |

## Phase 4 — the adversarial audit (see validation-history.md)

| # | decision | who | evidence |
|---|---|---|---|
| 4.1 | bound the reported σ per row (≤ 1.5 × true flux) | [HERMES] [MEASURED] | closed a 20× purchase that scored full marks |
| 4.2 | at most 6 of 40 rows may sit outside 3× their floor, failed rows included | [HERMES] [MEASURED] | closed the dilution of the aggregate pull statistic |
| 4.3 | cap the differential tolerance at 2× the relative floor | [HERMES] [MEASURED] | closed 50 % injection delivery |
| 4.4 | cap every archival row at 10× the catalogue magnitude, flag-independent | [HERMES] [MEASURED] | closed the last σ-priced escape |
| 4.5 | refuse three other defences that could not separate right from wrong | [HERMES] [MEASURED] | reference's own worst rows exceed the proposed bounds |

## Phase 5 — migration to a public package

| # | decision | who | evidence |
|---|---|---|---|
| 5.1 | repository `galaxy-spectroscopy`, package `gspec` | **[ASMD]** | — |
| 5.2 | ship the 20 case FITS plus the truth manifests; withhold `truth_real.json` pending provenance review | **[ASMD]** | it maps each case to plate/MJD/fibre |
| 5.3 | MIT for code, CC BY 4.0 for text and figures, SDSS provenance kept separate | **[ASMD]** | `DATA.md` |
| 5.4 | publish a condensed validation history; keep the full audit record private | **[ASMD]** | `docs/validation-history.md` |
| 5.5 | do not claim the absorption model improves Hβ; test it instead | **[ASMD]** | `validation/experiments/FINDINGS.md` — the claim did not survive |
| 5.6 | keep the original covariance estimator; investigate bound-pinned directions separately | **[ASMD]** | conditioning was tested and **rejected**: it degrades calibration (pull RMS 1.23 → 1.74) |
| 5.7 | include `m1b` as an explicitly labelled *alternative*, never as independent | **[ASMD]** | it shares the fitter and differs in three constants |
| 5.8 | port structurally, verify against a frozen baseline before refactoring | **[ASMD]** | the port reproduces all 80 rows to 0.0 |
| 5.9 | expose diagnostics (χ², dof, n_par, selected model, bound state) on every result | **[ASMD]** | the selected model is itself a scientific result |
| 5.10 | defer notebooks and machine learning to later stages | **[ASMD]** | — |

---

## Phase 6 — retrieval workflow (v0.3)

| # | decision | who | evidence |
|---|---|---|---|
| 6.1 | identifiers only; no catalogue query, so no selection function to declare | **[ASMD]** | design gate, HANDOVER/62 §F.8 |
| 6.2 | keep network retrieval **outside** the validated package; no `gspec fetch` subcommand | **[ASMD]** | the HTTP layer would enter the module tree the numbers came from, for a convenience the script already gives |
| 6.3 | the converter refuses an incompatible file rather than coercing it | [HERMES] | twelve contract checks; each has a test |
| 6.4 | settle `wdisp`'s unit from the data (pixels vs ångström), not from memory | [HERMES] [MEASURED] | it falls 1.258→0.696 across the band while Δλ doubles; R 1466→2650 matches the Legacy curve |
| 6.5 | correct the archive URL to the zero-padded plate | [HERMES] [MEASURED] | unpadded returns 404, padded returns 200 |
| 6.6 | **disclose the R06 identity** (plate 0294, MJD 51986, fibre 0036) | **[ASMD]** | it was already derivable from `sdss_source_picks.csv` + `archival_redshift.csv` (5 of 6, uniquely), and v0.3 makes the check one command |
| 6.7 | keep `truth_real.json` internal: five identities and all catalogue values stay unpublished | **[ASMD]** | inspected at the v0.3 review; nothing in it is needed by the workflow |

## Phase 7 — V0.3.1 CI remediation

| # | decision | who | evidence |
|---|---|---|---|
| 7.1 | split the suite into *portable* and *characterisation*; CI runs the portable half | **[ASMD]** | v0.3.0 CI failed 6/195; reproduced exactly on a wheel-based venv |
| 7.2 | tolerances from measurement, not from the failure: redshift 1e-6, flux 1e-3, sigma/error columns 1e-4 | [HERMES] [MEASURED] | worst observed cross-environment flux drift 6.4e-4 (I01 Hbeta, model-sweep case), so 1e-3 is 1.6x above it and 60x below the 6 % accuracy floor; real scientific changes move values by 1e-3..1e-1 |
| 7.3 | sigma is compared only for lines whose complex has no active bound | [HERMES] [MEASURED] | 36 of 80 sigma rows move up to 3.5x; the other 44 stay within 2e-6 |
| 7.4 | remove the exact-baseline hash from CI; keep it as a printed fingerprint | **[ASMD]** | it asserted portability the project never claimed |
| 7.5 | do **not** change the model to remove the flat direction | **[ASMD]** | that would change the numbers; recorded as O6 instead |

## Open questions

| # | question | owner |
|---|---|---|
| O1 | can a template-based continuum (pPXF) make Hβ defensible? | future work |
| O2 | what estimator is *both* environment-stable and calibrated for σ? | future work — both tested candidates fail one of the two |
| O3 | is the 3-line absorption approximation salvageable with a different parameterisation? | future work |
| O4 | broad-line AGN decomposition: untested, and R05 shows why | future work |
| O5 | would a reliability classifier beat the `σ/F < 3` rule? | planned, exploratory |
| O6 | a broad component pinned at zero amplitude keeps a free width, leaving an unconstrained direction that sigma charges for — remove the component, tie the width to the amplitude, or price it out of the covariance? | open, measured in L17; a model change, so not a maintenance-release fix |
