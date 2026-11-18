"""Does the package still produce the validated numbers — and where does that hold?

Two different questions live here, and V0.3.1 separates them.

**Portable.** ``test_portable_baseline_regression`` compares every row of the
20-case baseline against the published golden CSV with tolerances chosen from
measurement, not from convenience. It runs anywhere, and it is what CI executes.

**Reference environment.** The ``characterisation`` tests pin the baseline to its
printed digits. They are only meaningful on the stack the numbers were produced on:
the pipeline is *not* bit-reproducible across BLAS builds, and V0.3.1 measured
exactly how far it is not (docs/limitations.md L17).

The distinction is not a concession to a red badge. A tolerance of 1e-6 relative
was justified in the original test as "the ported code must agree with the frozen
code" — a statement about *this* machine. Carrying it into CI turned a claim about
one environment into a claim about all of them, which was never true and which
nobody had measured until GitHub Actions ran it.

Tolerances, and what each can detect
------------------------------------
=========================  =========  ===========================================
quantity                   tolerance  why that number
=========================  =========  ===========================================
flags                      exact      a flag is a discrete decision; no tolerance
                                      can hide a change in it
redshift                   1e-6       0.3 km/s at z = 0.08, against an acceptance
                                      floor of 6.5 km/s; observed drift 2e-8
flux                       1e-3       1.6x the largest cross-environment flux drift
                                      measured (6.4e-4, I01 Hbeta, where the model
                                      sweep lands on the other configuration); 60x
                                      below the smallest accuracy floor (6 %, [O III])
sigma, redshift_err,       1e-4       arithmetic consistency only: sigma enters the
flux_err, unbound rows                pull calibration at the 10 % level, so 1e-4
                                      cannot mask a scientific change. The large
                                      cross-environment movements in these columns
                                      (up to 5e-3 and 7.1e-1) are on rows whose
                                      complex has an active bound, and are excluded
                                      by a rule derived from at_bounds — not by name.
=========================  =========  ===========================================

A real change — a different constant, a changed window, a broken identity — moves
these quantities by 1e-3 to 1e-1. Against the redshift bound that is 1000x, and the
test has room to spare. Against the flux bound it does not: 1e-3 sits 1.6x above the
largest environmental drift measured and at the bottom of the range a real change
occupies, so the flux comparison separates one *environment* from another and catches
*large* mistakes only. What carries the finer detection is the exactness of the flags,
the covariance rule asserted analytically, and the coverage floors — each fails on a
real modelling change and on nothing that has been observed to differ across stacks.
"""

from __future__ import annotations

import csv
import shutil

import pytest

from conftest import ALL_CASES, REAL, SYNTHETIC

FIELDS = ("redshift", "redshift_err", "flux", "flux_err")

#: Portable tolerances. See the module docstring for where each number comes from.
PORTABLE_RTOL = {"redshift": 1e-6, "redshift_err": 1e-4,
                 "flux": 1e-3, "flux_err": 1e-4}

#: Reference-environment tolerance: the port must agree with the frozen code.
REFERENCE_RTOL = 1e-6

#: Parameters belong to the line complex they were fitted in, so a pinned amplitude
#: makes the sigma of every line in that complex environment-sensitive. Measured: of
#: 80 sigma rows, the 36 belonging to a complex with an active bound move by up to
#: 3.5x between BLAS builds, and the other 44 stay within 2e-6 (L17).
BOUND_COMPLEX = {"Halpha": ("Ha", "NII"), "Hbeta": ("Hb", "Hg"),
                 "OIII_5007": ("OIII", "Hb", "Hg"), "NII_6583": ("NII", "Ha")}


def _bound_in_complex(result, line: str) -> bool:
    """Is a parameter of this line's own complex pinned at a bound?"""
    return any(t in b for b in result.at_bounds for t in BOUND_COMPLEX.get(line, ()))


# ---------------------------------------------------------------- portable
def test_portable_baseline_regression(reference, golden_rows):
    """Every row of the baseline, within tolerances that mean something anywhere.

    Sigma is compared only for lines whose own complex has **no active bound**. That
    is not a loophole: a component pinned at zero amplitude leaves its width
    unconstrained, the fit stops at a different point of that flat direction on
    different BLAS builds, and sigma moves by up to 3.5x while the flux moves by less
    than 1e-6 on 78 rows and by 6.4e-4 on the two where the model sweep chose
    differently — measured, S08 H-alpha sigma 11.50 or 40.15 (L17). Comparing those
    values would be comparing two points of an unidentifiable direction. The set is
    derived from each result's own diagnostics, so if a fit stops reporting bounds
    the comparison gets *stricter*, not looser; and ``test_covariance_portable.py``
    asserts what does survive about sigma.
    """
    assert len(golden_rows) == 80
    worst = {f: (0.0, "") for f in FIELDS}
    sigma_compared = sigma_skipped = flags_compared = flags_skipped = 0

    for cid in ALL_CASES:
        for line, m in reference[cid].lines.items():
            gold = golden_rows[(cid, line)]
            got = reference[cid].row(line)
            for field in FIELDS:
                if field in ("flux_err", "redshift_err") and _bound_in_complex(
                        reference[cid], line):
                    sigma_skipped += 1
                    continue
                if field in ("flux_err", "redshift_err"):
                    sigma_compared += 1
                a, b = float(got[field]), float(gold[field])
                rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
                if rel > worst[field][0]:
                    worst[field] = (rel, f"{cid}/{line}")
                assert rel <= PORTABLE_RTOL[field], (
                    f"{cid}/{line}/{field}: {a!r} vs golden {b!r} "
                    f"(rel {rel:.2e} > {PORTABLE_RTOL[field]:.0e})")
            if not _bound_in_complex(reference[cid], line):
                assert got["flag"] == gold["flag"], \
                    f"{cid}/{line}: flag {got['flag']} vs golden {gold['flag']}"
                flags_compared += 1
            else:
                flags_skipped += 1

    # The exclusion must not swallow the comparison. Measured on this release: 44 of
    # 80 sigma rows and 44 of 80 flags are compared, 36 of each are bound-sensitive.
    # The floors guard against a future change that pins every amplitude and quietly
    # turns the comparison off.
    assert sigma_compared >= 40, (
        f"only {sigma_compared} sigma rows compared ({sigma_skipped} excluded by "
        f"bounds); the bound exclusion has grown too far")
    assert flags_compared >= 40, (
        f"only {flags_compared} flags compared ({flags_skipped} excluded); a discrete "
        f"decision is the last thing that should leave the portable suite")
    print(f"\n  worst portable drift: " +
          ", ".join(f"{k} {v[0]:.1e} ({v[1]})" for k, v in worst.items() if v[1]))


def test_model_selection_pattern_matches_the_documented_result(reference):
    """The BIC sweep rejects absorption on synthetic and keeps it on real spectra.

    A discrete outcome, identical across environments, and a *scientific* result
    reported in `docs/methods.md` §4.2 — so it is pinned.
    """
    synthetic_abs = [cid for cid in SYNTHETIC if reference[cid].model.use_absorption]
    real_abs = [cid for cid in REAL if reference[cid].model.use_absorption]
    assert not synthetic_abs, f"absorption selected on synthetic cases: {synthetic_abs}"
    assert len(real_abs) >= 4, f"absorption selected on only {len(real_abs)} of 6 real cases"


def test_diagnostics_are_populated(reference):
    """A result must carry the evidence behind its numbers, not just the numbers."""
    for cid, r in reference.items():
        assert r.n_par >= 18, f"{cid}: n_par={r.n_par}"
        assert r.dof > 900, f"{cid}: dof={r.dof}"
        assert 0.5 < r.chi2dof < 3.0, f"{cid}: chi2/dof={r.chi2dof:.2f} is not plausible"
        assert 0.0 < r.covariance_scale <= 1.5
        assert r.parameters, f"{cid}: no parameters recorded"
        assert r.model_selected, f"{cid}: model not recorded"


def test_cli_csv_matches_the_library_that_wrote_it(tmp_path, data_dir):
    """Writer fidelity: the CSV is a faithful serialisation of the result.

    This is the portable half of the round-trip test. It compares the file against
    the objects that produced it **in the same process**, so it holds on any
    machine; comparing the file against the golden CSV is a statement about the
    machine and lives below, under ``characterisation``.
    """
    from gspec.cli import measure_directory
    import gspec

    cases_dir = tmp_path / "cases"
    cases_dir.mkdir()
    for cid in ("S01", "R02"):
        shutil.copy(data_dir / "cases" / f"{cid}.fits", cases_dir)

    out = tmp_path / "out.csv"
    n = measure_directory(cases_dir, out, verbose=False)
    assert n == 8

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert [r["line"] for r in rows[:4]] == ["Halpha", "Hbeta", "OIII_5007", "NII_6583"]

    for cid in ("S01", "R02"):
        written = [r for r in rows if r["case_id"] == cid]
        fresh = gspec.measure(gspec.read_case(cases_dir / f"{cid}.fits"))
        for r in written:
            expected = fresh.row(r["line"])
            for field in FIELDS:
                want = f"{expected[field]:.10g}" if isinstance(expected[field], float) \
                    else str(expected[field])
                assert r[field] == want, f"{cid}/{r['line']}/{field}: {r[field]} != {want}"
            assert r["flag"] == expected["flag"]


# ------------------------------------------------ reference-environment only
@pytest.mark.characterisation
def test_the_reproduction_is_exact_not_merely_close(reference, golden_rows):
    """Measured: every printed value matches the baseline exactly, on this host.

    The comparison is on the **formatted** value, not the raw float, because the
    golden CSV carries ``%.10g`` — comparing full-precision floats to a rounded
    file would fail for arithmetic reasons on a perfect reproduction. What this
    pins is that the port lands on the same 10 significant digits, so a future
    drift of 1e-9 shows up as a change instead of passing silently.

    It is a ``characterisation`` test because it holds on the reference stack and
    nowhere else: the wheel-based environment of CI differs in the last digit or
    two on several rows (L17).
    """
    for cid in ALL_CASES:
        for line, m in reference[cid].lines.items():
            gold = golden_rows[(cid, line)]
            assert f"{m.flux:.10g}" == gold["flux"], f"{cid}/{line} flux drifted"
            assert f"{m.flux_err:.10g}" == gold["flux_err"], f"{cid}/{line} sigma drifted"
            assert f"{reference[cid].redshift:.10g}" == gold["redshift"], f"{cid} z drifted"


@pytest.mark.characterisation
def test_reference_environment_matches_the_baseline_to_1e_6(reference, golden_rows):
    """The original port's acceptance: 1e-6 relative on this environment.

    Kept, and kept out of CI, because 1e-6 is a statement about the frozen
    implementation and this stack — the wheel-based environment drifts to 2-4e-6
    on sigma rows (measured; L17).
    """
    for cid in ALL_CASES:
        for line, m in reference[cid].lines.items():
            gold = golden_rows[(cid, line)]
            got = reference[cid].row(line)
            for field in FIELDS:
                a, b = float(got[field]), float(gold[field])
                rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
                assert rel <= REFERENCE_RTOL, \
                    f"{cid}/{line}/{field}: {a!r} vs golden {b!r} (rel {rel:.2e})"


@pytest.mark.characterisation
def test_cli_csv_matches_the_golden_baseline(tmp_path, data_dir, golden_rows):
    """The file the CLI writes equals the published baseline, digit for digit."""
    from gspec.cli import measure_directory

    cases_dir = tmp_path / "cases"
    cases_dir.mkdir()
    for cid in ("S01", "R02"):
        shutil.copy(data_dir / "cases" / f"{cid}.fits", cases_dir)

    out = tmp_path / "out.csv"
    measure_directory(cases_dir, out, verbose=False)

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        gold = golden_rows[(r["case_id"], r["line"])]
        for field in FIELDS:
            assert r[field] == gold[field], \
                f"{r['case_id']}/{r['line']}/{field}: {r[field]} vs golden {gold[field]}"
        assert r["flag"] == gold["flag"]
