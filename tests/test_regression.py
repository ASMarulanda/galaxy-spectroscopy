"""Does the package still produce the validated baseline?

The golden file was produced by the *frozen pre-migration implementation* on this
same environment, so the comparison is tier 2 of the reproduction protocol: the
ported code must agree with the frozen code to floating-point precision, not
merely within the acceptance tolerances.

Why not byte equality of the files: the CSV is written with ``%.10g``, so byte
equality would test the formatter as much as the physics — and the pipeline is
explicitly *not* bit-reproducible across environments (one sigma row moves by
2.6x between a container and this host). What is asserted here is that the
reference run reproduces every value to 1e-6 relative **and** that no flag
changed, since a flag is a discrete decision that no tolerance can hide.
"""

from __future__ import annotations

import csv
import shutil

import pytest

from conftest import ALL_CASES, REAL, SYNTHETIC

FIELDS = ("redshift", "redshift_err", "flux", "flux_err")
RTOL = 1e-6


def test_every_value_matches_the_golden_baseline(reference, golden_rows):
    assert len(golden_rows) == 80
    worst = {}
    for line in FIELDS:
        worst[line] = 0.0
    for cid in ALL_CASES:
        for line, m in reference[cid].lines.items():
            gold = golden_rows[(cid, line)]
            got = reference[cid].row(line)
            for field in FIELDS:
                a, b = float(got[field]), float(gold[field])
                rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
                worst[field] = max(worst[field], rel)
                assert rel <= RTOL, f"{cid}/{line}/{field}: {a!r} vs golden {b!r} (rel {rel:.2e})"
            assert got["flag"] == gold["flag"], f"{cid}/{line}: flag {got['flag']} vs golden {gold['flag']}"


@pytest.mark.characterisation
def test_the_reproduction_is_exact_not_merely_close(reference, golden_rows):
    """Measured: every printed value matches the baseline exactly, on this host.

    The comparison is on the **formatted** value, not the raw float, because the
    golden CSV carries ``%.10g`` — comparing full-precision floats to a rounded
    file would fail for arithmetic reasons on a perfect reproduction. What this
    pins is that the port lands on the same 10 significant digits, so a future
    drift of 1e-9 shows up as a change instead of passing silently.
    """
    for cid in ALL_CASES:
        for line, m in reference[cid].lines.items():
            gold = golden_rows[(cid, line)]
            assert f"{m.flux:.10g}" == gold["flux"], f"{cid}/{line} flux drifted"
            assert f"{m.flux_err:.10g}" == gold["flux_err"], f"{cid}/{line} sigma drifted"
            assert f"{reference[cid].redshift:.10g}" == gold["redshift"], f"{cid} z drifted"


def test_csv_export_round_trips_through_the_cli(tmp_path, data_dir, golden_rows):
    """The optional CLI writes the same numbers the library computed."""
    from gspec.cli import measure_directory

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
    for r in rows:
        gold = golden_rows[(r["case_id"], r["line"])]
        for field in FIELDS:
            assert float(r[field]) == float(gold[field]), \
                f"{r['case_id']}/{r['line']}/{field} differs from the baseline"
        assert r["flag"] == gold["flag"]


def test_model_selection_pattern_matches_the_documented_result(reference):
    """The BIC sweep rejects absorption on synthetic and keeps it on real spectra.

    This is a *scientific* result reported in `docs/methods.md` §4.2, so it is
    pinned: if it changes, the documentation is wrong.
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
