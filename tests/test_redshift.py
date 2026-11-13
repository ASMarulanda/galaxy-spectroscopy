"""C1, the redshift half, and the internal consistency of the reported value.

Two separate claims:

* the recovered redshift matches analytic truth to better than 150 km/s — 2.5x
  the instrumental resolution and about 7x the worst error measured;
* one spectrum has one systemic velocity, so its four rows repeat one redshift.
  That is a *definition* the pipeline guarantees by construction, and the test
  exists to catch a future refactor that fits a redshift per line.
"""

from __future__ import annotations

import statistics as st

import pytest

from conftest import INJECTED, REAL, SYNTHETIC

C_KMS = 299792.458
TOL_KMS = 150.0


def dz_kms(z, true_z):
    return (z - true_z) * C_KMS / (1.0 + true_z)


def test_synthetic_redshifts_are_within_150_km_s(reference, synthetic_truth):
    bad = []
    for cid in SYNTHETIC:
        dv = dz_kms(reference[cid].redshift, synthetic_truth[cid]["z"])
        if abs(dv) > TOL_KMS:
            bad.append(f"{cid}: dv={dv:+.1f} km/s")
    assert not bad, "redshift outside 150 km/s:\n  " + "\n  ".join(bad)


@pytest.mark.characterisation
def test_measured_redshift_accuracy(reference, synthetic_truth):
    """Median 1.7 km/s, max 18.5 km/s on this sample (validated: 1.8 / 20.0)."""
    dvs = [abs(dz_kms(reference[cid].redshift, synthetic_truth[cid]["z"])) for cid in SYNTHETIC]
    assert st.median(dvs) < 5.0
    assert max(dvs) < 40.0


def test_every_row_of_a_case_repeats_one_redshift(reference):
    for cid, r in reference.items():
        zs = {r.row(line)["redshift"] for line in r.lines}
        assert len(zs) == 1, f"{cid} reports {len(zs)} different redshifts"


def test_redshifts_are_physically_plausible(reference):
    for cid, r in reference.items():
        assert 0.008 < r.redshift < 0.102, f"{cid}: z={r.redshift} outside the search range"
        assert 0 < r.redshift_err < 0.05, f"{cid}: sigma_z={r.redshift_err}"


def test_injected_copies_share_their_base_redshift(reference):
    """An injection adds a line at the object's own redshift: z must not move much."""
    for iid in INJECTED:
        base = "R" + iid[1:]
        dv = dz_kms(reference[iid].redshift, reference[base].redshift)
        assert abs(dv) < TOL_KMS, f"{iid} vs {base}: redshift moved by {dv:+.1f} km/s"


@pytest.mark.parametrize("cid", REAL)
def test_real_redshifts_land_in_the_expected_neighbourhood(reference, cid):
    """No exact truth is published for the real cases; this is a sanity range.

    The archive redshifts are 0.0206-0.0812 (DATA.md). The wide band is
    deliberate: it catches a gross misidentification (the ~950 km/s doublet
    alias), not a small error.
    """
    z = reference[cid].redshift
    assert 0.015 < z < 0.095, f"{cid}: z={z} is outside the archive range"
