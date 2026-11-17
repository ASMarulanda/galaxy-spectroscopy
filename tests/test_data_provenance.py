"""The claims DATA.md makes about where the data came from, checked against the data.

`DATA.md` §1.1 discloses that ``R06`` is SDSS plate 0294, MJD 51986, fibre 0036,
and that this was *derivable* from two shipped files rather than revealed by the
disclosure. That is a checkable statement about shipped artefacts, so it is checked
here: if either table is ever regenerated, the published reasoning must still hold.

None of this touches the network, and none of it asserts anything about the fitted
numbers — identification is not validation, which is the distinction §1.1 draws.
"""

from __future__ import annotations

import csv
import hashlib

import pytest

from conftest import DATA, ROOT

ARCHIVAL = ROOT / "validation" / "results" / "tables" / "archival_redshift.csv"
PICKS = DATA / "provenance" / "sdss_source_picks.csv"
HASHES = DATA / "provenance" / "sha256_data.txt"

DISCLOSED = (294, 51986, 36)


def _catalogue_redshifts() -> dict[str, float]:
    with open(ARCHIVAL, newline="") as fh:
        return {r["case_id"]: float(r["catalogue_z"]) for r in csv.DictReader(fh)}


def _picks() -> list[dict]:
    with open(PICKS, newline="") as fh:
        return list(csv.DictReader(fh))


def test_the_disclosed_identity_is_the_one_the_shipped_tables_derive():
    """§1.1: matching catalogue_z against picks.zspec identifies R06 uniquely."""
    z_r06 = _catalogue_redshifts()["R06"]
    matches = [p for p in _picks() if abs(float(p["zspec"]) - z_r06) < 1e-6]
    assert len(matches) == 1, f"the match must be unique; got {len(matches)}"
    p = matches[0]
    assert (int(p["plate"]), int(p["mjd"]), int(p["fiberID"])) == DISCLOSED


def test_five_of_six_identities_are_derivable_and_r05_is_not():
    """§1.1 and §6 finding 3: only R05 is genuinely unidentified."""
    picks = _picks()
    derivable, withheld = [], []
    for case_id, z in sorted(_catalogue_redshifts().items()):
        hits = [p for p in picks if abs(float(p["zspec"]) - z) < 1e-6]
        (derivable if len(hits) == 1 else withheld).append(case_id)
    assert derivable == ["R01", "R02", "R03", "R04", "R06"]
    assert withheld == ["R05"]


def test_the_disclosed_archive_hash_is_recorded_in_a_shipped_file():
    """The SHA-256 quoted in §1.1 is checkable offline: it is in sha256_data.txt."""
    name = "spec-0294-51986-0036.fits"
    recorded = {line.split()[1].rsplit("/", 1)[-1]: line.split()[0]
                for line in HASHES.read_text().splitlines() if line.strip()}
    assert recorded[name] == \
        "12e069c2cceb53459b49204f3f508aef3acdd704fd14a136c0b6449e8a118554"


def test_the_case_file_is_the_shipped_checklisted_one():
    """The published `data/checksums.sha256` covers every file under `data/`."""
    listed = {}
    for line in (DATA / "checksums.sha256").read_text().splitlines():
        if not line.strip():
            continue
        digest, path = line.split(maxsplit=1)
        listed[path.strip().lstrip("*").removeprefix("./")] = digest
    key = next(k for k in listed if k.endswith("cases/R06.fits"))
    assert hashlib.sha256((ROOT / key).read_bytes()).hexdigest() == listed[key]


def test_the_disclosure_is_stated_in_the_data_record():
    """A claim in a test is worthless if the document stops making it."""
    text = (ROOT / "DATA.md").read_text()
    assert "### 1.1 The R06 identity, disclosed" in text
    assert "plate 0294, MJD 51986, fibre 0036" in text
    # and it must keep drawing the distinction the test file is named after
    assert "**What it does not establish.**" in text
