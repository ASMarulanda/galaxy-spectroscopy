"""Produce the golden reference CSV from the FROZEN original implementation.

This script is deliberately *not* part of the gspec package and never will be:
its whole purpose is to run the pre-migration code so the baseline is
independent of the thing it is used to check. It reads the frozen tree
read-only and writes only under validation/.

Usage:
    python3 /tmp/kp0/make_golden.py <repo> <frozen-solution-dir> <cases-dir>
"""
import csv
import glob
import hashlib
import json
import math
import os
import platform
import statistics as st
import sys
import time

import numpy as np
import scipy
import astropy

REPO = sys.argv[1]
FROZEN = sys.argv[2]
CASES = sys.argv[3]

OUT_CSV = os.path.join(REPO, "validation/results/golden_reference.csv")
OUT_META = os.path.join(REPO, "validation/results/golden_reference.meta.json")
TRUTH = os.path.join(REPO, "data/truth/truth_synthetic.json")

sys.path.insert(0, FROZEN)
import common as C          # noqa: E402
from methods import fit_spectrum   # noqa: E402

LINES = ("Halpha", "Hbeta", "OIII_5007", "NII_6583")
COLS = ["case_id", "line", "redshift", "redshift_err", "flux", "flux_err", "flag"]

files = sorted(glob.glob(os.path.join(CASES, "*.fits")))
t0 = time.time()
rows = []
per_case = {}
for fn in files:
    case = C.read_case(fn)
    cid = case["case_id"] or os.path.basename(fn)[:-5]
    ts = time.time()
    r = fit_spectrum(case)
    per_case[cid] = time.time() - ts
    for line in LINES:
        rec = r.get("lines", {}).get(line)
        if rec is None:
            rows.append(dict(case_id=cid, line=line, redshift=float("nan"),
                             redshift_err=float("nan"), flux=float("nan"),
                             flux_err=float("nan"), flag="failed"))
        else:
            rows.append(dict(case_id=cid, line=line, redshift=r["z"],
                             redshift_err=r["z_err"], flux=rec["flux"],
                             flux_err=rec["flux_err"], flag=rec["flag"]))
wall = time.time() - t0

with open(OUT_CSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    for row in rows:
        w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in row.items()})

sha_csv = hashlib.sha256(open(OUT_CSV, "rb").read()).hexdigest()

# ---------------------------------------------------------------- summary
syn = json.load(open(TRUTH))
FL = {"Halpha": 0.07, "Hbeta": 0.20, "OIII_5007": 0.06, "NII_6583": 0.10}
acc = {}
for line in LINES:
    errs, within = [], 0
    for cid in sorted(syn):
        t = syn[cid]["lines"][line]["flux"]
        row = next(r for r in rows if r["case_id"] == cid and r["line"] == line)
        if row["flag"] == "failed" or t == 0:
            continue
        e = (row["flux"] - t) / t
        errs.append(100 * e)
        within += abs(e) <= FL[line]
    acc[line] = dict(n=len(errs), median_pct=round(st.median(errs), 4),
                     worst_pct=round(max(errs, key=abs), 4),
                     within_floor=int(within), floor_pct=100 * FL[line])

dz = []
for cid in sorted(syn):
    tz = syn[cid]["z"]
    z = next(r for r in rows if r["case_id"] == cid)["redshift"]
    dz.append((z - tz) * 299792.458 / (1 + tz))

flags = {}
for r in rows:
    flags[r["flag"]] = flags.get(r["flag"], 0) + 1

meta = dict(
    produced="golden reference produced by the FROZEN pre-migration implementation",
    frozen_solution_dir=FROZEN,
    frozen_methods_sha256=hashlib.sha256(
        open(os.path.join(FROZEN, "methods.py"), "rb").read()).hexdigest(),
    frozen_common_sha256=hashlib.sha256(
        open(os.path.join(FROZEN, "common.py"), "rb").read()).hexdigest(),
    cases_dir=CASES,
    csv_sha256=sha_csv,
    csv_rows=len(rows),
    environment=dict(python=platform.python_version(), numpy=np.__version__,
                     scipy=scipy.__version__, astropy=astropy.__version__,
                     platform=platform.platform(), machine=platform.machine()),
    wall_seconds=round(wall, 2),
    per_case_seconds={k: round(v, 2) for k, v in sorted(per_case.items())},
    median_case_seconds=round(st.median(per_case.values()), 2),
    synthetic_accuracy_vs_truth=acc,
    redshift_kms=dict(median_abs=round(st.median([abs(x) for x in dz]), 2),
                      max_abs=round(max(abs(x) for x in dz), 2), n=len(dz)),
    flag_counts=flags,
)
json.dump(meta, open(OUT_META, "w"), indent=1, sort_keys=True)

print(f"wrote {OUT_CSV}")
print(f"  sha256 {sha_csv}")
print(f"  {len(rows)} rows, {wall:.1f}s, median {st.median(per_case.values()):.2f}s/case")
for line, a in acc.items():
    print(f"  {line:10s} median {a['median_pct']:+8.2f}%  worst {a['worst_pct']:+8.2f}%  "
          f"within floor {a['within_floor']}/{a['n']}")
print(f"  flags {flags}")
