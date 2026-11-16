#!/usr/bin/env python3
"""Run the independent analysis exactly as it produced results/out_independent.csv.

This is a thin driver: the analysis itself is measure_independent2.py, copied
unmodified from the authoring workspace.  Nothing is re-implemented here.

    python3 run_independent.py                       # defaults, from this directory
    python3 run_independent.py --cases ../../data/cases \
                              --out ../results/out_independent_reproduced.csv

The default --cases path is ../../data/cases, which resolves to the shipped case
files.  The original run took about 45 minutes for the 20 cases on 2 CPUs; about
25 minutes on a modern desktop.
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import measure_independent2 as MI  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default=os.path.join(HERE, "..", "..", "data", "cases"),
                    help="directory of .fits case files")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "out_independent_reproduced.csv"))
    args = ap.parse_args()

    cases = os.path.abspath(args.cases)
    out = os.path.abspath(args.out)
    if not os.path.isdir(cases):
        raise SystemExit(f"no case directory at {cases}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    print(f"cases: {cases}\nout:   {out}", flush=True)
    MI.main(cases, out)


if __name__ == "__main__":
    main()
