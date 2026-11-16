#!/usr/bin/env python3
"""Regenerate every figure and record what produced it.

    .venv/bin/python validation/figures/make_all.py      # ~6 min on 2 CPUs

Runs the five figure scripts in order, then writes
``validation/figures/FIGURE_MANIFEST.json`` containing:

* the software versions the figures were drawn with,
* the SHA-256 of every input file the figures consume,
* the SHA-256 of every output, so that a re-run can be compared byte for byte,
* the command that produced each figure.

If a run reproduces the committed hashes, the figures are reproducible in the
narrow sense that matters: same inputs, same code, same bytes.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

SCRIPTS = [
    ("fig1_s01_fit.py", "Figure 1 — S01 fit and residuals"),
    ("fig2_hbeta_configurations.py", "Figure 2 — Hβ model dependence"),
    ("fig3_pull_calibration.py", "Figure 3 — uncertainty calibration"),
    ("fig4_injection_recovery.py", "Figure 4 — injection recovery"),
    ("fig5_diagnostics.py", "Figure 5 — line-ratio diagnostics"),
]

INPUTS = [
    "data/cases/*.fits",
    "data/truth/truth_synthetic.json",
    "data/truth/truth_injections.json",
    "validation/results/golden_reference.csv",
    "validation/results/out_independent.csv",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def versions() -> dict:
    import matplotlib
    import numpy
    import scipy
    import astropy
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "astropy": astropy.__version__,
        "matplotlib": matplotlib.__version__,
    }


def main() -> int:
    for script, label in SCRIPTS:
        print(f"\n=== {label} ===")
        rc = subprocess.call([sys.executable, str(HERE / script)], cwd=ROOT)
        if rc != 0:
            print(f"FAILED: {script} (exit {rc})")
            return rc

    outputs = sorted([q for ext in ("png", "pdf") for q in HERE.glob(f"fig*.{ext}")])
    inputs = []
    for pattern in INPUTS:
        inputs.extend(sorted(ROOT.glob(pattern)))

    manifest = {
        "generated_by": "validation/figures/make_all.py",
        "versions": versions(),
        "scripts": [dict(script=s, purpose=p) for s, p in SCRIPTS],
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "outputs": {p.name: sha256(p) for p in outputs},
    }
    out = HERE / "FIGURE_MANIFEST.json"
    out.write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}  "
          f"({len(manifest['inputs'])} inputs, {len(manifest['outputs'])} outputs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
