"""E-Cov: what happens to sigma_F when a parameter sits on its bound?

The question
------------
``cov = pinv(JtJ) * min(chi2/dof, 1.5)`` is the shipped estimator and is **not**
changed by this experiment. The question is narrower: when the BIC sweep pins an
amplitude at zero (the broad-Hbeta component, in most cases), the estimator keeps
a direction in the covariance whose curvature is essentially undefined. How much
of ``sigma_F`` comes from such a direction, and what happens if those directions
are dropped instead?

The two estimators
------------------
* **shipped** — pinv over all parameters. This is what the package does and what
  every number in validation/results/ was produced with.
* **conditioned** — the bound-pinned parameters are removed from the problem
  before the covariance is formed, i.e. the uncertainty is computed *given* that
  the model holds them fixed. Mathematically that is conditioning on the boundary,
  which is what "this component is zero" asserts.

Measuring the share
-------------------
For each row the variance is decomposed over the eigen-directions of ``JtJ``:

    sigma_F^2 = sum_k (v_k . g)^2 / lambda_k

The direction with the largest share is reported with its eigenvalue and its
position relative to the pseudoinverse's truncation cutoff
(``rcond * lambda_max``, rcond = 1e-15). A direction close to that cutoff is one
whose curvature the data barely constrain.

Run:  python validation/experiments/e_cov_bound_pinned.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

import gspec
from gspec import constants as K
from gspec.covariance import covariance, covariance_scale, gradient
from gspec.fit import build_complexes
from gspec.flags import classify
from gspec.fluxes import line_total
from gspec.pipeline import measure, measure_fixed_model

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = Path(__file__).resolve().parent / "results"

CASES = tuple(sorted(f"S{i:02d}" for i in range(1, 11))
               + sorted(f"R{i:02d}" for i in range(1, 7))
               + sorted(f"I{i:02d}" for i in range(1, 5)))
RCOND = 1e-15


def sigma_with(cov: np.ndarray, g: np.ndarray, keep: np.ndarray | None = None) -> float:
    if keep is None:
        return math.sqrt(max(float(g @ cov @ g), 0.0))
    gg = g[keep]
    return math.sqrt(max(float(gg @ cov @ gg), 0.0))


def main() -> int:
    rows = []
    for cid in CASES:
        case = gspec.read_case(DATA / "cases" / f"{cid}.fits")
        shipped = measure(case)
        # Refit the *selected* model at the *reported* redshift to obtain J.
        refit = measure_fixed_model(case, shipped.model, z_known=shipped.redshift,
                                    keep_jacobian=True)
        if refit.jacobian is None:
            print(f"  {cid}: no Jacobian, skipped")
            continue

        names = list(shipped.parameters)
        p = np.array([shipped.parameters[n] for n in names])
        idx = {n: i for i, n in enumerate(names)}
        jac = refit.jacobian
        scale = covariance_scale(refit.chi2, refit.dof)
        cov = covariance(jac, refit.chi2, refit.dof)

        pinned = set(refit.at_bounds)
        keep = np.array([n not in pinned for n in names])
        jr = jac[:, keep]
        cov_cond = np.linalg.pinv(jr.T @ jr) * scale

        lam, vec = np.linalg.eigh(jac.T @ jac)
        cutoff = RCOND * lam.max()
        complexes = build_complexes(refit.model)

        for line in K.LINE_NAMES:
            fn = (lambda pv, _l=line: line_total(pv, idx, case, complexes, _l))
            f0 = fn(p)
            g = gradient(fn, p)
            s_ship = sigma_with(cov, g)
            s_cond = sigma_with(cov_cond, g, keep)

            proj2 = (vec.T @ g) ** 2
            contrib = proj2 / np.maximum(lam, 1e-300) * scale
            kept = lam > cutoff
            total = contrib[kept].sum()
            if total > 0:
                k = int(np.argmax(np.where(kept, contrib, 0.0)))
                share = float(contrib[k] / total)
                dom = names[int(np.argmax(np.abs(vec[:, k])))]
                lam_k = float(lam[k])
            else:
                share, dom, lam_k = float("nan"), "", float("nan")

            rows.append(dict(
                case=cid, line=line,
                n_par=refit.n_par, pinned=",".join(sorted(pinned)),
                flux=round(f0, 6),
                sigma_shipped=round(s_ship, 6), sigma_conditioned=round(s_cond, 6),
                ratio=round(s_ship / s_cond, 4) if s_cond else None,
                flag_shipped=shipped.lines[line].flag,
                flag_conditioned=classify(f0, s_cond,
                                          shipped.lines[line].broad_amplitude),
                dominant_param=dom, dominant_share_pct=round(100 * share, 2),
                lambda_dominant=lam_k,
                lambda_max=float(lam.max()),
                lambda_over_cutoff=round(lam_k / cutoff, 3) if cutoff else None,
            ))
        print(f"  {cid}: pinned={sorted(pinned)}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "e_cov_bound_pinned.json").write_text(json.dumps(rows, indent=1))

    # ------------------------------------------------------------- report
    print("\n=== sigma_F: shipped vs conditioned on the bound (all rows) ===")
    ratios = [r["ratio"] for r in rows if r["ratio"]]
    print(f"  rows: {len(rows)}   ratio shipped/conditioned: "
          f"min {min(ratios):.3f}  median {sorted(ratios)[len(ratios)//2]:.3f}  max {max(ratios):.3f}")

    print("\n  rows where the two estimators differ by more than 10 %:")
    print(f"    {'case':5} {'line':11} {'shipped':>10} {'conditioned':>12} {'ratio':>7} "
          f"{'dominant param':16} {'share':>7} {'lam/cutoff':>11}")
    movers = [r for r in rows if r["ratio"] and abs(r["ratio"] - 1) > 0.10]
    for r in sorted(movers, key=lambda x: -abs(x["ratio"] - 1)):
        print(f"    {r['case']:5} {r['line']:11} {r['sigma_shipped']:10.4f} "
              f"{r['sigma_conditioned']:12.4f} {r['ratio']:7.3f} {r['dominant_param']:16} "
              f"{r['dominant_share_pct']:6.1f}% {r['lambda_over_cutoff']:11.2e}")
    print(f"  ({len(movers)} of {len(rows)} rows move by more than 10 %)")

    flips = [r for r in rows if r["flag_shipped"] != r["flag_conditioned"]]
    print(f"\n  flag changes if the covariance were conditioned: {len(flips)}")
    for r in flips:
        print(f"    {r['case']} {r['line']}: {r['flag_shipped']} -> {r['flag_conditioned']}")

    print("\n=== how much of sigma_F^2 comes from a weakly constrained direction? ===")
    print(f"  {'case':5} {'line':11} {'share of sigma^2':>16} {'lambda/lambda_max':>18} "
          f"{'lambda/cutoff':>14}  dominant")
    for cid in ("S01", "S04", "S08"):
        for r in rows:
            if r["case"] == cid and r["line"] in ("Halpha", "Hbeta"):
                print(f"  {r['case']:5} {r['line']:11} {r['dominant_share_pct']:15.1f}% "
                      f"{r['lambda_dominant']/r['lambda_max']:18.2e} "
                      f"{r['lambda_over_cutoff']:14.2e}  {r['dominant_param']}")

    print(f"\nwrote {OUT / 'e_cov_bound_pinned.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
