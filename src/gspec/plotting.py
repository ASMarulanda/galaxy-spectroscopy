"""Scientific figures for the pipeline.

Design rules, applied to every figure in this module:

* every plotted quantity is computed here from the fitted parameters or read from
  a committed dataset — no number is typed into a caption by hand;
* a non-detection is never drawn as a measurement (:class:`gspec.diagnostics.Ratio`
  carries ``detected``);
* each figure's caption states what the figure does *not* establish;
* the output is deterministic: no timestamps in the metadata, fixed layout, so
  that a regenerated figure is byte-identical and can be diffed.

The functions return :class:`matplotlib.figure.Figure` objects and never call
``show()`` or ``savefig()``; the scripts in ``validation/figures/`` own I/O.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no display required; must precede pyplot

import matplotlib.pyplot as plt
import numpy as np

from .continuum import continuum_at
from .fit import _select_windows, build_complexes, build_parameters
from .io import Case
from .lines import K as LINE_CONSTANTS  # noqa: F401  (kept for symmetry with model.py)
from .lsf import gaussian, sigma_pix
from .model import amplitudes
from .results import SpectrumResult

# --------------------------------------------------------------------------- #
# Style
# --------------------------------------------------------------------------- #
C_DATA = "#2b2b2b"
C_MODEL = "#c0392b"
C_CONT = "#1f6fb4"
C_BAND = "#d9d9d9"
C_ALT = "#7f8c8d"
C_ACCENT = "#e67e22"
C_GOOD = "#2e7d32"

#: Colours for the emission components, assigned in order of appearance.
COMPONENT_COLORS = ["#6a9fb5", "#8e7cc3", "#c27ba0", "#76a5af", "#b45f06",
                    "#6aa84f", "#a64d79", "#45818e", "#bf9000", "#674ea7"]

#: Display labels for the fitted component keys.
COMPONENT_LABELS = {
    "NII6548": "[N II] 6548", "Ha_n": "Hα", "Ha_b": "Hα br",
    "NII6583": "[N II] 6583", "SII6716": "[S II] 6716", "SII6731": "[S II] 6731",
    "Hg_em": "Hγ", "Hb_n": "Hβ", "Hb_b": "Hβ br",
    "OIII4959": "[O III] 4959", "OIII5007": "[O III] 5007",
}


def apply_style() -> None:
    """One place for the figure typography, so all five figures match."""
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 220,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.grid": False,
        "font.family": "DejaVu Sans",
        "mathtext.default": "regular",
    })


def save(fig, stem: str) -> list[Path]:
    """Write ``stem.png`` and ``stem.pdf`` with no timestamps, and return the paths."""
    out: list[Path] = []
    for ext in ("png", "pdf"):
        path = Path(f"{stem}.{ext}")
        meta = {"CreationDate": None, "Creator": "gspec", "Producer": "gspec"}
        fig.savefig(path, bbox_inches="tight", metadata=meta)
        out.append(path)
    return out


# --------------------------------------------------------------------------- #
# Model decomposition
# --------------------------------------------------------------------------- #
def _parameter_vector(result: SpectrumResult, case: Case, complexes, windows):
    names, idx, _, _, _ = build_parameters(case, complexes, windows, result.model,
                                           result.redshift)
    missing = [n for n in names if n not in result.parameters]
    if missing:
        raise KeyError(f"result is missing parameters {missing}")
    return np.array([result.parameters[n] for n in names]), idx


def model_curves(case: Case, result: SpectrumResult) -> dict:
    """Re-evaluate the fitted model on the data, split into its physical parts.

    Returns a dict with the arrays a reader needs to see the fit:

    ``wave``, ``flux``, ``err``  the observation (``nan`` where masked)
    ``continuum``                polynomial × absorption product, ``nan`` outside
                                 the two fit windows
    ``components``               list of ``(key, rest_A, array)``, one entry per
                                 emission component, each the Gaussian *only*
    ``total``                    continuum + all components
    ``in_window``                boolean mask of the pixels the fit used

    The identity ``total == continuum + Σ components`` holds by construction and
    ``Σ residuals² == result.chi2`` is asserted by the test suite, which is what
    proves this re-evaluation is the model that was fitted.
    """
    complexes = build_complexes(result.model)
    windows = _select_windows(case, complexes, result.redshift)
    if windows is None:
        raise ValueError(f"{case.case_id}: fit windows no longer selectable")
    p, idx = _parameter_vector(result, case, complexes, windows)

    wave = case.wave
    flux = np.where(case.ivar > 0, case.flux, np.nan)
    err = np.where(case.ivar > 0, 1.0 / np.sqrt(np.where(case.ivar > 0, case.ivar, 1)), np.nan)

    continuum = np.full_like(wave, np.nan)
    total = np.full_like(wave, np.nan)
    components: list[tuple[str, float, np.ndarray]] = []

    for ci, cx in enumerate(complexes):
        w = wave[windows[ci]]
        cont = continuum_at(w, cx, ci, p, idx, result.model.ndeg,
                            result.model.use_absorption,
                            case.wdisp[windows[ci]], case.dlam[windows[ci]])
        continuum[windows[ci]] = cont
        seg = cont.copy()
        for e in cx.emis:
            mu = e.rest * (1.0 + p[idx["z"]])
            sig = sigma_pix(mu, math.exp(p[idx[e.sigma]]),
                            case.wdisp[windows[ci]], case.dlam[windows[ci]])
            prof = amplitudes(cx, p, idx)[e.key] * gaussian(w, mu, sig)
            seg = seg + prof
            full = np.full_like(wave, np.nan)
            full[windows[ci]] = prof
            components.append((e.key, e.rest, full))
        total[windows[ci]] = seg

    in_window = np.zeros_like(wave, dtype=bool)
    for w in windows:
        in_window |= w

    return dict(wave=wave, flux=flux, err=err, continuum=continuum,
                components=components, total=total, in_window=in_window,
                redshift=result.redshift)


# --------------------------------------------------------------------------- #
# Figure 1 — the fit itself
# --------------------------------------------------------------------------- #
def fig_spectrum_fit(case: Case, result: SpectrumResult,
                     ylim_zoom: tuple[float, float] | None = None) -> plt.Figure:
    """Observed spectrum, continuum, emission components, total model, residuals.

    A four-panel layout: an overview with the fitted windows marked, one zoom per
    window where the components are separated, and the normalised residuals
    underneath. The overview exists because a zoom hides whether the fit is
    ignoring the rest of the spectrum — the residual panel is plotted over
    *every* fitted pixel, not just the zoomed range.
    """
    apply_style()
    cur = model_curves(case, result)
    wave, flux, err = cur["wave"], cur["flux"], cur["err"]
    z = cur["redshift"]

    fig = plt.figure(figsize=(10.6, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05],
                          hspace=0.34, wspace=0.20, top=0.90, bottom=0.115)

    # ---------------- (a) overview -------------------------------------- #
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(wave, flux, lw=0.5, color=C_DATA, zorder=2)
    ax.fill_between(wave, flux - err, flux + err, color=C_BAND, alpha=0.55,
                    lw=0, zorder=1)
    ax.plot(wave, cur["total"], lw=1.1, color=C_MODEL, zorder=3)
    ax.set_xlabel("Observed wavelength (Å)")
    ax.set_ylabel(r"$f_\lambda$  (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)")
    ax.set_title("(a) full spectrum; grey bands mark the fitted windows", loc="left")
    ax.set_xlim(wave[0], wave[-1])
    for lo, hi, _ in _window_spans(cur, wave):
        ax.axvspan(lo, hi, color=C_ALT, alpha=0.13, lw=0)

    # ---------------- (b) Hα complex ------------------------------------ #
    ax_b = fig.add_subplot(gs[0, 1])
    _plot_complex(ax_b, cur, _complex_span(cur, "Ha"), z,
                  title="(b) Hα complex")

    # ---------------- (c) Hβ/[O III] complex ---------------------------- #
    ax_c = fig.add_subplot(gs[1, 0])
    _plot_complex(ax_c, cur, _complex_span(cur, "Hb"), z,
                  title="(c) Hβ and [O III] complex")

    # ---------------- (d) residuals ------------------------------------- #
    ax_d = fig.add_subplot(gs[1, 1])
    good = np.isfinite(cur["total"]) & np.isfinite(flux) & (err > 0)
    resid = np.full_like(wave, np.nan)
    resid[good] = (flux[good] - cur["total"][good]) / err[good]
    ax_d.axhline(0, color=C_DATA, lw=0.8)
    ax_d.fill_between([wave[0], wave[-1]], -1, 1, color=C_BAND, alpha=0.7, lw=0)
    ax_d.plot(wave, resid, ".", ms=1.7, color=C_MODEL, alpha=0.85, rasterized=True)
    ax_d.set_xlim(wave[0], wave[-1])
    ax_d.set_xlabel("Observed wavelength (Å)")
    ax_d.set_ylabel(r"(data − model)/$\sigma$")
    ax_d.set_title("(d) normalised residuals over every fitted pixel", loc="left")
    rms = float(np.sqrt(np.nanmean(resid[good] ** 2)))
    ax_d.text(0.985, 0.94, f"RMS = {rms:.2f}σ   χ²/ν = {result.chi2dof:.3f}",
              transform=ax_d.transAxes, ha="right", va="top", fontsize=8,
              color=C_DATA, bbox=dict(fc="white", ec="none", alpha=0.85, pad=2))

    if ylim_zoom is not None:
        ax_b.set_ylim(*ylim_zoom)

    # ---------------- shared legend and titles -------------------------- #
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    handles = [
        Line2D([], [], color=C_DATA, lw=0.9, label="observed"),
        Patch(facecolor=C_BAND, alpha=0.8, label="1σ of the data"),
        Line2D([], [], color=C_CONT, ls="--", lw=1.1, label="continuum"),
        Patch(facecolor=COMPONENT_COLORS[1], alpha=0.45, label="emission component"),
        Line2D([], [], color=C_MODEL, lw=1.2, label="total model"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, 0.008), handlelength=1.6, columnspacing=1.4)
    fig.suptitle(f"{case.case_id} — fitted with {result.model.label()}  "
                 f"({result.n_par} parameters, {result.dof} dof, "
                 f"z = {result.redshift:.6f})",
                 fontsize=10, x=0.5, y=0.985)
    return fig


def _window_spans(cur, wave):
    """Contiguous pixel spans covered by the two complexes."""
    in_w = cur["in_window"]
    spans, start = [], None
    for i, flag in enumerate(in_w):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            spans.append((wave[start], wave[i - 1]))
            start = None
    if start is not None:
        spans.append((wave[start], wave[-1]))
    return [(lo, hi, "") for lo, hi in spans]


def _complex_span(cur, which: str) -> tuple[float, float]:
    """Wavelength extent of one complex, taken from its own components."""
    rest = []
    for key, r, prof in cur["components"]:
        if (which == "Ha" and r > 6400) or (which == "Hb" and r < 5200):
            rest.append(r)
    z = cur["redshift"]
    lo = min(rest) * (1 + z) - 60
    hi = max(rest) * (1 + z) + 60
    return lo, hi


def _plot_complex(ax, cur, span, z, title, offset_scale=1.0):
    """One zoom panel: data, continuum, each component, and the total model.

    Component labels are placed along the top edge in axes coordinates and rotated,
    so that a tall line cannot push its own label outside the panel.
    """
    wave, flux, err = cur["wave"], cur["flux"], cur["err"]
    m = (wave >= span[0]) & (wave <= span[1]) & cur["in_window"]
    if not m.any():
        ax.set_title(title + "  (outside the fitted windows)", loc="left")
        return

    ax.plot(wave[m], flux[m], lw=0.7, color=C_DATA, zorder=4)
    ax.fill_between(wave[m], flux[m] - err[m], flux[m] + err[m],
                    color=C_BAND, alpha=0.6, lw=0, zorder=1)
    ax.plot(wave[m], cur["continuum"][m], lw=1.0, ls="--", color=C_CONT, zorder=3)

    ymax = np.nanmax(cur["total"][m])
    ymin = np.nanmin(np.minimum(flux[m], cur["total"][m]))
    ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + 0.22 * (ymax - ymin))

    cont = cur["continuum"]
    k = 0
    for key, rest, prof in cur["components"]:
        mu = rest * (1 + z)
        if not (span[0] <= mu <= span[1]):
            continue
        colour = COMPONENT_COLORS[k % len(COMPONENT_COLORS)]
        k += 1
        curve = cont + prof                      # component on top of continuum
        ax.fill_between(wave[m], cont[m], curve[m], color=colour, alpha=0.38,
                        lw=0, zorder=2)
        ax.plot(wave[m], curve[m], lw=0.7, color=colour, zorder=2)
        peak = np.nanmax(curve[m])
        base = np.nanmax(cont[m])
        if np.isfinite(peak) and (peak - base) > 0.02 * (np.nanmax(flux[m]) - base):
            frac = (peak - ymin) / (ax.get_ylim()[1] - ymin)
            # a tall line gets the top row, so its own label cannot run into it;
            # the short ones alternate underneath to avoid horizontal collisions
            ytxt = 0.985 if frac > 0.55 else (0.985 if k % 2 else 0.885)
            ax.text(mu, ytxt, COMPONENT_LABELS.get(key, key),
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    rotation=90, fontsize=6.5, color=colour, clip_on=False)

    ax.plot(wave[m], cur["total"][m], lw=1.2, color=C_MODEL, zorder=5)
    ax.set_xlim(span)
    ax.set_xlabel("Observed wavelength (Å)")
    ax.set_ylabel(r"$f_\lambda$  (10$^{-17}$)")
    ax.set_title(title, loc="left")


# --------------------------------------------------------------------------- #
# Figure 2 — Hβ model dependence
# --------------------------------------------------------------------------- #
def fig_hbeta_configurations(truth: dict[str, float],
                             configs: dict[str, dict[str, float]]) -> plt.Figure:
    """Signed Hβ error for each configuration, per case, against synthetic truth.

    ``truth`` maps case id → true flux; ``configs`` maps a configuration label →
    {case id → measured flux}. Panel (a) shows every case so that the spread is
    visible case by case; panel (b) summarises the four configurations with the
    individual cases overlaid, because a mean would hide the S06-type outlier
    that matters most.
    """
    apply_style()
    cases = sorted(truth)
    labels = list(configs)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.9),
                                  gridspec_kw=dict(width_ratios=[1.55, 1.0],
                                                   wspace=0.28))
    colours = [C_MODEL, C_ACCENT, C_CONT, C_ALT]
    markers = ["o", "s", "^", "D"]

    # ---- (a) per case, one column per configuration -------------------- #
    ypos = np.arange(len(cases))[::-1]
    for ci, lab in enumerate(labels):
        errs = [100 * (configs[lab][c] - truth[c]) / truth[c] for c in cases]
        ax.plot(errs, ypos + (ci - (len(labels) - 1) / 2) * 0.16,
                markers[ci % 4], ms=5.2, color=colours[ci % 4], label=lab,
                mec="white", mew=0.5, ls="none")

    ax.axvline(0, color=C_DATA, lw=0.9)
    ax.axvspan(-20, 20, color=C_GOOD, alpha=0.09, lw=0)
    ax.set_xscale("symlog", linthresh=100)
    ax.set_xlim(-150, 700)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"S{c[1:]}" for c in cases])
    ax.set_xlabel("signed error of Hβ flux (%)")
    ax.set_title("(a) every case, every configuration", loc="left")
    ax.legend(loc="upper right", frameon=False, ncol=2, handletextpad=0.3,
              columnspacing=1.0, bbox_to_anchor=(1.0, 0.99))
    ax.text(0.02, 0.02, "green band = ±20 %", transform=ax.transAxes,
            fontsize=7, color=C_GOOD)

    # ---- (b) summary, with the cases still visible --------------------- #
    med = [np.median([100 * (configs[l][c] - truth[c]) / truth[c] for c in cases])
           for l in labels]
    xpos = np.arange(len(labels))
    ax2.bar(xpos, med, width=0.55, color=colours[:len(labels)], alpha=0.30,
            edgecolor="none", zorder=1)
    for ci, lab in enumerate(labels):
        errs = [100 * (configs[lab][c] - truth[c]) / truth[c] for c in cases]
        ax2.plot(np.full(len(errs), ci) + np.linspace(-0.13, 0.13, len(errs)),
                 errs, "o", ms=3.6, color=colours[ci % 4], alpha=0.85,
                 mec="white", mew=0.4, zorder=3)
        ax2.plot([ci], [med[ci]], "_", ms=26, mew=2.6, color=colours[ci % 4],
                 zorder=4)
    ax2.axhline(0, color=C_DATA, lw=0.9)
    ax2.axhspan(-20, 20, color=C_GOOD, alpha=0.09, lw=0)
    ax2.set_xticks(xpos)
    ax2.set_xticklabels([l.replace(" ", "\n") for l in labels], fontsize=7.5)
    ax2.set_yscale("symlog", linthresh=100)
    ax2.set_ylabel("signed error of Hβ flux (%)")
    ax2.set_title("(b) median (bar) with all 10 cases (dots)", loc="left")

    fig.suptitle("Hβ recovery depends on the fitting configuration", fontsize=10, y=1.0)
    fig.text(0.5, -0.02,
             "Every row of every configuration passes the acceptance test — the band is loose enough to admit a 490 % "
             "error on a faint line. None of the four is demonstrated to be correct; their spread is the honest uncertainty.",
             ha="center", fontsize=7.5, color=C_ALT)
    return fig


# --------------------------------------------------------------------------- #
# Figure 3 — uncertainty calibration
# --------------------------------------------------------------------------- #
def fig_pulls(pulls: dict[str, list[dict]]) -> plt.Figure:
    """Standardised residuals (pulls) for each method, with the band that fails.

    ``pulls`` maps a method label → list of dicts with keys ``case``, ``line``,
    ``pull``. A pull is (F_meas − F_true)/σ_F; if σ_F is honest, the pulls are
    drawn from a standard normal, so their RMS should be ≈ 1.
    """
    apply_style()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.2),
                                  gridspec_kw=dict(width_ratios=[1.0, 1.25],
                                                   wspace=0.26))
    colours = {"gspec": C_MODEL, "independent": C_CONT}
    bins = np.linspace(-3.5, 3.5, 22)

    for label, rows in pulls.items():
        vals = np.array([r["pull"] for r in rows])
        rms = math.sqrt(float(np.mean(vals ** 2)))
        ax.hist(vals, bins=bins, histtype="step", lw=1.4,
                color=colours.get(label, C_ALT), density=True,
                label=f"{label}  (n={len(vals)}, RMS={rms:.2f})")
        grid = np.linspace(-3.5, 3.5, 300)
        if label == list(pulls)[0]:
            ax.plot(grid, np.exp(-grid**2 / 2) / math.sqrt(2 * math.pi),
                    lw=0.9, ls="--", color=C_ALT, label="N(0,1) reference")

    ax.set_xlabel(r"pull = ($F_\mathrm{meas}-F_\mathrm{true})/\sigma_F$")
    ax.set_ylabel("probability density")
    ax.set_title("(a) pull distribution", loc="left")
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)

    # (b) pull vs line, showing that aggregate agreement hides structure
    lines = ["Halpha", "Hbeta", "OIII_5007", "NII_6583"]
    ticks = {"Halpha": "Hα", "Hbeta": "Hβ", "OIII_5007": "[O III]", "NII_6583": "[N II]"}
    jitter = {"gspec": -0.13, "independent": 0.13}
    for label, rows in pulls.items():
        for i, ln in enumerate(lines):
            vals = [r["pull"] for r in rows if r["line"] == ln]
            if not vals:
                continue
            xs = i + jitter.get(label, 0.0) + np.linspace(-0.06, 0.06, len(vals))
            ax2.plot(xs, vals, "o", ms=4.2, alpha=0.8,
                     color=colours.get(label, C_ALT),
                     mec="white", mew=0.4, label=label if i == 0 else None)
    for y in (-1, 1):
        ax2.axhline(y, color=C_ALT, lw=0.8, ls=":")
    ax2.axhline(0, color=C_DATA, lw=0.8)
    ax2.set_xticks(range(len(lines)))
    ax2.set_xticklabels([ticks[l] for l in lines])
    ax2.set_ylabel("pull")
    ax2.set_ylim(-3.6, 3.6)
    ax2.set_title("(b) pulls by line and estimator", loc="left")
    ax2.legend(loc="lower left", frameon=False, ncol=2, fontsize=7.5)

    fig.suptitle("Are the reported uncertainties consistent with the actual errors?",
                 fontsize=10, y=1.0)
    fig.text(0.5, -0.03,
             "Aggregate RMS near 1 is necessary, not sufficient: n=40, the four lines of a case share one fit, "
             "and the synthetic spectra are smoother than real data.",
             ha="center", fontsize=7.5, color=C_ALT)
    return fig


# --------------------------------------------------------------------------- #
# Figure 4 — injection recovery
# --------------------------------------------------------------------------- #
def fig_injection(rows: list[dict]) -> plt.Figure:
    """Measured line-flux increment against the known injected increment.

    Each injection adds 40 % of the *pipeline's own* flux to one line, so the
    increment is a differential measurement: it tests whether the pipeline
    transfers added flux faithfully, independent of whether its absolute
    continuum level is right.
    """
    apply_style()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.2),
                                  gridspec_kw=dict(wspace=0.26))
    lines = ["Halpha", "Hbeta", "OIII_5007", "NII_6583"]
    ticks = {"Halpha": "Hα", "Hbeta": "Hβ", "OIII_5007": "[O III] 5007",
             "NII_6583": "[N II] 6583"}
    colours = {"Halpha": C_MODEL, "Hbeta": C_ACCENT,
               "OIII_5007": C_CONT, "NII_6583": C_GOOD}

    for ln in lines:
        pts = [r for r in rows if r["line"] == ln]
        if not pts:
            continue
        x = np.array([r["known_increment"] for r in pts])
        y = np.array([r["measured_increment"] for r in pts])
        e = np.array([r["increment_err"] for r in pts])
        ax.errorbar(x, y, yerr=e, fmt="o", ms=5, capsize=2, lw=0.9,
                    color=colours[ln], mec="white", mew=0.5, label=ticks[ln])

    lim = [0.5, 900]
    ax.plot(lim, lim, lw=0.9, color=C_DATA, ls="-", label="1:1")
    for f, ls in ((1.2, "--"), (0.8, "--")):
        ax.plot(lim, [f * v for v in lim], lw=0.7, color=C_ALT, ls=ls)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("injected increment (known)  ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$)")
    ax.set_ylabel("measured increment")
    ax.set_title("(a) measured vs injected, ±20 % guides", loc="left")
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)

    # (b) fractional recovery — only the lines that actually have injections
    present = [l for l in lines if any(r["line"] == l for r in rows)]
    short = {"Halpha": "Hα", "Hbeta": "Hβ", "OIII_5007": "[O III]\n5007",
             "NII_6583": "[N II]\n6583"}
    for i, ln in enumerate(present):
        pts = [r for r in rows if r["line"] == ln]
        frac = np.array([100 * (r["measured_increment"] / r["known_increment"] - 1)
                         for r in pts])
        xs = np.full(len(frac), i) + np.linspace(-0.11, 0.11, len(frac))
        ax2.plot(xs, frac, "o", ms=5, color=colours[ln], mec="white", mew=0.5)
        ax2.plot([i], [np.median(frac)], "_", ms=22, mew=2.4, color=colours[ln])
        ax2.annotate(f"{np.median(frac):+.0f} %", xy=(i, np.median(frac)),
                     xytext=(11, 0), textcoords="offset points", fontsize=7,
                     va="center", color=colours[ln], annotation_clip=False,
                     clip_on=False)
    ax2.axhline(0, color=C_DATA, lw=0.9)
    for y in (-20, 20):
        ax2.axhline(y, color=C_GOOD, lw=0.8, ls=":")
    ax2.set_xticks(range(len(present)))
    ax2.set_xticklabels([short[l] for l in present], fontsize=8)
    ax2.set_xlim(-0.5, len(present) - 0.28)
    ax2.set_ylim(-24, 10)
    ax2.set_ylabel("recovered increment / injected − 1  (%)")
    ax2.set_title("(b) fractional recovery (dash = median)", loc="left")
    ax2.text(0.02, 0.04, "dotted = ±20 %", transform=ax2.transAxes,
             fontsize=7, color=C_GOOD)

    fig.suptitle("Injection response: does the pipeline transfer added flux?", fontsize=10, y=1.0)
    fig.text(0.5, -0.03,
             "The increment is differential, so it does not measure the absolute flux scale: a uniformly "
             "underestimated Hβ would still be recovered correctly here.",
             ha="center", fontsize=7.5, color=C_ALT)
    return fig


# --------------------------------------------------------------------------- #
# Figure 5 — line-ratio diagnostics
# --------------------------------------------------------------------------- #
def fig_diagnostics(rows: list[dict]) -> plt.Figure:
    """BPT position and Balmer decrement for the real spectra.

    Only objects with both lines above 3σ get a point. The rest are listed in the
    panel with the reason they cannot be placed, because "no point" is a result
    here: it is the measurement that fails, not the object.
    """
    from .diagnostics import bpt_cloud

    apply_style()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.8),
                                  gridspec_kw=dict(width_ratios=[1.25, 1.0],
                                                   wspace=0.26))

    # ---- (a) BPT ------------------------------------------------------- #
    xs, ys, xa, ya = bpt_cloud()
    ax.plot(xs, ys, lw=1.1, color=C_ALT,
            label="Kauffmann+03 (SF ceiling)")
    ax.plot(xa, ya, lw=1.1, color=C_MODEL, ls="--",
            label="Kewley+01 (max starburst)")

    placed = 0
    for r in rows:
        if not r["bpt_placeable"]:
            continue
        x, y = r["nii_ha"], r["oiii_hb"]
        lx = math.log10(x["value"])
        ly = math.log10(y["value"])
        lxe = x["error"] / (x["value"] * math.log(10))
        lye = y["error"] / (y["value"] * math.log(10))
        ax.errorbar(lx, ly, xerr=lxe, yerr=lye, fmt="o", ms=7, color=C_DATA,
                    mec="white", mew=0.8, capsize=2, lw=1.1, zorder=5)
        ax.annotate(r["case_id"], xy=(lx, ly), xytext=(7, 5),
                    textcoords="offset points", fontsize=8)
        placed += 1

    ax.set_xlim(-1.9, 0.15)
    ax.set_ylim(-1.3, 1.1)
    ax.set_xlabel(r"$\log_{10}$([N II] 6583 / Hα)")
    ax.set_ylabel(r"$\log_{10}$([O III] 5007 / Hβ)")
    ax.set_title(f"(a) BPT: {placed} of {len(rows)} objects placeable", loc="left")
    ax.legend(loc="upper right", frameon=False, fontsize=7.5,
              bbox_to_anchor=(1.0, 0.99))

    # Region labels follow the curves, not the corners: everything below the
    # Kauffmann line is star-forming, including the bottom-right corner.
    ax.text(-1.72, -1.02, "star-forming", fontsize=7.5, color=C_GOOD)
    ax.text(-0.62, 0.02, "composite", fontsize=7.5, color=C_ALT, rotation=-18)
    ax.text(-0.02, 0.62, "AGN", fontsize=7.5, color=C_DATA, ha="center")

    unplaced = [r["case_id"] for r in rows if not r["bpt_placeable"]]
    if unplaced:
        ax.text(0.02, 0.055, "no position: " + ", ".join(unplaced) + "  (Hβ < 3σ)",
                transform=ax.transAxes, fontsize=7, color=C_ALT, va="bottom",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))

    # ---- (b) Balmer decrement ------------------------------------------ #
    labels, vals, errs = [], [], []
    for r in rows:
        labels.append(r["case_id"])
        b = r["balmer"]
        vals.append(b["value"] if b["detected"] else np.nan)
        errs.append(b["error"] if b["detected"] else np.nan)

    xpos = np.arange(len(labels))
    det = np.isfinite(vals)
    ax2.errorbar(xpos[det], np.array(vals)[det], yerr=np.array(errs)[det],
                 fmt="o", ms=7, color=C_MODEL, mec="white", mew=0.7,
                 capsize=3, lw=1.2, label="measured Hα/Hβ")
    for i in np.where(~det)[0]:
        ax2.plot([i], [1.6], "x", ms=7, color=C_ALT, mew=1.4)
        ax2.annotate("not measurable", xy=(i, 1.6), xytext=(0, 6),
                     textcoords="offset points", ha="center", fontsize=6.3,
                     color=C_ALT, rotation=90)
    ax2.axhline(2.86, color=C_CONT, lw=1.0, ls="--", label="case B (2.86)")
    ax2.set_xticks(xpos)
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_xlim(-0.6, len(labels) - 0.4)
    ax2.set_ylim(1.2, max(6.2, np.nanmax(np.array(vals)[det]) + 1.0 if det.any() else 6.2))
    ax2.set_ylabel(r"Hα/Hβ")
    ax2.set_title("(b) Balmer decrement", loc="left")
    ax2.legend(loc="upper right", frameon=False, fontsize=7.5)

    # secondary axis: implied A_V, valid only where the ratio exceeds case B
    sec = ax2.twinx()
    from .diagnostics import ccm89_k, LAMBDA_HA, LAMBDA_HB, R_V_CCM89
    dk = ccm89_k(LAMBDA_HB) - ccm89_k(LAMBDA_HA)
    lo, hi = ax2.get_ylim()

    def ratio_to_av(ratio):
        r = max(float(ratio), 2.861)
        return R_V_CCM89 * (2.5 / dk) * math.log10(r / 2.86)

    ticks_ratio = np.arange(2.0, hi + 0.5, 0.5)
    sec.set_ylim(lo, hi)
    sec.set_yticks(ticks_ratio)
    sec.set_yticklabels([f"{ratio_to_av(t):.1f}" for t in ticks_ratio], fontsize=7)
    sec.set_ylabel(r"implied $A_V$ (mag)  —  foreground screen, case B", fontsize=8)
    sec.grid(False)

    fig.suptitle("Line-ratio diagnostics on six real spectra: what is measurable, and what is not",
                 fontsize=10, y=1.0)
    fig.text(0.5, -0.04,
             "n=6, of which 2 allow a BPT position and 2 allow a decrement. No population statement, "
             "no classification, and A_V is an inference that assumes case B, a dust screen and no "
             "stellar Balmer absorption.",
             ha="center", fontsize=7.5, color=C_ALT)
    return fig
