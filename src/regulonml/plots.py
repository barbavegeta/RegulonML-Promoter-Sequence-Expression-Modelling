"""Figures for the README / results folder."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#b9b8b2", "#e7e6e2", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def per_locus(pe: pd.DataFrame, full_key, ctx_key, path):
    """Held-out Spearman per locus: sequence context only vs context + motif features."""
    a = pe[pe.model == full_key].groupby("locus").spearman.mean()
    b = pe[pe.model == ctx_key].groupby("locus").spearman.mean()
    order = a.sort_values().index
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(order) + 1.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    ax.grid(axis="y", visible=False)
    ax.hlines(y, b[order], a[order], color=MUTED, lw=1.5, zorder=1)
    ax.scatter(b[order], y, s=40, c=MUTED, edgecolors=SURFACE, linewidths=1.2, zorder=2, label="context only")
    ax.scatter(a[order], y, s=48, c=BLUE, edgecolors=SURFACE, linewidths=1.2, zorder=3, label="context + TF motifs")
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(y, order, fontsize=8.5, color=INK)
    ax.set_xlabel("Spearman ρ, predicted vs measured effect (locus held out)", color=INK2, fontsize=9.5)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right", labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def summary(metrics: pd.DataFrame, path):
    """Median held-out Spearman and AUROC for each model; the leaky random split shown for contrast."""
    m = metrics.copy()
    m["label"] = np.where(m.model == "Mean baseline", "Mean baseline", m.model + "\n" + m.features)
    m.loc[m.split.str.startswith("random"), "label"] = "Gradient boosting\ncontext + motifs\nRANDOM split (leaky)"
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    for ax, col, title, ref in ((axes[0], "median_element_spearman", "Median Spearman ρ per element", 0),
                                (axes[1], "pooled_AUROC_sig", "AUROC: significant variants (p < 1e-5)", 0.5)):
        _style(ax)
        d = m.dropna(subset=[col])
        cols = [MUTED if "leaky" in lbl else (INK2 if lbl == "Mean baseline" else BLUE) for lbl in d.label]
        hatch = ["//" if "leaky" in lbl else "" for lbl in d.label]
        bars = ax.barh(np.arange(len(d)), d[col] - ref, left=ref, color=cols, height=0.6, edgecolor=SURFACE)
        for b, h in zip(bars, hatch):
            b.set_hatch(h)
        for i, v in enumerate(d[col]):
            ax.text(max(v, ref) + 0.005, i, f"{v:.2f}", va="center", fontsize=8.5, color=INK)
        ax.set_yticks(np.arange(len(d)), d.label, fontsize=8, color=INK)
        ax.invert_yaxis()
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.axvline(ref, color=INK2, lw=1)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def tert_case(F: pd.DataFrame, pred_key, path, element="TERT-GBM", hotspots=(1295113, 1295135)):
    """Measured effects along the TERT promoter, with the two cancer hotspots and their ETS-site gain."""
    t = F[F.Element == element]
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.2), dpi=150, sharex=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        _style(ax)
    hs = t.Pos.isin(hotspots) & (t.Alt == "A")
    axes[0].scatter(t.Pos[~hs], t.effect[~hs], s=8, c=MUTED, edgecolors="none", label="all variants")
    axes[0].scatter(t.Pos[hs], t.effect[hs], s=60, c=ORANGE, edgecolors=SURFACE, linewidths=1.2, zorder=3,
                    label="C228T / C250T (G>A on + strand)")
    axes[0].axhline(0, color=INK2, lw=0.8)
    axes[0].set_ylabel("measured effect\n(log2, MPRA)", color=INK2, fontsize=9)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
    axes[0].set_title(f"{element}: saturation mutagenesis of the TERT promoter", loc="left", fontsize=10, color=INK)
    axes[1].scatter(t.Pos[~hs], t["gain__Ets_related"][~hs], s=8, c=MUTED, edgecolors="none")
    axes[1].scatter(t.Pos[hs], t["gain__Ets_related"][hs], s=60, c=ORANGE, edgecolors=SURFACE, linewidths=1.2, zorder=3)
    axes[1].set_ylabel("ETS-family site created\n(relative motif score gain)", color=INK2, fontsize=9)
    axes[1].set_xlabel("Position (GRCh38, chr5)", color=INK2, fontsize=9)
    axes[1].ticklabel_format(useOffset=False, style="plain", axis="x")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
