"""Sequence-only features for each variant. No read counts, no element identity."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .motifs import IDX, variant_deltas

# Scores are put on each motif's relative scale, (score - min) / (max - min), so motifs of
# different length and information content are comparable.
GAIN_SITE = 0.90   # a gain counts only if the variant creates a strong site (alt >= 0.90)
LOSS_SITE = 0.85   # a loss counts only if the variant hits an existing site (ref >= 0.85)
DELTA_MIN = 0.05   # relative change counted in n_gain / n_loss
CONTEXT = 5            # bases either side one-hot encoded


def _slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")[:40]


def motif_features(d: pd.DataFrame, seqs: dict, motifs) -> pd.DataFrame:
    """Per-variant motif gain/loss, aggregated by TF family plus global summaries."""
    families = sorted({m.family for m in motifs})
    fam_idx = {f: k for k, f in enumerate(families)}
    n = len(d)
    gain = np.zeros((n, len(families)))
    loss = np.zeros((n, len(families)))
    g_max = np.zeros(n); l_min = np.zeros(n); n_gain = np.zeros(n); n_loss = np.zeros(n)
    # reference-only annotation: is this base inside a strong predicted binding site?
    n_sites = np.zeros(n); best_site = np.zeros(n)
    top_gain_motif = np.array([""] * n, dtype=object)

    cache = {}
    for e, rows in d.groupby("Element").groups.items():
        start, seq = seqs[e]
        sub = d.loc[rows]
        pos = (sub.Pos.values - start).astype(int)
        alt = sub.Alt.values
        is_del = alt == "-"
        col = np.array([IDX.get(a, 0) for a in alt])
        out_rows = d.index.get_indexer(rows)
        for m in motifs:
            key = (seq, m.matrix_id)
            if key not in cache:
                cache[key] = variant_deltas(seq, m)
            snv, best_snv, dl, best_dl, ref_best = cache[key]
            delta = np.where(is_del, dl[pos], snv[pos, col])
            rng = m.max_score - m.min_score
            ref_rel = np.nan_to_num((ref_best[pos] - m.min_score) / rng, neginf=0.0)
            alt_rel = np.nan_to_num((ref_best[pos] + delta - m.min_score) / rng, neginf=0.0, nan=0.0)
            d_rel = alt_rel - ref_rel
            dlt = np.where(d_rel > 0, np.where(alt_rel >= GAIN_SITE, d_rel, 0.0),
                           np.where(ref_rel >= LOSS_SITE, d_rel, 0.0))
            n_sites[out_rows] += ref_rel >= LOSS_SITE
            best_site[out_rows] = np.maximum(best_site[out_rows], ref_rel)
            f = fam_idx[m.family]
            gain[out_rows, f] = np.maximum(gain[out_rows, f], np.maximum(dlt, 0))
            loss[out_rows, f] = np.minimum(loss[out_rows, f], np.minimum(dlt, 0))
            better = dlt > g_max[out_rows]
            top_gain_motif[out_rows[better]] = m.name
            g_max[out_rows] = np.maximum(g_max[out_rows], dlt)
            l_min[out_rows] = np.minimum(l_min[out_rows], dlt)
            n_gain[out_rows] += dlt >= DELTA_MIN
            n_loss[out_rows] += dlt <= -DELTA_MIN
    cols = {}
    for f, k in fam_idx.items():
        cols[f"gain__{_slug(f)}"] = gain[:, k]
        cols[f"loss__{_slug(f)}"] = loss[:, k]
    feats = pd.DataFrame(cols, index=d.index)
    feats = feats.loc[:, feats.abs().sum() > 0]  # drop families never touched
    feats["motif_max_gain"] = g_max
    feats["motif_max_loss"] = l_min
    feats["motif_n_gain"] = n_gain
    feats["motif_n_loss"] = n_loss
    feats["ref_n_sites"] = n_sites
    feats["ref_best_site"] = best_site
    # predicted sites per base over a 21-bp window (reference only): marks dense TF clusters
    per_pos = (pd.DataFrame({"Element": d.Element.values, "Pos": d.Pos.values, "s": n_sites})
               .drop_duplicates(["Element", "Pos"]).sort_values(["Element", "Pos"]))
    per_pos["dens"] = per_pos.groupby("Element").s.transform(
        lambda v: v.rolling(21, center=True, min_periods=1).mean())
    feats["ref_site_density_21bp"] = d[["Element", "Pos"]].merge(
        per_pos[["Element", "Pos", "dens"]], on=["Element", "Pos"], how="left")["dens"].values
    return feats, pd.Series(top_gain_motif, index=d.index, name="top_gain_motif")


def context_features(d: pd.DataFrame, seqs: dict) -> pd.DataFrame:
    rows = []
    for e, sub in d.groupby("Element"):
        start, seq = seqs[e]
        L = len(seq)
        for idx, r in sub.iterrows():
            i = r.Pos - start
            f = {"_idx": idx, "rel_position": i / max(L - 1, 1),
                 "gc_25bp": (lambda w: (w.count("G") + w.count("C")) / max(len(w), 1))(seq[max(0, i - 25): i + 26])}
            for k in range(-CONTEXT, CONTEXT + 1):
                if k == 0:
                    continue
                b = seq[i + k] if 0 <= i + k < L else "N"
                for base in "ACGT":
                    f[f"ctx{k:+d}_{base}"] = int(b == base)
            ref, alt = r.Ref, r.Alt
            prev_b = seq[i - 1] if i > 0 else "N"
            next_b = seq[i + 1] if i + 1 < L else "N"
            f["cpg_ref"] = int((ref == "C" and next_b == "G") or (ref == "G" and prev_b == "C"))
            f["cpg_alt"] = int((alt == "C" and next_b == "G") or (alt == "G" and prev_b == "C"))
            rows.append(f)
    c = pd.DataFrame(rows).set_index("_idx").reindex(d.index)
    for b in "ACGT":
        c[f"ref_{b}"] = (d.Ref == b).astype(int)
        c[f"alt_{b}"] = (d.Alt == b).astype(int)
    c["is_deletion"] = d.is_deletion
    trans = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    c["is_transition"] = [int((r, a) in trans) for r, a in zip(d.Ref, d.Alt)]
    c["is_promoter"] = (d.region_type == "promoter").astype(int)
    return c


def build_features(d, seqs, motifs):
    ctx = context_features(d, seqs)
    mot, top = motif_features(d, seqs, motifs)
    return ctx, mot, top
