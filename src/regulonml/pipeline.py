"""End-to-end run: download -> sequences -> features -> leave-one-locus-out evaluation -> figures.

    python -m regulonml.pipeline --outdir results
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

from . import plots
from .data import download, element_sequences, load_variants
from .evaluate import evaluate
from .features import build_features
from .motifs import load_jaspar

LEAKY = {"Barcodes", "DNA", "RNA", "effect", "p"}  # measurement columns: never features


def feature_sets(ctx: pd.DataFrame, mot: pd.DataFrame):
    c = [x for x in ctx.columns if x != "is_deletion"] + ["is_deletion"]
    sets = {"context only": c, "context + motifs": c + list(mot.columns)}
    for cols in sets.values():
        assert not LEAKY & set(cols), "measurement columns must not be used as features"
    return sets


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", default="data/raw/kircher_elements.tsv.gz")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--min-barcodes", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    raw = download(a.raw)
    d = load_variants(raw, min_barcodes=a.min_barcodes)
    seqs = element_sequences(d)
    print(f"{len(d):,} variants, {d.Element.nunique()} datasets, {d.locus.nunique()} loci")

    feat_file = out / "features.parquet"
    if feat_file.exists():
        F = pd.read_parquet(feat_file)
        ctx_cols = [c for c in F if c.startswith(("ctx", "ref_A", "ref_C", "ref_G", "ref_T", "alt_", "cpg",
                                                  "gc_", "rel_position", "is_"))]
        mot_cols = [c for c in F if c.startswith(("gain__", "loss__", "motif_", "ref_n_sites",
                                                  "ref_best_site", "ref_site_density"))]
        sets = feature_sets(F[ctx_cols], F[mot_cols])
    else:
        motifs = load_jaspar()
        print(f"scoring {len(motifs)} JASPAR motifs ...")
        ctx, mot, top = build_features(d, seqs, motifs)
        F = pd.concat([d, ctx.drop(columns="is_deletion"), mot, top], axis=1)
        F.to_parquet(feat_file, index=False)
        sets = feature_sets(ctx, mot)
    print(f"features: {len(sets['context only'])} context, {len(sets['context + motifs'])} total "
          f"({time.time() - t0:.0f} s)")

    metrics, per_el, preds = evaluate(F, sets, seed=a.seed)
    metrics.to_csv(out / "metrics.csv", index=False)
    per_el.to_csv(out / "per_element_spearman.csv", index=False)
    pd.concat([F[["Element", "locus", "Pos", "Ref", "Alt", "effect", "p"]], preds], axis=1).to_csv(
        out / "oof_predictions.csv.gz", index=False)

    full, ctx_key = "Gradient boosting | context + motifs", "Gradient boosting | context only"
    plots.per_locus(per_el, full, ctx_key, out / "per_locus_spearman.png")
    plots.summary(metrics, out / "model_summary.png")
    plots.tert_case(F, full, out / "tert_case_study.png")

    with open(out / "run_summary.json", "w") as fh:
        json.dump({"variants": len(d), "datasets": int(d.Element.nunique()), "loci": int(d.locus.nunique()),
                   "min_barcodes": a.min_barcodes, "cv": "leave-one-locus-out",
                   "seconds": round(time.time() - t0)}, fh, indent=2)
    print(metrics.round(3).to_string(index=False))
    print(f"results in {out}/ ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
