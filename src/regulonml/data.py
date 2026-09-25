"""Load the Kircher et al. (2019) saturation-mutagenesis MPRA data and rebuild element sequences.

Every position of every element was mutated to all alternatives, so each
element's reference sequence can be reconstructed from the ``Ref`` alleles
without downloading a genome.
"""
from __future__ import annotations

import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

URL = "https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/master/data/elements.tsv.gz"

# Several elements are the same DNA tested in different cells, time points or
# orientations. They must stay together when splitting train/test, otherwise the
# model is tested on a sequence it has already seen.
LOCUS = {
    "TERT-GAa": "TERT", "TERT-GBM": "TERT", "TERT-GSc": "TERT", "TERT-HEK": "TERT",
    "PKLR-24h": "PKLR", "PKLR-48h": "PKLR",
    "LDLR": "LDLR", "LDLR.2": "LDLR",
    "SORT1": "SORT1", "SORT1.2": "SORT1", "SORT1-flip": "SORT1",
    "ZRSh-13": "ZRS", "ZRSh-13h2": "ZRS",
}
PROMOTERS = {"F9", "FOXE1", "GP1BA", "HBB", "HBG1", "HNF4A", "LDLR", "MSMB", "PKLR", "TERT"}


def download(path, url=URL):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.stat().st_size == 0:
        urllib.request.urlretrieve(url, path)
    return path


def load_variants(path, release="GRCh38", min_barcodes=10) -> pd.DataFrame:
    """One row per (element, position, alt) with the measured log2 effect.

    Read counts (DNA, RNA, barcodes) are kept only for filtering and QC: the
    effect estimate is computed from them, so they must never be model features.
    """
    d = pd.read_csv(path, sep="\t", dtype={"Chrom": str}, low_memory=False)
    d = d[(d.Release == release) & d.Pos.notna() & (d.Barcodes >= min_barcodes)].copy()
    d["Pos"] = d.Pos.astype(int)
    d = d.rename(columns={"Coefficient": "effect", "pValue": "p"})
    d["locus"] = d.Element.map(LOCUS).fillna(d.Element)
    d["region_type"] = np.where(d.locus.isin(PROMOTERS), "promoter", "enhancer")
    d["is_deletion"] = (d.Alt == "-").astype(int)
    return d.reset_index(drop=True)


def element_sequences(d: pd.DataFrame) -> dict[str, tuple[int, str]]:
    """Reconstruct each element's reference sequence: {element: (start position, sequence)}."""
    out = {}
    for e, x in d.groupby("Element"):
        refs = x.groupby("Pos").Ref.first()
        start, end = int(refs.index.min()), int(refs.index.max())
        seq = "".join(refs.get(p, "N") for p in range(start, end + 1))
        out[e] = (start, seq.upper())
    return out
