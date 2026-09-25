"""Tests for sequence reconstruction, motif scoring, feature hygiene and the CV split."""
import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import LeaveOneGroupOut

from regulonml.data import LOCUS, element_sequences
from regulonml.features import build_features
from regulonml.motifs import IDX, Motif, brute_force_delta, variant_deltas
from regulonml.pipeline import LEAKY, feature_sets

RNG = np.random.default_rng(7)


def _random_motif(L=8, name="M", family="F"):
    p = RNG.dirichlet(np.ones(4) * 0.5, size=L).T
    return Motif(name, name, family, np.log2(p / 0.25))


def _consensus_motif(consensus, name="ETS", family="Ets-related"):
    p = np.full((4, len(consensus)), 0.02)
    for j, b in enumerate(consensus):
        p[IDX[b], j] = 0.94
    return Motif(name, name, family, np.log2(p / 0.25))


def test_vectorised_deltas_match_brute_force():
    seq = "".join(RNG.choice(list("ACGT"), 80))
    for m in [_random_motif(L) for L in (5, 8, 12)]:
        snv, _, dl, _, _ = variant_deltas(seq, m)
        for i in (0, 1, 7, 40, 78, 79):
            for k, b in enumerate("ACGT"):
                if b != seq[i]:
                    assert snv[i, k] == pytest.approx(brute_force_delta(seq, m, i, b))
            bf = brute_force_delta(seq, m, i, "-")
            assert (np.isinf(bf) and np.isinf(dl[i])) or dl[i] == pytest.approx(bf)


def test_site_creation_and_destruction_on_both_strands():
    m = _consensus_motif("CCGGAAGT")
    seq = "TTTTCCGGTAGTTTTT"             # one base away from CCGGAAGT at index 4
    snv, *_ = variant_deltas(seq, m)
    assert snv[8, IDX["A"]] > 5            # T>A creates the site (gain)
    rc = "AAAAACTACCGGAAAA"               # the same site on the reverse strand
    snv_rc, *_ = variant_deltas(rc, m)
    assert snv_rc[7, IDX["T"]] > 5         # A>T creates it on the minus strand
    site = "TTTTCCGGAAGTTTTT"
    snv_s, *_ = variant_deltas(site, m)
    assert snv_s[7, IDX["C"]] < -5         # breaking the core GGAA is a loss


def _toy_data():
    rows = []
    for e, s in {"TERT-GBM": "ACGTTGCAACGT", "TERT-HEK": "ACGTTGCAACGT", "HBB": "GGGCCCAAATTT"}.items():
        for i, ref in enumerate(s):
            for alt in "ACGT-":
                if alt != ref:
                    rows.append({"Element": e, "Pos": 100 + i, "Ref": ref, "Alt": alt, "effect": RNG.normal(),
                                 "p": RNG.uniform(), "Barcodes": 50, "DNA": 1, "RNA": 1})
    d = pd.DataFrame(rows)
    d["locus"] = d.Element.map(LOCUS).fillna(d.Element)
    d["region_type"] = "promoter"
    d["is_deletion"] = (d.Alt == "-").astype(int)
    return d


def test_sequence_reconstruction():
    seqs = element_sequences(_toy_data())
    assert seqs["TERT-GBM"] == (100, "ACGTTGCAACGT")
    assert seqs["HBB"][1] == "GGGCCCAAATTT"


def test_features_exclude_measurements_and_split_keeps_loci_together():
    d = _toy_data()
    ctx, mot, _ = build_features(d, element_sequences(d), [_random_motif(6, "A", "F1"), _random_motif(9, "B", "F2")])
    for cols in feature_sets(ctx, mot).values():
        assert not LEAKY & set(cols)
    for tr, te in LeaveOneGroupOut().split(d, groups=d.locus):
        assert set(d.locus.iloc[tr]).isdisjoint(d.locus.iloc[te])
    assert d[d.Element.str.startswith("TERT")].locus.nunique() == 1
