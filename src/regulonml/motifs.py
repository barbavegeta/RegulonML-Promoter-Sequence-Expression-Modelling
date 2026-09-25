"""Score how each variant changes transcription-factor binding-site strength (JASPAR PWMs).

For a variant at position i the relevant binding sites are the motif windows
that overlap i, on either strand. The change is
``max(alt window scores) - max(ref window scores)`` in log2 odds units.
SNVs are handled fully vectorised; single-base deletions use prefix/suffix
sums so the shifted windows are scored without re-scanning the sequence.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BASES = "ACGT"
IDX = {b: i for i, b in enumerate(BASES)}


@dataclass
class Motif:
    matrix_id: str
    name: str
    family: str
    pwm: np.ndarray  # 4 x L log2-odds (A, C, G, T)

    @property
    def length(self):
        return self.pwm.shape[1]

    @property
    def min_score(self):
        return float(self.pwm.min(axis=0).sum())

    @property
    def max_score(self):
        return float(self.pwm.max(axis=0).sum())


def load_jaspar(release="JASPAR2026", collection="CORE", tax_group="vertebrates", pseudocount=0.5):
    """JASPAR PWMs from the pyjaspar package (the database ships with it: no download)."""
    from pyjaspar import jaspardb

    ms = jaspardb(release=release).fetch_motifs(collection=collection, tax_group=[tax_group])
    out = []
    for m in ms:
        pssm = m.counts.normalize(pseudocounts=pseudocount).log_odds()
        pwm = np.array([pssm[b] for b in BASES], dtype=float)
        fam = (m.tf_family[0] if m.tf_family else (m.tf_class[0] if m.tf_class else "Unclassified"))
        out.append(Motif(m.matrix_id, m.name, fam, pwm))
    return out


def encode(seq: str) -> np.ndarray:
    """Sequence -> int array (A0 C1 G2 T3, N -> -1)."""
    return np.array([IDX.get(c, -1) for c in seq.upper()], dtype=int)


def _window_scores(x: np.ndarray, pwm: np.ndarray):
    """Per-window cumulative contributions C[s, t] and totals for windows fully inside x."""
    L = pwm.shape[1]
    n = len(x) - L + 1
    if n <= 0:
        return np.zeros((0, L)), np.zeros(0)
    pos = np.arange(n)[:, None] + np.arange(L)[None, :]
    bases = x[pos]
    contrib = np.where(bases >= 0, pwm[np.clip(bases, 0, 3), np.arange(L)[None, :]], pwm.min(axis=0)[None, :])
    C = np.cumsum(contrib, axis=1)
    return C, C[:, -1]


def variant_deltas(seq: str, motif: Motif):
    """Motif score change for every SNV and single-base deletion in ``seq``.

    Returns (snv_delta [n, 4], snv_best [n, 4], del_delta [n], del_best [n], ref_best [n]) where
    *_best is max(ref, alt) window score (used to decide whether a site is present) and
    ref_best is the best reference window covering each position.
    Rows are positions; SNV columns are alt bases A, C, G, T (the ref column is 0).
    """
    x = encode(seq)
    n = len(x)
    out_snv = np.full((n, 4), -np.inf)
    best_snv = np.full((n, 4), -np.inf)
    out_del = np.full(n, -np.inf)
    best_del = np.full(n, -np.inf)
    first = True
    for pwm in (motif.pwm, motif.pwm[::-1, ::-1]):  # forward and reverse-complement
        L = pwm.shape[1]
        C, tot = _window_scores(x, pwm)
        nwin = len(tot)
        if nwin == 0:
            continue
        # covering windows for position i: start s = i - j, offset j = 0..L-1
        i = np.arange(n)[:, None]
        j = np.arange(L)[None, :]
        s = i - j
        valid = (s >= 0) & (s < nwin)
        s_c = np.clip(s, 0, nwin - 1)
        ref_sc = np.where(valid, tot[s_c], -np.inf)                       # [n, L]
        ref_best = ref_sc.max(axis=1)                                        # [n]
        refb = np.clip(x, 0, 3)[:, None]
        # alt score = ref window score + pwm[alt, j] - pwm[ref, j]
        alt_sc = ref_sc[:, :, None] + pwm.T[j[0]][None, :, :] - pwm[refb, j][:, :, None]  # [n, L, 4]
        alt_best = alt_sc.max(axis=1)                                         # [n, 4]
        # single-base deletion at i: windows starting at a = i - m (m = 1..L-1) that span the junction
        m = np.arange(1, L)[None, :]
        a = i - m
        ok = (a >= 0) & (a + 1 < nwin)
        a_c = np.clip(a, 0, nwin - 2 if nwin > 1 else 0)
        pre = C[a_c, m - 1]                                                   # offsets < m from window a
        suf = tot[np.clip(a_c + 1, 0, nwin - 1)] - C[np.clip(a_c + 1, 0, nwin - 1), m - 1]
        del_sc = np.where(ok, pre + suf, -np.inf)
        del_best_s = del_sc.max(axis=1)
        if first:
            ref_all, alt_all, del_all = ref_best, alt_best, del_best_s
            first = False
        else:
            ref_all = np.maximum(ref_all, ref_best)
            alt_all = np.maximum(alt_all, alt_best)
            del_all = np.maximum(del_all, del_best_s)
    if first:
        return out_snv, best_snv, out_del, best_del, np.full(n, -np.inf)
    finite = np.isfinite(ref_all)
    out_snv[finite] = alt_all[finite] - ref_all[finite, None]
    best_snv[finite] = np.maximum(alt_all[finite], ref_all[finite, None])
    okd = finite & np.isfinite(del_all)
    out_del[okd] = del_all[okd] - ref_all[okd]
    best_del[okd] = np.maximum(del_all[okd], ref_all[okd])
    # the reference base itself is not a variant
    rows = np.where(x >= 0)[0]
    out_snv[rows, x[rows]] = 0.0
    return out_snv, best_snv, out_del, best_del, ref_all


def brute_force_delta(seq: str, motif: Motif, i: int, alt: str) -> float:
    """Slow reference implementation used by the tests (alt '-' = deletion)."""
    def best_covering(s, lo, hi):
        best = -np.inf
        for pwm in (motif.pwm, motif.pwm[::-1, ::-1]):
            L = pwm.shape[1]
            for st in range(max(0, lo - L + 1), min(hi, len(s) - L) + 1):
                if st + L - 1 < lo or st > hi:
                    continue
                xx = encode(s[st:st + L])
                if len(xx) < L:
                    continue
                best = max(best, sum(pwm[b, k] if b >= 0 else pwm[:, k].min() for k, b in enumerate(xx)))
        return best
    ref = best_covering(seq, i, i)
    if alt == "-":
        s2 = seq[:i] + seq[i + 1:]
        # windows must span the junction (contain both i-1 and i in the new string)
        best = -np.inf
        for pwm in (motif.pwm, motif.pwm[::-1, ::-1]):
            L = pwm.shape[1]
            for st in range(max(0, i - L + 1), i):
                if st + L <= len(s2) and st + L - 1 >= i:
                    xx = encode(s2[st:st + L])
                    best = max(best, sum(pwm[b, k] for k, b in enumerate(xx)))
        altb = best
    else:
        altb = best_covering(seq[:i] + alt + seq[i + 1:], i, i)
    return altb - ref
