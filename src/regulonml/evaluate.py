"""Leave-one-locus-out evaluation of variant-effect models.

Each fold holds out every element measured at one locus (e.g. all four TERT
cell-line datasets), so the model is always scored on DNA it has never seen.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import average_precision_score, r2_score, roc_auc_score
from sklearn.model_selection import KFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SIG_P = 1e-5  # significance threshold used by Kircher et al. 2019


def regressors(seed=0):
    return {
        "Ridge": make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 13))),
        "Gradient boosting": HistGradientBoostingRegressor(max_iter=250, learning_rate=0.08, max_leaf_nodes=31,
                                                           min_samples_leaf=40, l2_regularization=1.0,
                                                           random_state=seed),
    }


def classifiers(seed=0):
    return {
        "Logistic": make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=500)),
        "Gradient boosting": HistGradientBoostingClassifier(max_iter=250, learning_rate=0.08, max_leaf_nodes=31,
                                                            min_samples_leaf=40, l2_regularization=1.0,
                                                            random_state=seed),
    }


def _fit_predict(model_factory, X, y, tr, te, proba):
    m = model_factory()
    m.fit(X.iloc[tr], y.iloc[tr])
    return te, (m.predict_proba(X.iloc[te])[:, 1] if proba else m.predict(X.iloc[te]))


def oof_predict(model_factory, X, y, groups=None, cv=None, proba=False, n_jobs=-1):
    """Out-of-fold predictions; ``groups`` -> leave-one-group-out. Folds run in parallel."""
    from joblib import Parallel, delayed

    pred = np.full(len(y), np.nan)
    splitter = LeaveOneGroupOut() if groups is not None else cv
    jobs = Parallel(n_jobs=n_jobs)(delayed(_fit_predict)(model_factory, X, y, tr, te, proba)
                                   for tr, te in splitter.split(X, y, groups))
    for te, p in jobs:
        pred[te] = p
    return pred


def per_element_spearman(d, pred):
    out = []
    for e, idx in d.groupby("Element").groups.items():
        rows = d.index.get_indexer(idx)
        rho = spearmanr(d.effect.values[rows], pred[rows]).statistic
        out.append({"Element": e, "locus": d.locus.iloc[rows[0]], "n": len(rows), "spearman": rho})
    return pd.DataFrame(out)


def per_element_auroc(d, score):
    out = []
    for e, idx in d.groupby("Element").groups.items():
        rows = d.index.get_indexer(idx)
        y = (d.p.values[rows] < SIG_P).astype(int)
        if 0 < y.sum() < len(y):
            out.append({"Element": e, "auroc": roc_auc_score(y, score[rows]), "prevalence": y.mean()})
    return pd.DataFrame(out)


def evaluate(d, feature_sets: dict, seed=0):
    """Returns (summary table, per-element table, out-of-fold predictions)."""
    y = d.effect
    y_sig = (d.p < SIG_P).astype(int)
    groups = d.locus.values
    rows, per_el, preds = [], [], {}

    base = np.full(len(d), np.nan)
    for tr, te in LeaveOneGroupOut().split(d, y, groups):
        base[te] = y.iloc[tr].mean()
    preds["Mean baseline"] = base
    rows.append({"features": "-", "model": "Mean baseline", "split": "leave-one-locus-out",
                 "pooled_R2": r2_score(y, base), "median_element_spearman": np.nan,
                 "pooled_AUROC_sig": 0.5, "pooled_AUPRC_sig": y_sig.mean()})

    for fs_name, cols in feature_sets.items():
        X = d[cols]
        for mname, _ in regressors(seed).items():
            factory = lambda m=mname: regressors(seed)[m]  # noqa: E731
            p = oof_predict(factory, X, y, groups=groups)
            key = f"{mname} | {fs_name}"
            preds[key] = p
            pe = per_element_spearman(d, p).assign(model=key)
            per_el.append(pe)
            cname = "Gradient boosting" if mname == "Gradient boosting" else "Logistic"
            cf = lambda m=cname: classifiers(seed)[m]  # noqa: E731
            ps = oof_predict(cf, X, y_sig, groups=groups, proba=True)
            preds[f"{key} | P(sig)"] = ps
            rows.append({"features": fs_name, "model": mname, "split": "leave-one-locus-out",
                         "pooled_R2": r2_score(y, p),
                         "median_element_spearman": pe.spearman.median(),
                         "pooled_AUROC_sig": roc_auc_score(y_sig, ps),
                         "pooled_AUPRC_sig": average_precision_score(y_sig, ps)})

    # optimistic reference: random split, the same locus appears in train and test
    full = list(feature_sets.values())[-1]
    p = oof_predict(lambda: regressors(seed)["Gradient boosting"], d[full], y,
                    cv=KFold(5, shuffle=True, random_state=seed))
    pe = per_element_spearman(d, p)
    rows.append({"features": list(feature_sets)[-1], "model": "Gradient boosting",
                 "split": "random 5-fold (leaky, for comparison)", "pooled_R2": r2_score(y, p),
                 "median_element_spearman": pe.spearman.median(), "pooled_AUROC_sig": np.nan,
                 "pooled_AUPRC_sig": np.nan})
    return pd.DataFrame(rows), pd.concat(per_el, ignore_index=True), pd.DataFrame(preds, index=d.index)
