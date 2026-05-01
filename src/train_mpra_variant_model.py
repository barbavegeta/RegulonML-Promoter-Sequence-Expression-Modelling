import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from mpra_feature_engineering import encode_allele_features

def rmse(y_true, y_pred):
    return float(mean_squared_error(y_true, y_pred) ** 0.5)

def split_data(X, y, df, group_col, test_size, random_state):
    if group_col:
        splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
        train_idx, test_idx = next(splitter.split(X, y, groups=df[group_col]))
        return X.iloc[train_idx], X.iloc[test_idx], y.iloc[train_idx], y.iloc[test_idx], df.iloc[train_idx], df.iloc[test_idx]
    return train_test_split(X, y, df, test_size=test_size, random_state=random_state)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="results/real_mpra")
    parser.add_argument("--target", default="Value")
    parser.add_argument("--group-col", default="Element")
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    y = pd.to_numeric(df[args.target], errors="coerce")
    mask = y.notna()
    df = df.loc[mask].reset_index(drop=True)
    y = y.loc[mask].reset_index(drop=True)

    X = encode_allele_features(df)

    X_train, X_test, y_train, y_test, df_train, df_test = split_data(
        X, y, df, args.group_col, args.test_size, args.random_state
    )

    models = {
        "mean_baseline": DummyRegressor(strategy="mean"),
        "elastic_net": Pipeline([
            ("scale", StandardScaler()),
            ("model", ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=10000, random_state=args.random_state)),
        ]),
        "random_forest": RandomForestRegressor(
            n_estimators=300,
            min_samples_leaf=5,
            random_state=args.random_state,
            n_jobs=-1,
        ),
    }

    metrics = []
    fitted = {}

    for name, model in models.items():
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        metrics.append({
            "model": name,
            "r2": r2_score(y_test, pred),
            "rmse": rmse(y_test, pred),
            "mae": mean_absolute_error(y_test, pred),
        })
        fitted[name] = model

    metrics_df = pd.DataFrame(metrics).sort_values("rmse")
    metrics_df.to_csv(output / "metrics.csv", index=False)

    best_name = metrics_df.iloc[0]["model"]
    best_model = fitted[best_name]
    best_pred = best_model.predict(X_test)

    pred_df = df_test.copy()
    pred_df["observed_value"] = y_test.values
    pred_df["predicted_value"] = best_pred
    pred_df["best_model"] = best_name
    pred_df.to_csv(output / "predictions.csv", index=False)

    rf = fitted["random_forest"]
    pd.DataFrame({
        "feature": X.columns,
        "importance": rf.feature_importances_,
    }).sort_values("importance", ascending=False).to_csv(output / "feature_importance.csv", index=False)

    perm = permutation_importance(
        best_model, X_test, y_test,
        n_repeats=10,
        random_state=args.random_state,
        n_jobs=-1,
        scoring="neg_root_mean_squared_error",
    )
    pd.DataFrame({
        "feature": X.columns,
        "importance_mean": perm.importances_mean,
        "importance_std": perm.importances_std,
    }).sort_values("importance_mean", ascending=False).to_csv(output / "permutation_importance.csv", index=False)

    summary = {
        "input": args.input,
        "target": args.target,
        "group_col": args.group_col,
        "train_rows": int(len(X_train)),
        "test_rows": int(len(X_test)),
        "n_features": int(X.shape[1]),
        "best_model": str(best_name),
    }
    with open(output / "run_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(metrics_df.to_string(index=False))
    print(f"Best model: {best_name}")
    print(f"Results written to: {output}")

if __name__ == "__main__":
    main()
