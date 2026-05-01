import pandas as pd

ALLELES = ["A", "C", "G", "T", "-"]

def encode_allele_features(df: pd.DataFrame) -> pd.DataFrame:
    """Create model-ready features from MPRA variant-level data."""
    features = pd.DataFrame(index=df.index)

    numeric_cols = [
        "Tags", "DNA", "RNA", "is_deletion", "log_tags", "log_dna",
        "log_rna", "rna_dna_ratio", "position_scaled_within_element"
    ]

    for col in numeric_cols:
        if col in df.columns:
            features[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    for allele in ALLELES:
        features[f"ref_{allele}"] = (df["Ref"].astype(str) == allele).astype(int)
        features[f"alt_{allele}"] = (df["Alt"].astype(str) == allele).astype(int)

    # Basic substitution class features.
    ref = df["Ref"].astype(str)
    alt = df["Alt"].astype(str)
    features["is_transition"] = (
        ((ref == "A") & (alt == "G")) |
        ((ref == "G") & (alt == "A")) |
        ((ref == "C") & (alt == "T")) |
        ((ref == "T") & (alt == "C"))
    ).astype(int)
    features["is_transversion"] = ((alt != "-") & (features["is_transition"] == 0) & (ref != alt)).astype(int)

    return features.fillna(0.0)
