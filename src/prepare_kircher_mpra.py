import argparse
from pathlib import Path
import numpy as np
import pandas as pd

PROMOTER_ELEMENTS = {
    "F9", "FOXE1", "GP1BA", "HBB", "HBG1", "HNF4A", "LDLR.2", "LDLR",
    "MSMB", "TERT-GBM", "PKLR-24h", "PKLR-48h", "TERT-GAa", "TERT-GSc", "TERT-HEK"
}

COLUMN_RENAME = {
    "Chrom": "Chromosome",
    "Pos": "Position",
    "Coefficient": "Value",
    "pValue": "P-Value",
    "Barcodes": "Tags",
}

def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={k: v for k, v in COLUMN_RENAME.items() if k in df.columns})
    return df

def prepare_mpra(df: pd.DataFrame, release: str, min_tags: int, promoters_only: bool) -> pd.DataFrame:
    df = normalise_columns(df)

    if "Release" in df.columns:
        df = df[df["Release"] == release].copy()
    else:
        df["Release"] = release

    if promoters_only and "Element" in df.columns:
        df = df[df["Element"].isin(PROMOTER_ELEMENTS)].copy()

    if "Tags" in df.columns:
        df = df[df["Tags"] >= min_tags].copy()

    required = ["Chromosome", "Position", "Ref", "Alt", "Tags", "DNA", "RNA", "Value", "P-Value", "Element"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing expected MPRA columns: {missing}. Present columns: {list(df.columns)}")

    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df["Tags"] = pd.to_numeric(df["Tags"], errors="coerce")
    df["DNA"] = pd.to_numeric(df["DNA"], errors="coerce")
    df["RNA"] = pd.to_numeric(df["RNA"], errors="coerce")
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
    df["P-Value"] = pd.to_numeric(df["P-Value"], errors="coerce")

    df = df.dropna(subset=["Position", "Tags", "DNA", "RNA", "Value", "P-Value", "Element"])

    df["is_deletion"] = (df["Alt"] == "-").astype(int)
    df["abs_value"] = df["Value"].abs()
    df["log_tags"] = np.log1p(df["Tags"])
    df["log_dna"] = np.log1p(df["DNA"])
    df["log_rna"] = np.log1p(df["RNA"])
    df["rna_dna_ratio"] = (df["RNA"] + 1) / (df["DNA"] + 1)

    # Relative position within each element. This is a proxy for local regulatory position.
    grouped = df.groupby("Element")["Position"]
    element_min = grouped.transform("min")
    element_max = grouped.transform("max")
    span = (element_max - element_min).replace(0, np.nan)
    df["position_scaled_within_element"] = ((df["Position"] - element_min) / span).fillna(0.0)

    keep = [
        "Chromosome", "Position", "Ref", "Alt", "Tags", "DNA", "RNA", "Value", "P-Value",
        "Element", "Release", "is_deletion", "abs_value", "log_tags", "log_dna", "log_rna",
        "rna_dna_ratio", "position_scaled_within_element"
    ]
    return df[keep].sort_values(["Element", "Position", "Ref", "Alt"])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default="data/real/kircher_promoter_mpra_variants.csv")
    parser.add_argument("--release", default="GRCh38", choices=["GRCh37", "GRCh38"])
    parser.add_argument("--min-tags", type=int, default=10)
    parser.add_argument("--promoters-only", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", compression="infer")
    prepared = prepare_mpra(df, release=args.release, min_tags=args.min_tags, promoters_only=args.promoters_only)
    prepared.to_csv(output_path, index=False)

    print(f"Wrote {len(prepared):,} rows to {output_path}")
    print(f"Elements: {prepared['Element'].nunique()}")

if __name__ == "__main__":
    main()
