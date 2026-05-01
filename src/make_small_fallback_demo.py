from pathlib import Path
import pandas as pd

def main():
    out = Path("data/demo/fallback_small_mpra_demo.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame([
        {"Chromosome": "chr11", "Position": 5227079, "Ref": "A", "Alt": "G", "Tags": 20, "DNA": 120, "RNA": 180, "Value": 0.45, "P-Value": 0.001, "Element": "HBB", "Release": "GRCh38"},
        {"Chromosome": "chr11", "Position": 5227080, "Ref": "C", "Alt": "T", "Tags": 15, "DNA": 90, "RNA": 60, "Value": -0.58, "P-Value": 0.003, "Element": "HBB", "Release": "GRCh38"},
        {"Chromosome": "chr19", "Position": 11089240, "Ref": "G", "Alt": "A", "Tags": 32, "DNA": 300, "RNA": 450, "Value": 0.62, "P-Value": 0.0005, "Element": "LDLR", "Release": "GRCh38"},
    ])
    df.to_csv(out, index=False)
    print(f"Wrote fallback demo to {out}")

if __name__ == "__main__":
    main()
