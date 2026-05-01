import sys
from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from mpra_feature_engineering import encode_allele_features

def test_encode_allele_features_basic_columns():
    df = pd.DataFrame({
        "Ref": ["A", "C"],
        "Alt": ["G", "-"],
        "Tags": [10, 20],
        "DNA": [100, 200],
        "RNA": [120, 150],
        "is_deletion": [0, 1],
        "log_tags": [2.3, 3.0],
        "log_dna": [4.6, 5.3],
        "log_rna": [4.8, 5.0],
        "rna_dna_ratio": [1.2, 0.75],
        "position_scaled_within_element": [0.1, 0.9],
    })
    X = encode_allele_features(df)
    assert "ref_A" in X.columns
    assert "alt_-" in X.columns
    assert X.loc[0, "is_transition"] == 1
    assert X.loc[1, "is_deletion"] == 1
