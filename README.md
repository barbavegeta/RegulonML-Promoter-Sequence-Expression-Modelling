# RegulonML: Real MPRA Regulatory Variant-Effect Modelling

**RegulonML** is a reproducible bioinformatics and machine-learning workflow for modelling regulatory activity from MPRA data.

This version upgrades the original synthetic promoter-expression demo into a **real-data MPRA workflow** using the public Kircher et al. saturation-mutagenesis MPRA dataset. The project now focuses on predicting measured MPRA variant effects from regulatory-variant features and experimental count metadata.

> This is no longer just a toy promoter-expression project. It is now a real MPRA variant-effect modelling framework.

---

## What changed in this version

The original version used synthetic promoter-like sequences and simulated expression values. That was useful for demonstrating a workflow, but it was not biological evidence.

This upgraded version adds a real public MPRA data path:

- downloads Kircher et al. saturation-mutagenesis MPRA data
- filters promoter elements
- prepares a model-ready MPRA variant-effect table
- trains models to predict measured variant effects
- uses element-level group splitting to reduce leakage
- keeps the synthetic dataset only as a fallback/demo

---

## Biological context

Massively parallel reporter assays, or MPRAs, measure the regulatory activity of large numbers of DNA sequences or sequence variants in parallel. The Kircher et al. dataset tested saturation-mutagenesis variants across disease-associated promoters and enhancers and reported variant-level regulatory effects.

This repository uses the promoter subset to ask:

> Can measured MPRA regulatory effects be modelled from variant-level features such as reference allele, alternate allele, position, DNA/RNA counts, barcode support, and regulatory element identity?

This is not the same as full promoter sequence design. It is a more realistic intermediate project: **regulatory variant-effect prediction from real MPRA data**.

---

## Dataset source

The real-data workflow uses:

**Kircher M, Xiong C, Martin B, Schubach M, Inoue F, Bell RJA, Costello JF, Shendure J, Ahituv N.**  
*Saturation mutagenesis of twenty disease-associated regulatory elements at single base-pair resolution.*  
Nature Communications, 2019.

The MPRA data access portal provides variant files with columns including chromosome, position, reference allele, alternate allele, barcode/tag support, DNA count, RNA count, fitted variant-effect value, p-value, and element name.

The raw data are downloaded from the public GitHub repository backing the MPRA data portal:

```text
https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/master/data/elements.tsv.gz
```

---

## Repository structure

```text
RegulonML_real_MPRA/
├── data/
│   ├── raw/
│   │   └── downloaded raw MPRA file goes here
│   ├── real/
│   │   └── prepared real MPRA tables go here
│   └── demo/
│       └── small fallback demo file
├── docs/
│   ├── data_source_kircher_mpra.md
│   └── real_data_workflow.md
├── src/
│   ├── download_kircher_mpra.py
│   ├── prepare_kircher_mpra.py
│   ├── mpra_feature_engineering.py
│   ├── train_mpra_variant_model.py
│   └── make_small_fallback_demo.py
├── tests/
│   └── test_mpra_feature_engineering.py
├── results/
│   └── example_real_mpra/
├── run_real_mpra.sh
├── requirements.txt
├── LICENSE
├── .gitignore
└── README.md
```

---

## Quick start

Create an environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run the real-data workflow:

```bash
bash run_real_mpra.sh
```

This will:

1. Download the Kircher MPRA raw table.
2. Prepare a promoter-only MPRA variant-effect table.
3. Train baseline, Elastic Net, and Random Forest models.
4. Evaluate performance using element-level group splitting.
5. Export metrics, predictions, and feature-importance files.

---

## Manual commands

### 1. Download real MPRA data

```bash
python src/download_kircher_mpra.py \
  --output data/raw/kircher_elements.tsv.gz
```

### 2. Prepare promoter-only MPRA table

```bash
python src/prepare_kircher_mpra.py \
  --input data/raw/kircher_elements.tsv.gz \
  --output data/real/kircher_promoter_mpra_variants.csv \
  --release GRCh38 \
  --min-tags 10 \
  --promoters-only
```

### 3. Train the model

```bash
python src/train_mpra_variant_model.py \
  --input data/real/kircher_promoter_mpra_variants.csv \
  --output results/real_mpra \
  --target Value \
  --group-col Element
```

---

## Input columns expected after preparation

The prepared real MPRA table contains:

| Column | Meaning |
|---|---|
| `Chromosome` | Chromosome |
| `Position` | Variant genomic position |
| `Ref` | Reference allele |
| `Alt` | Alternate allele or deletion |
| `Tags` | Number of unique tags/barcodes |
| `DNA` | DNA count used in fitting |
| `RNA` | RNA count used in fitting |
| `Value` | Fitted log2 MPRA variant-effect coefficient |
| `P-Value` | P-value for fitted effect |
| `Element` | Regulatory element/promoter name |
| `Release` | Genome release |
| `is_deletion` | Whether the variant is a deletion |
| `position_scaled_within_element` | Relative position within element |
| `abs_value` | Absolute MPRA effect size |

---

## Feature engineering

The workflow uses real MPRA-derived features:

- reference allele one-hot encoding
- alternate allele one-hot encoding
- deletion indicator
- relative position within regulatory element
- barcode/tag support
- DNA count
- RNA count
- count ratios/log-counts
- element-level group labels

The target is the fitted MPRA effect value, `Value`.

---

## Leakage-aware validation

A random split can put variants from the same regulatory element into both train and test sets. That can inflate performance because neighbouring variants from the same promoter share local context and assay behaviour.

This project therefore supports **group splitting by element**:

```bash
--group-col Element
```

This tests whether the model generalises better across regulatory elements, not merely across nearby variants from the same element.

---

## Outputs

Training exports:

| File | Description |
|---|---|
| `metrics.csv` | Model performance comparison |
| `predictions.csv` | Observed vs predicted MPRA effect values |
| `feature_importance.csv` | Random Forest feature importance |
| `permutation_importance.csv` | Model-agnostic permutation importance |
| `run_summary.json` | Run metadata |

---

## Why this is better than the synthetic version

The old synthetic dataset proved that the code could run.

This version proves something more useful:

- you can find a real public regulatory-genomics dataset
- document its provenance
- prepare raw MPRA outputs into a modelling table
- account for leakage
- train and benchmark interpretable models
- avoid overstating the biological meaning

That is much stronger for a bioinformatics application.

---

## Still not a full promoter-design model

Be careful. This is **not** yet a full promoter sequence-expression model.

It does not yet use:

- full promoter sequences as model input
- MPRA construct-level sequence activity labels
- one-hot encoded whole sequences
- CNNs trained on full regulatory sequences
- transformer embeddings from DNA language models
- external validation across independent MPRA libraries

Those would be the next stage.

The honest project claim is:

> Real MPRA regulatory variant-effect modelling using public saturation-mutagenesis data.

Not:

> I built a promoter-design AI model.

---

## CV bullet

Built **RegulonML**, a real-data MPRA variant-effect modelling workflow using public saturation-mutagenesis reporter-assay data, including raw-data download, promoter filtering, MPRA feature engineering, leakage-aware element-level splitting, Elastic Net and Random Forest regression, permutation importance, and reproducible documentation of data provenance and limitations.

---

## GitHub description

Real-data MPRA regulatory variant-effect modelling workflow using public saturation-mutagenesis reporter-assay data, with promoter filtering, feature engineering, leakage-aware validation, Random Forest/Elastic Net benchmarking, and interpretable outputs.

---

## Suggested GitHub topics

```text
bioinformatics
mpra
regulatory-genomics
promoter
variant-effect-prediction
machine-learning
synthetic-biology
genomics
python
random-forest
elastic-net
```

---

## Author

Salvatore Barbagallo  
Bioinformatics Analyst | NGS, Computational Genomics & Data Pipelines
