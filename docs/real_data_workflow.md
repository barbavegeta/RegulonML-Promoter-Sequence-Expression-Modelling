# Real-data workflow

## Step 1: Download

```bash
python src/download_kircher_mpra.py --output data/raw/kircher_elements.tsv.gz
```

## Step 2: Prepare promoter subset

```bash
python src/prepare_kircher_mpra.py \
  --input data/raw/kircher_elements.tsv.gz \
  --output data/real/kircher_promoter_mpra_variants.csv \
  --release GRCh38 \
  --min-tags 10 \
  --promoters-only
```

## Step 3: Train

```bash
python src/train_mpra_variant_model.py \
  --input data/real/kircher_promoter_mpra_variants.csv \
  --output results/real_mpra \
  --target Value \
  --group-col Element
```

## Warning

This is a real MPRA project, but it is not a full promoter-design model. It predicts measured variant effects from tabular variant/count features.
