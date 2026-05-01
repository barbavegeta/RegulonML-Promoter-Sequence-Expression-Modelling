#!/usr/bin/env bash
set -euo pipefail

python src/download_kircher_mpra.py \
  --output data/raw/kircher_elements.tsv.gz

python src/prepare_kircher_mpra.py \
  --input data/raw/kircher_elements.tsv.gz \
  --output data/real/kircher_promoter_mpra_variants.csv \
  --release GRCh38 \
  --min-tags 10 \
  --promoters-only

python src/train_mpra_variant_model.py \
  --input data/real/kircher_promoter_mpra_variants.csv \
  --output results/real_mpra \
  --target Value \
  --group-col Element
