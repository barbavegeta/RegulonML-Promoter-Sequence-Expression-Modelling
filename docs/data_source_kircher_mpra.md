# Real data source: Kircher saturation-mutagenesis MPRA

This project uses the public Kircher et al. saturation-mutagenesis MPRA dataset.

## Source

Kircher M, Xiong C, Martin B, Schubach M, Inoue F, Bell RJA, Costello JF, Shendure J, Ahituv N.  
*Saturation mutagenesis of twenty disease-associated regulatory elements at single base-pair resolution.*  
Nature Communications, 2019.

## Data portal

The public portal describes MPRA saturation mutagenesis across disease-associated promoters and enhancers. The downloadable table includes variant-level columns such as chromosome, position, reference allele, alternate allele, tags/barcodes, DNA count, RNA count, fitted value, p-value, and element.

## Why this dataset

It is real reporter-assay data and is much stronger than the original synthetic demonstration dataset. It is not, however, full promoter sequence-expression data. It is variant-level MPRA effect data.

## Modelling target

The target is the fitted MPRA variant effect (`Coefficient` in the raw file: log2 change in reporter
activity caused by the variant). The DNA, RNA and barcode counts are the measurements that
effect is estimated from, so they are used only for quality filtering (at least 10 barcodes)
and never as model features.

## Loci

The 29 datasets come from 21 genomic loci. Some loci were tested more than once (TERT in four
cell lines, PKLR at two time points, LDLR and SORT1 in replicate or flipped constructs, ZRS in two
versions). Train/test splits are made by locus so the same DNA never appears on both sides.
