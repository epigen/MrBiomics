# Quick Start

This is the fastest end-to-end example in `MrBiomics`.

It runs the `enrichment_analysis` module on three small "mystery" genomic region sets and produces a single summary plot that is easy to interpret. The feature sets come from Corces ATAC-seq cell-type signatures, but they are intentionally relabeled as `setA`, `setB`, and `setC` so the quick start feels like a small real analysis instead of a canned demo.

## What it does

- Input: three BED files with genomic regions plus one shared background BED file
- Module: `enrichment_analysis`
- Database: `Azimuth_2023`
- Expected runtime: a few seconds to a few minutes, depending on your environment
- Main output: a clustered bubble summary plot for the three mystery sets

## Files used by the quick start

- Annotation table: `config/quickstart/quickstart_enrichment_analysis_annotation.csv`
- Module config: `config/quickstart/quickstart_enrichment_analysis_config.yaml`
- Input BED files: `docs/quickstart/data/`
- Database: `resources/Azimuth_2023.gmt`

## Run it

From the repository root:

```bash
snakemake --software-deployment-method conda --cores 1
```

This uses the active default target in `workflow/Snakefile`, which currently points to the quick-start enrichment analysis.

## Where to look after it finishes

The quick-start results are written to:

```text
results/quickstart/
```

The most useful output to inspect first is the grouped Azimuth summary produced for `mystery_cell_types`, which should clearly separate the three sets into recognizable hematopoietic identities.

## Mystery set key

If you want to keep the exercise spoiler-free, stop reading here.

- `setA`: monocyte-derived regions
- `setB`: B-cell-derived regions
- `setC`: erythroid-derived regions
