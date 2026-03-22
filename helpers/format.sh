#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mrbiomics-formatting

snakefmt .
black .
npx prettier . --write
Rscript -e 'styler::style_dir("workflow/scripts")'
