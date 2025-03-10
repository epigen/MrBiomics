[![MrBiomics](https://img.shields.io/badge/MrBiomics-red)](https://github.com/epigen/MrBiomics/)
[![DOI](https://zenodo.org/badge/XXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXX) 
[![](https://tokei.rs/b1/github/epigen/XXXX?category=code)]() 
[![](https://tokei.rs/b1/github/epigen/XXXX?category=files)]()
[![GitHub license](https://img.shields.io/github/license/epigen/XXXX)](https://github.com/epigen/XXXX/blob/master/LICENSE)
![GitHub Release](https://img.shields.io/github/v/release/epigen/XXXXX)
[![Snakemake](https://img.shields.io/badge/Snakemake->=8.20.1-green)](https://snakemake.readthedocs.io/en/stable/)

# Analysis & Visualization Workflow Using/Powered by PACKAGE for DATA/MODALITY
A [Snakemake 8](https://snakemake.readthedocs.io/en/stable/) workflow for performing and visualizing analyses of data (e.g., ...) powered by the package [package](https://www.packageURL.org).

> [!NOTE]  
> This workflow adheres to the module specifications of [MrBiomics](https://github.com/epigen/MrBiomics), an effort to augment research by modularizing (biomedical) data science. For more details, instructions, and modules check out the project's repository.
>
> ⭐️ **Star and share modules you find valuable** 📤 - help others discover them, and guide our future work!

> [!IMPORTANT]  
> **If you use this workflow in a publication, please don't forget to give credit to the authors by citing it using this DOI [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

# 🖋️ Authors
- [MrBiomics](https://github.com/epigen/MrBiomics)


# 💿 Software
This project wouldn't be possible without the following software and their dependencies.

| Software | Reference (DOI) |
| :---: | :---: |
| Snakemake | https://doi.org/10.12688/f1000research.29032.2 |
| packageA | https://doi.org/10.AAA/ |
| packageB | https://doi.org/10.BBBB/ |


# 🔬 Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (`workflow/envs/*.yaml file`) or post-execution in the result directory (`{module}/envs/*.yaml`). Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g., [X].

__Analysis.__ Analysis was performed...

__Visualization.__ The results were visualized...


**The analysis and visualizations described here were performed using a publicly available Snakemake (ver) [ref] workflow [ref - cite this workflow here].**

# 🚀 Features
The workflow performs the following steps that produce the outlined results:

- Analysis
  - ...
  - (optional) ...
- Visualizations
  - ...
- Limitations
  - ...


# 🛠️ Usage
Here are some tips for the usage of this workflow:
- ...

# ⚙️ Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# 📖 Examples
 Explore detailed examples showcasing module usage in comprehensive end-to-end analyses (including data, configuration, annotation and results) in our [MrBiomics Recipes](https://github.com/epigen/MrBiomics?tab=readme-ov-file#-recipes):
- [ATACseq Analysis Recipe](https://github.com/epigen/MrBiomics/wiki/ATACseq-Analysis-Recipe)
- ...

# 🔗 Links
- [GitHub Repository](https://github.com/user/module/)
- [GitHub Page](https://user.github.io/module/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=user/module)

# 📚 Resources
- Recommended compatible [MrBiomics Modules](https://github.com/epigen/MrBiomics/#-modules) for up-/downstream analyses:
  - [Unsupervised Analysis](https://github.com/epigen/unsupervised_analysis) to understand and visualize similarities and variations between cells/samples, including dimensionality reduction and cluster analysis. Useful for all tabular data including single-cell and bulk sequencing data.
  - [<ins>Sp</ins>lit, F<ins>ilter</ins>, Norma<ins>lize</ins> and <ins>Integrate</ins> Sequencing Data](https://github.com/epigen/spilterlize_integrate/) after count quantification.
  - [Differential Analysis with limma](https://github.com/epigen/dea_limma) to identify and visualize statistically significantly different features (e.g., genes or genomic regions) between sample groups.
  - [Enrichment Analysis](https://github.com/epigen/enrichment_analysis) for biomedical interpretation of (differential) analysis results using prior knowledge.
  - [Genome Browser Track Visualization](https://github.com/epigen/genome_tracks/) for quality control and visual inspection/analysis of genomic regions/genes of interest or top hits.
  - [ATAC-seq Data Processing & Quantification Pipeline](https://github.com/epigen/atacseq_pipeline) for processing, quantification and annotation of chromatin accessibility.
  - [scRNA-seq Data Processing & Visualization](https://github.com/epigen/scrnaseq_processing_seurat) for processing (multimodal) single-cell transcriptome data.
  - [Differential Analysis using Seurat](https://github.com/epigen/dea_seurat) to identify and visualize statistically significantly different features (e.g., genes or proteins) between groups.
  - [Perturbation Analysis using Mixscape from Seurat](https://github.com/epigen/mixscape_seurat) to identify perturbed cells from pooled (multimodal) CRISPR screens with sc/snRNA-seq read-out (scCRISPR-seq).

# 📑 Publications
The following publications successfully used this module for their analyses.
- [FirstAuthors et al. (202X) Journal Name - Paper Title.](https://doi.org/10.XXX/XXXX)
- ...

# ⭐ Star History

[![Star History Chart](https://api.star-history.com/svg?repos=epigen/XXX&type=Date)](https://star-history.com/#epigen/XXX&Date)
