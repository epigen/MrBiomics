# 🚀 MR.PARETO - Modules & Recipes for Pragmatic Augmentation of Research Efficiency Towards Optimum
> _"For many outcomes, roughly 80% of consequences come from 20% of causes (the "vital few")."_ - The Pareto Principle by Vilfredo Pareto

Get 80% of all standard (biomedical) data science analyses done semi-automated with 20% of the effort, by leveraging [Snakemake's](https://snakemake.github.io/) module functionality to use and combine pre-existing workflows into arbitrarily complex analyses.

# ⏳ TL;DR - More Time for Science!
> _"Programming is about trying to make the future less painful. It’s about making things easier for our teammates."_ from _The Pragmatic Programmer_ by Andy Hunt & Dave Thomas
- **Why**: Time is the most precious resource. By taking care of efficiency (i.e., maximum output with limited resources) scientists can re-distribute their time to focus on effectiveness (i.e., the biggest impact possible).
- **How**: Leverage the latest developments in workflow management to (re-)use and combine independent computational modules into arbitrarily complex analyses in combination with modern innovation methods (e.g., fast prototyping, design thinking, and agile concepts).
- **What**: Independent computational **Modules** implemented as Snakemake workflows, encoding best practices and standard approaches, are used to scale, automate, and parallelize analyses. Snakemake's module functionality enables arbitrarily complex combinations of pre-existing modules for any **Project**. **Recipes** combine modules into the most conceivable standard analyses, thereby accelerating projects to the point of the unknown. 

Altogether this enables complex, portable, transparent, reproducible, and documented analyses of biomedical data at scale.

# 🧩 Modules
> _"Is it functional, multifunctional, durable, well-fitted, simple, easy to maintain, and thoroughly tested? Does it provide added value, and doesn't cause unnecessary harm? Can it be simpler? Is it an innovation?"_ - Patagonia Design Principles

**Modules** are Snakemake workflows, consisting of **Rules** for multi-step analyses, that are independent and self-contained. A `{module}` can be general-purpose (e.g., Unsupervised Analysis) or modality-specific (e.g., ATAC-seq processing). Currently, the following **nine modules** are available, ordered by their applicability from general to specific:

| Module | Type (Data Modality) | DOI | Stars |
| :---: | :---: | :---: | :---: |
| [Unsupervised Analysis](https://github.com/epigen/unsupervised_analysis) | General Purpose<br>(tabular data) | [![DOI](https://zenodo.org/badge/475465311.svg)](https://zenodo.org/badge/latestdoi/475465311) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/unsupervised_analysis?style=plastic"> |
| [Split, Filter, Normalize and Integrate Sequencing Data](https://github.com/epigen/spilterlize_integrate) | Bioinformatics<br>(NGS counts) | [![DOI](https://zenodo.org/badge/659800258.svg)](https://zenodo.org/badge/latestdoi/659800258) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/spilterlize_integrate?style=plastic"> |
| [Differential Analysis with limma](https://github.com/epigen/dea_limma) | Bioinformatics<br>(NGS data) | [![DOI](https://zenodo.org/badge/524038188.svg)](https://zenodo.org/badge/latestdoi/524038188) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/dea_limma?style=plastic">
| [Enrichment Analysis](https://github.com/epigen/enrichment_analysis) | Bioinformatics<br>(genes/genomic regions) | [![DOI](https://zenodo.org/badge/377527671.svg)](https://zenodo.org/badge/latestdoi/377527671) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/enrichment_analysis?style=plastic"> |
| [Genome Track Visualization](https://github.com/epigen/genome_tracks) | Bioinformatics<br>(aligned BAM files) | [![DOI](https://zenodo.org/badge/438573546.svg)](https://zenodo.org/doi/10.5281/zenodo.10849097) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/genome_tracks?style=plastic"> |
| [ATAC-seq Processing](https://github.com/epigen/atacseq_pipeline) | Bioinformatics<br>(ATAC-seq) | [![DOI](https://zenodo.org/badge/350342694.svg)](https://zenodo.org/doi/10.5281/zenodo.6323634) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/atacseq_pipeline?style=plastic"> |
| [scRNA-seq Processing using Seurat](https://github.com/epigen/scrnaseq_processing_seurat) | Bioinformatics<br>(scRNA-seq) | [![DOI](https://zenodo.org/badge/468752490.svg)](https://zenodo.org/doi/10.5281/zenodo.10679327) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/scrnaseq_processing_seurat?style=plastic"> |
| [Differential Analysis using Seurat](https://github.com/epigen/dea_seurat) | Bioinformatics<br>(scRNA-seq) | [![DOI](https://zenodo.org/badge/483638364.svg)](https://zenodo.org/doi/10.5281/zenodo.10689139) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/dea_seurat?style=plastic"> |
| [Perturbation Analysis using Mixscape from Seurat](https://github.com/epigen/mixscape_seurat) | Bioinformatics<br>(scCRISPR-seq) | [![DOI](https://zenodo.org/badge/481635018.svg)](https://zenodo.org/badge/latestdoi/481635018) | <img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/epigen/mixscape_seurat?style=plastic"> |

> [!NOTE]  
> ⭐️ **Star and share modules you find valuable** 📤 — help others discover them, and guide our focus for future work!

> [!TIP]
> For detailed instructions on the [installation](../../wiki/Installation), [configuration](../../wiki/Configuration), and [execution](../../wiki/Execution) of modules, you can check out the [wiki](https://github.com/epigen/mr.pareto/wiki). Generic instructions are also shown in the modules' respective [Snakmake workflow catalog entry](https://snakemake.github.io/snakemake-workflow-catalog).

# 📋 Projects using multiple Modules
> _“Absorb what is useful. Discard what is not. Add what is uniquely your own.”_ - Bruce Lee

You can [(re-)use and combine pre-existing workflows](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows) within your projects by loading them as [**Modules**](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules) since [Snakemake 6](https://slides.com/johanneskoester/snakemake-6#/8). The combination of multiple modules into projects that analyze multiple datasets represents the overarching vision and power of MR.PARETO.

> [!NOTE]
> When applied to multiple datasets within a project, each dataset should have its own result directory within the project directory.

Three components are required to use a module within your Snakemake workflow (i.e., a project).
- Configuration: The [`config/config.yaml`](./config/config.yaml) file has to point to the respective configuration files per dataset and workflow.
  ```python
  #### Datasets and Workflows to include ###
  workflows:
      MyData:
          other_workflow: "config/MyData/MyData_other_workflow_config.yaml"
  ```
- Snakefile: Within the main Snakefile ([`workflow/Snakefile`](.workflow/Snakefile)) we have to:
  - load all configurations;
  - include the snakefiles that contain the dataset-specific loaded modules and rules (see next point);
  - and add all modules' outputs to the target's rule `input`.
- Modules: Load the required modules and their rules within separate snakefiles (`*.smk`) in the `rule/` folder. Recommendation: Use one snakefile per dataset.
  ```python
    module MyData_other_workflow:
      # here, plain paths, URLs and the special markers for code hosting providers (e.g., github) are possible.
      snakefile: "other_workflow/Snakefile"
      config: config["MyData_other_workflow"]

    use rule * from MyData_other_workflow as MyData_other_workflow_*
    ```

> [!TIP]
> A full [tutorial](../../wiki/Module-Usage-in-Projects) is available on the [wiki](https://github.com/epigen/mr.pareto/wiki).

# 📜 Recipes
> _"Civilization advances by extending the number of important operations which we can perform without thinking of them."_ - Alfred North Whitehead, author of _Principia Mathematica_

**Recipes** are combinations of existing modules into end-to-end best practice analyses. They can be used as templates for standard analyses by leveraging existing modules, thereby enabling fast iterations and progression into the unknown. Every recipe is described and presented using a [wiki](https://github.com/epigen/mr.pareto/wiki) page by application to a public data set.

> [!TIP]
> Process each dataset module by module. Check the results of each module to inform the configuration of the next module. This iterative method allows for quick initial completion, followed by refinement in subsequent iterations based on feedback from yourself or collaborators. Adjustments in later iterations are straightforward, requiring only changes to individual configurations or annotations. Ultimately you end up with a reproducible and readable end-to-end analysis for each dataset.

| Recipe | Description | # Modules | Results |
| :---: | :---: | :---: | :---: |
| [ATAC-seq Analysis](../../wiki/ATACseq-Analysis-Recipe) | From raw BAM files to enrichemnts of differentially accessible regions. | 6(-7) | ... |
| [RNA-seq Analysis](../../wiki/RNAseq-Analysis-Recipe) | From raw BAM files to enrichemnts of differentially expressed genes. | 6(-7) | ... |
| [Integrative ATAC-seq & RNA-seq Analysis](../../wiki/Integrative-Analysis-Recipe) | From count matrices to epigenetic potential and relative transcriptional abundance. | 7(-8) | ... |
| [scRNA-seq Analysis](../../wiki/scRNAseq-Analysis-Recipe) | From count matrix to enrichemnts of differentially expressed genes. | 5(-6) | ... |
| [scCRISPR-seq Analysis](../../wiki/scCRISPRseq-Analysis-Recipe) | From count matrix to knockout phenotype enrichemnts. | 6(-7) | ... |

> [!NOTE]  
> ⭐️ **Star this repository and share recipes you find valuable** 📤 — help others discover them, and guide our focus for future work!

# 📚 Resources
- [MR.PARETO Wiki for instructions & tutorials](https://github.com/epigen/mr.pareto/wiki)
- [GitHub list of MR.PARETO modules](https://github.com/stars/sreichl/lists/mr-pareto)
- [My Data Science Setup - A tutorial for developing Snakemake workflows and beyond](https://bit.ly/TAP-data-science-setup)
- [GitHub Page of this repository](https://epigen.github.io/mr.pareto/)
- Curated and published workflows that could be used as modules:
    - [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
    - [WorkflowHub](https://workflowhub.eu/)
    - [SnakePipes](https://snakepipes.readthedocs.io/en/latest/index.html)
    - [Seq2Science](https://vanheeringen-lab.github.io/seq2science/index.html)
- Software
    - [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
    - [Conda](https://docs.conda.io/en/latest/)
    - [Docker](https://www.docker.com/)
    - [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/index.html#) 

# ⭐ Star History of Modules
[![Star History Chart](https://api.star-history.com/svg?repos=epigen/atacseq_pipeline,epigen/unsupervised_analysis,epigen/enrichment_analysis,epigen/spilterlize_integrate,epigen/dea_limma,epigen/dea_seurat,epigen/scrnaseq_processing_seurat,epigen/mixscape_seurat,epigen/genome_tracks&type=Date)](https://star-history.com/#epigen/atacseq_pipeline&epigen/unsupervised_analysis&epigen/enrichment_analysis&epigen/spilterlize_integrate&epigen/dea_limma&epigen/dea_seurat&epigen/scrnaseq_processing_seurat&epigen/mixscape_seurat&epigen/genome_tracks&Date)
