# scCRISPR-seq Analysis Recipe applied to ECCITE-seq dataset of 25 knockouts generated from stimulated THP-1 cell line.
# from Papalexi et al. 2021 Nature Genetics (https://www.nature.com/articles/s41588-021-00778-2)

### Papalexi2021scCRISPR - Get & Prepare Data (custom rule) ####
rule Papalexi2021scCRISPR_get_data:
    output:
        output_dir = directory(os.path.join("data/Papalexi2021scCRISPR")),
        matrix = os.path.join("data/Papalexi2021scCRISPR/matrix.mtx"),
        barcodes = os.path.join("data/Papalexi2021scCRISPR/barcodes.tsv"),
        features = os.path.join("data/Papalexi2021scCRISPR/features.tsv"),
        metadata = os.path.join("data/Papalexi2021scCRISPR/metadata.csv"),
        metadata_selected = os.path.join("data/Papalexi2021scCRISPR/metadata_selected.csv"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","Papalexi2021scCRISPR_get_data.log"),
    script:
        "../scripts/Papalexi2021scCRISPR_get_data.R"

#### Papalexi2021scCRISPR - scRNA-seq Processing ####
module Papalexi2021scCRISPR_scrnaseq_processing_seurat:
    snakefile:
        github("epigen/scrnaseq_processing_seurat", path="workflow/Snakefile", tag="v3.0.5")
    config:
        config_wf["Papalexi2021scCRISPR_scrnaseq_processing_seurat"]

use rule * from Papalexi2021scCRISPR_scrnaseq_processing_seurat as Papalexi2021scCRISPR_scrnaseq_processing_seurat_*

#### Papalexi2021scCRISPR - Mixscape Perturbation Analysis #### 
module Papalexi2021scCRISPR_mixscape_seurat:
    snakefile:
        github("epigen/mixscape_seurat", path="workflow/Snakefile", tag="v2.0.2")
    config:
        config_wf["Papalexi2021scCRISPR_mixscape_seurat"]

use rule * from Papalexi2021scCRISPR_mixscape_seurat as Papalexi2021scCRISPR_mixscape_seurat_*

#### Papalexi2021scCRISPR - Unsupervised Analysis #### 
module Papalexi2021scCRISPR_unsupervised_analysis:
    snakefile:
        github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.3")
    config:
        config_wf["Papalexi2021scCRISPR_unsupervised_analysis"]

use rule * from Papalexi2021scCRISPR_unsupervised_analysis as Papalexi2021scCRISPR_unsupervised_analysis_*

#### Papalexi2021scCRISPR - Differential Expression Analysis #### 
module Papalexi2021scCRISPR_dea_seurat:
    snakefile:
        github("epigen/dea_seurat", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["Papalexi2021scCRISPR_dea_seurat"]

use rule * from Papalexi2021scCRISPR_dea_seurat as Papalexi2021scCRISPR_dea_seurat_*

#### Papalexi2021scCRISPR - Enrichment Analysis #### 
module Papalexi2021scCRISPR_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["Papalexi2021scCRISPR_enrichment_analysis"]

use rule * from Papalexi2021scCRISPR_enrichment_analysis as Papalexi2021scCRISPR_enrichment_analysis_*

### Papalexi2021scCRISPR - Functional Similarity Graph (custom rule) ####
rule Papalexi2021scCRISPR_functional_similarity:
    input:
        data = os.path.join("results/Papalexi2021scCRISPR/mixscape_seurat/merged/FILTERED_PRTB_data.csv"),
        metadata = os.path.join("results/Papalexi2021scCRISPR/mixscape_seurat/merged/FILTERED_metadata.csv"),
        feature_annotation = [],
    output:
        adjacency_matrix = os.path.join("results/Papalexi2021scCRISPR/special_analyses/crossprediction/adjacency_matrix.csv"),
        top_features = os.path.join("results/Papalexi2021scCRISPR/special_analyses/crossprediction/top_features.csv"),
        graph = os.path.join("results/Papalexi2021scCRISPR/special_analyses/crossprediction/graph.png"),
    params:
        group_var = "gene",
        group_rm = "NT",
        top_features_n = 5,
        prune_th = 0.2,
        feature_annotation_var = "",
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 8
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","Papalexi2021scCRISPR_functional_similarity.log"),
    script:
        "../scripts/crossprediction.py"
