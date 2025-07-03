# RNA-seq Analysis Recipe for the healthy hematopoeitic RNA-seq samples 
# from Corces et al. 2016 Nature Genetics (https://www.nature.com/articles/ng.3646)
# Subset of GEO SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75384
# Subset of GEO SuperSeries for **unstranded** RNA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246

### CorcesRNA - Fetch NGS ####
module CorcesRNA_fetch_ngs:
    snakefile:
       github("epigen/fetch_ngs", path="workflow/Snakefile", tag="v1.0.5")
    config:
        config_wf["CorcesRNA_fetch_ngs"]

use rule * from CorcesRNA_fetch_ngs as CorcesRNA_fetch_ngs_*

### CorcesRNA - RNA-seq processing ####
module CorcesRNA_rnaseq_pipeline:
    snakefile:
        github("epigen/rnaseq_pipeline", path="workflow/Snakefile", tag="v1.2.1")
    config:
        config_wf["CorcesRNA_rnaseq_pipeline"]

use rule * from CorcesRNA_rnaseq_pipeline as CorcesRNA_rnaseq_pipeline_*

#### CorcesRNA - Genome Tracks #### 
module CorcesRNA_genome_tracks:
    snakefile:
        github("epigen/genome_tracks", path="workflow/Snakefile", tag="v2.0.5")
    config:
        config_wf["CorcesRNA_genome_tracks"]

use rule * from CorcesRNA_genome_tracks as CorcesRNA_genome_tracks_*

#### CorcesRNA - Spilterlize & Integrate #### 
module CorcesRNA_spilterlize_integrate:
    snakefile:
        github("epigen/spilterlize_integrate", path="workflow/Snakefile", tag="v3.0.1")
    config:
        config_wf["CorcesRNA_spilterlize_integrate"]

use rule * from CorcesRNA_spilterlize_integrate as CorcesRNA_spilterlize_integrate_*

### CorcesRNA - Unsupervised Analysis ####
module CorcesRNA_unsupervised_analysis:
    snakefile:
        github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.3")
    config:
        config_wf["CorcesRNA_unsupervised_analysis"]

use rule * from CorcesRNA_unsupervised_analysis as CorcesRNA_unsupervised_analysis_*

#### CorcesRNA - Differential Expression Analysis #### 
module CorcesRNA_dea_limma:
    snakefile:
        github("epigen/dea_limma", path="workflow/Snakefile", branch="main")
    config:
        config_wf["CorcesRNA_dea_limma"]

use rule * from CorcesRNA_dea_limma as CorcesRNA_dea_limma_*

#### CorcesRNA - Enrichment Analysis ####
module CorcesRNA_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["CorcesRNA_enrichment_analysis"]

use rule * from CorcesRNA_enrichment_analysis as CorcesRNA_enrichment_analysis_*

#### CorcesRNA - Lineage Reconstruction (custom rule) ####
rule CorcesRNA_reconstruct_lineage:
    input:
        data = os.path.join("results/CorcesRNA/spilterlize_integrate/all/normCQN_integrated_HVF.csv"),
        metadata = os.path.join("results/CorcesRNA/spilterlize_integrate/all/annotation.csv"),
        feature_annotation = os.path.join("results/CorcesRNA/rnaseq_pipeline/counts/gene_annotation.csv"),
    output:
        adjacency_matrix = os.path.join("results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"),
        top_features = os.path.join("results/CorcesRNA/special_analyses/crossprediction/top_features.csv"),
        graph = os.path.join("results/CorcesRNA/special_analyses/crossprediction/graph.png"),
    params:
        group_var = "cell_type",
        top_features_n = 5,
        prune_th = 0.2,
        feature_annotation_var = "external_gene_name",
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","CorcesRNA_reconstruct_lineage.log"),
    script:
        "../scripts/crossprediction.py"

#### CorcesRNA - Copy selected plots for documentation and visualization in wiki (custom rule) ####
# Define the mapping of input to output files
CorcesRNA_plots_map = {
    "sample_annotation.png": "results/CorcesRNA/rnaseq_pipeline/report/sample_annotation.png",
    "CD34.svg": "results/CorcesRNA/genome_tracks/tracks/CD34.svg",
    "MS4A1.svg": "results/CorcesRNA/genome_tracks/tracks/MS4A1.svg",
    "filtered.png": "results/CorcesRNA/spilterlize_integrate/all/plots/filtered.png",
    "normCQN_integrated.png": "results/CorcesRNA/spilterlize_integrate/all/plots/normCQN_integrated.png",
    "filtered_CFA.png": "results/CorcesRNA/spilterlize_integrate/all/plots/filtered_CFA.png",
    "normCQN_integrated_CFA.png": "results/CorcesRNA/spilterlize_integrate/all/plots/normCQN_integrated_CFA.png",
    "normCQN_integrated_PCA.png": "results/CorcesRNA/unsupervised_analysis/normCQN_integrated/PCA/plots/PCA_auto_0.9_2/metadata/cell_type.png",
    "normCQN_integrated_HVF_PCA.png": "results/CorcesRNA/unsupervised_analysis/normCQN_integrated_HVF/PCA/plots/PCA_auto_0.9_2/metadata/cell_type.png",
    "normCQN_integrated_UMAP.png": "results/CorcesRNA/unsupervised_analysis/normCQN_integrated/UMAP/plots/UMAP_correlation_15_0.1_2/metadata/cell_type.png",
    "normCQN_integrated_HVF_UMAP.png": "results/CorcesRNA/unsupervised_analysis/normCQN_integrated_HVF/UMAP/plots/UMAP_correlation_15_0.1_2/metadata/cell_type.png",
    "dea_stats.png": "results/CorcesRNA/dea_limma/normCQN_OvA_cell_type/plots/stats.png",
    "markerGenes.png": "results/CorcesRNA/dea_limma/normCQN_OvA_cell_type/plots/heatmap/markerGenes.png",
    "Bcell_ReactomePathways.png": "results/CorcesRNA/enrichment_analysis/Bcell/preranked_GSEApy/ReactomePathways/Bcell_ReactomePathways.png",
    "Bcell_Azimuth_2023.png": "results/CorcesRNA/enrichment_analysis/Bcell/preranked_GSEApy/Azimuth_2023/Bcell_Azimuth_2023.png",
    "cell_types_Azimuth_2023_summary.png": "results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_summary.png",
    "cell_types_ReactomePathways_summary.png": "results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/ReactomePathways/cell_types_ReactomePathways_summary.png",
    "crossprediction_graph.png": "results/CorcesRNA/special_analyses/crossprediction/graph.png",
}

# Copy input to outputs to include the plots in the repo and wiki
# This rule can only be used for docs after all results (including untracked ones: volcano and unsuervised analysis plots)
# have been generated (i.e. leave commented in the Snakefile's target rule until the end/last iteration)
rule CorcesRNA_plots:
    input:
        [CorcesRNA_plots_map[plot] for plot in CorcesRNA_plots_map]
    output:
        [f"docs/CorcesRNA/{plot}" for plot in CorcesRNA_plots_map]
    resources:
        mem_mb="1000",
    threads: config.get("threads", 1)
    log:
        "logs/rules/CorcesRNA_plots.log",
    run:
        for i, o in zip(input, output):
            shell(f"cp {i} {o}")
