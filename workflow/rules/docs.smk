# Snakemake file for creating plots for the recipes and documentation in the wiki


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
# This rule can only be used for docs after all results (including untracked ones: volcano and unsupervised analysis plots)
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

#### CorcesATAC - Copy selected plots for documentation and visualization in wiki (custom rule) ####
# Define the mapping of input to output files
CorcesATAC_plots_map = {
    "sample_annotation.png": "results/CorcesATAC/atacseq_pipeline/report/sample_annotation.png",
    "CD34.svg": "results/CorcesATAC/genome_tracks/tracks/CD34.svg",
    "MS4A1.svg": "results/CorcesATAC/genome_tracks/tracks/MS4A1.svg",
    "filtered.png": "results/CorcesATAC/spilterlize_integrate/all/plots/filtered.png",
    "normCQN_integrated.png": "results/CorcesATAC/spilterlize_integrate/all/plots/normCQN_integrated.png",
    "filtered_CFA.png": "results/CorcesATAC/spilterlize_integrate/all/plots/filtered_CFA.png",
    "normCQN_integrated_CFA.png": "results/CorcesATAC/spilterlize_integrate/all/plots/normCQN_integrated_CFA.png",
    "normCQN_integrated_PCA.png": "results/CorcesATAC/unsupervised_analysis/normCQN_integrated/PCA/plots/PCA_auto_0.9_2/metadata/cell_type.png",
    "normCQN_integrated_HVF_PCA.png": "results/CorcesATAC/unsupervised_analysis/normCQN_integrated_HVF/PCA/plots/PCA_auto_0.9_2/metadata/cell_type.png",
    "normCQN_integrated_UMAP.png": "results/CorcesATAC/unsupervised_analysis/normCQN_integrated/UMAP/plots/UMAP_correlation_15_0.1_2/metadata/cell_type.png",
    "normCQN_integrated_HVF_UMAP.png": "results/CorcesATAC/unsupervised_analysis/normCQN_integrated_HVF/UMAP/plots/UMAP_correlation_15_0.1_2/metadata/cell_type.png",
    "dea_stats.png": "results/CorcesATAC/dea_limma/normCQN_OvA_cell_type/plots/stats.png",
    "markerGenes.png": "results/CorcesATAC/dea_limma/normCQN_OvA_cell_type/plots/heatmap/markerGenes.png",
    "Bcell_up_ReactomePathways.png": "results/CorcesATAC/enrichment_analysis/Bcell_up/GREAT/ReactomePathways/Bcell_up_ReactomePathways.png",
    "Bcell_up_Azimuth_2023.png": "results/CorcesATAC/enrichment_analysis/Bcell_up/GREAT/Azimuth_2023/Bcell_up_Azimuth_2023.png",
    "cell_types_Azimuth_2023_summary.png": "results/CorcesATAC/enrichment_analysis/cell_types_up/GREAT/Azimuth_2023/cell_types_up_Azimuth_2023_summary.png",
    "cell_types_ReactomePathways_summary.png": "results/CorcesATAC/enrichment_analysis/cell_types_up/GREAT/ReactomePathways/cell_types_up_ReactomePathways_summary.png",
    "crossprediction_graph.png": "results/CorcesATAC/special_analyses/crossprediction/graph.png",
}

# Copy input to outputs to include the plots in the repo and wiki
# This rule can only be used for docs after all results (including untracked ones: volcano and unsupervised analysis plots)
# have been generated (i.e. leave commented in the Snakefile's target rule until the end/last iteration)
rule CorcesATAC_plots:
    input:
        [CorcesATAC_plots_map[plot] for plot in CorcesATAC_plots_map]
    output:
        [f"docs/CorcesATAC/{plot}" for plot in CorcesATAC_plots_map]
    resources:
        mem_mb="1000",
    threads: config.get("threads", 1)
    log:
        "logs/rules/CorcesATAC_plots.log",
    run:
        for i, o in zip(input, output):
            shell(f"cp {i} {o}")

#### CorcesINT - Copy selected plots for documentation and visualization in wiki (custom rule) ####
# Define the mapping of input to output files
CorcesINT_plots_map = {
    "filtered.png": "results/CorcesINT/spilterlize_integrate/all/plots/filtered.png",
    "normupperquartile_integrated.png": "results/CorcesINT/spilterlize_integrate/all/plots/normupperquartile_integrated.png",
    "filtered_CFA.png": "results/CorcesINT/spilterlize_integrate/all/plots/filtered_CFA.png",
    "normupperquartile_integrated_CFA.png": "results/CorcesINT/spilterlize_integrate/all/plots/normupperquartile_integrated_CFA.png",
    "normupperquartile_integrated_PCA_cell_type.png": "results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/PCA/plots/PCA_auto_0.9_2/metadata/cell_type.png",
    "normupperquartile_integrated_UMAP_cell_type.png": "results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/UMAP/plots/UMAP_correlation_15_0.1_2/metadata/cell_type.png",
    "normupperquartile_integrated_PCA_modality.png": "results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/PCA/plots/PCA_auto_0.9_2/metadata/modality.png",
    "normupperquartile_integrated_UMAP_modality.png": "results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/UMAP/plots/UMAP_correlation_15_0.1_2/metadata/modality.png",
    "dea_stats.png": "results/CorcesINT/dea_limma/normupperquartile_integrated/plots/stats.png",
    "markerGenes.png": "results/CorcesINT/dea_limma/normupperquartile_integrated/plots/heatmap/markerGenes.png",
    "HSC_correlation.png": "results/CorcesINT/special_analysis/correlation_plots/HSC_correlation.png",
    "Mono_correlation.png": "results/CorcesINT/special_analysis/correlation_plots/Mono_correlation.png",
    "cell_types_GO_Biological_Process_2025_summary.png": "results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/GO_Biological_Process_2025/cell_types_GO_Biological_Process_2025_summary.png",
    "cell_types_ReactomePathways_summary.png": "results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/ReactomePathways/cell_types_ReactomePathways_summary.png",
    "Mono_EP.png": "results/CorcesINT/enrichment_analysis/Mono_EP/RcisTarget/hg38_500bp_up_100bp_down_v10clust/Mono_EP_hg38_500bp_up_100bp_down_v10clust.png",
    "Mono_TA.png": "results/CorcesINT/enrichment_analysis/Mono_TA/RcisTarget/hg38_500bp_up_100bp_down_v10clust/Mono_TA_hg38_500bp_up_100bp_down_v10clust.png",
}

# Copy input to outputs to include the plots in the repo and wiki
# This rule can only be used for docs after all results (including untracked ones: volcano and unsuervised analysis plots)
# have been generated (i.e. leave commented in the Snakefile's target rule until the end/last iteration)
rule CorcesINT_plots:
    input:
        [CorcesINT_plots_map[plot] for plot in CorcesINT_plots_map]
    output:
        [f"docs/CorcesINT/{plot}" for plot in CorcesINT_plots_map]
    resources:
        mem_mb="1000",
    threads: config.get("threads", 1)
    log:
        "logs/rules/CorcesINT_plots.log",
    run:
        for i, o in zip(input, output):
            shell(f"cp {i} {o}")
