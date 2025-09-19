# Snakemake file for creating the figure for the manuscript
#### CorcesRNA - Figure panels (custom rule) ####
rule CorcesRNA_figures:
    input:
        umap_coords = os.path.join("results/CorcesRNA/unsupervised_analysis/normCQN_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"),
        dea_ova = os.path.join("results/CorcesRNA/dea_limma/normCQN_OvA_cell_type/results.csv"),
        enrichment_results = os.path.join("results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_all.csv"),
        crossprediction_adj_mtx = os.path.join("results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"),
    output:
        umap_plot = os.path.join("paper/CorcesRNA/umap.pdf"),
        dea_heatmap_plot = os.path.join("paper/CorcesRNA/differential_heatmap.pdf"),
        enrichment_plot = os.path.join("paper/CorcesRNA/enrichment.pdf"),
        crossprediction_plot = os.path.join("paper/CorcesRNA/crossprediction.pdf"),
        crossprediction_coordinates = os.path.join("paper/CorcesRNA/crossprediction_coordinates.csv"),
    params:
        figure_theme_path = workflow.source_path("../scripts/figure_theme.R"),
        figure_utils_path = workflow.source_path("../scripts/figure_utils.R"),
        # enrichment analysis
        fdr_threshold = 0.05,
        log2FC_threshold = 3,
        # lineage reconstructions
        lineage_tree_cut_off = 0.07,
        hierarchy_coordinates = True
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","CorcesRNA_figures.log"),
    script:
        "../scripts/CorcesRNA_figures.R"

#### CorcesATAC - Figure panels (custom rule) ####
rule CorcesATAC_figures:
    input:
        umap_coords = os.path.join("results/CorcesATAC/unsupervised_analysis/normCQN_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"),
        dea_ova = os.path.join("results/CorcesATAC/dea_limma/normCQN_OvA_cell_type/results.csv"),
        enrichment_results = os.path.join("results/CorcesATAC/enrichment_analysis/cell_types_up/GREAT/Azimuth_2023/cell_types_up_Azimuth_2023_all.csv"),
        crossprediction_adj_mtx = os.path.join("results/CorcesATAC/special_analyses/crossprediction/adjacency_matrix.csv"),
    output:
        umap_plot = os.path.join("paper/CorcesATAC/umap.pdf"),
        dea_heatmap_plot = os.path.join("paper/CorcesATAC/differential_heatmap.pdf"),
        enrichment_plot = os.path.join("paper/CorcesATAC/enrichment.pdf"),
        crossprediction_plot = os.path.join("paper/CorcesATAC/crossprediction.pdf"),
        crossprediction_coordinates = os.path.join("paper/CorcesATAC/crossprediction_coordinates.csv"),
    params:
        figure_theme_path = workflow.source_path("../scripts/figure_theme.R"),
        figure_utils_path = workflow.source_path("../scripts/figure_utils.R"),
        fdr_threshold = 0.05,
        log2FC_threshold = 3,
        lineage_tree_cut_off = 0.2,
        hierarchy_coordinates = True
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","CorcesATAC_figures.log"),
    script:
        "../scripts/CorcesATAC_figures.R"

# #### CorcesINT - Figure panels (custom rule) ####
rule CorcesINT_figures:
    input:
        unintegrated_umap_coords = os.path.join("results/CorcesINT/unsupervised_analysis/normupperquartile/UMAP/UMAP_correlation_15_0.1_2_data.csv"),
        integrated_umap_coords = os.path.join("results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"),
        unintegrated_cfa = os.path.join("results/CorcesINT/spilterlize_integrate/all/normupperquartile_CFA.csv"),
        integrated_cfa = os.path.join("results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated_CFA.csv"),
        norm_counts = os.path.join("results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated.csv"),
        metadata = os.path.join("results/CorcesINT/spilterlize_integrate/all/annotation.csv"),
        dea_results = os.path.join("results/CorcesINT/dea_limma/normupperquartile_integrated/results.csv"),
        gene_annotation = os.path.join("results/CorcesRNA/rnaseq_pipeline/counts/gene_annotation.csv"),
        GO_enrichment_results = os.path.join("results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/GO_Biological_Process_2025/cell_types_GO_Biological_Process_2025_all.csv"),
        Reactome_enrichment_results = os.path.join("results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/ReactomePathways/cell_types_ReactomePathways_all.csv"),
        Mono_TF_EP = os.path.join("results/CorcesINT/enrichment_analysis/Mono_EP/RcisTarget/hg38_500bp_up_100bp_down_v10clust/Mono_EP_hg38_500bp_up_100bp_down_v10clust.csv"),
        Mono_TF_TA = os.path.join("results/CorcesINT/enrichment_analysis/Mono_TA/RcisTarget/hg38_500bp_up_100bp_down_v10clust/Mono_TA_hg38_500bp_up_100bp_down_v10clust.csv"),
        HSC_TF_EP = os.path.join("results/CorcesINT/enrichment_analysis/HSC_EP/RcisTarget/hg38_500bp_up_100bp_down_v10clust/HSC_EP_hg38_500bp_up_100bp_down_v10clust.csv"),
        HSC_TF_TA = os.path.join("results/CorcesINT/enrichment_analysis/HSC_TA/RcisTarget/hg38_500bp_up_100bp_down_v10clust/HSC_TA_hg38_500bp_up_100bp_down_v10clust.csv"),
    output:
        unintegrated_cfa_plot = os.path.join("paper/CorcesINT/unintegrated_cfa.pdf"),
        integrated_cfa_plot = os.path.join("paper/CorcesINT/integrated_cfa.pdf"),
        integrated_umap_plot = os.path.join("paper/CorcesINT/integrated_umap.pdf"),
        unintegrated_umap_plot = os.path.join("paper/CorcesINT/unintegrated_umap.pdf"),
        epigenetic_scatter_dir = directory(os.path.join("paper/CorcesINT/correlation_plots")),
        GO_enrichment_plot = os.path.join("paper/CorcesINT/GO_enrichment.pdf"),
        Reactome_enrichment_plot = os.path.join("paper/CorcesINT/Reactome_enrichment.pdf"),
        TF_plot = os.path.join("paper/CorcesINT/TF_lollipop.pdf"),
    params:
        figure_theme_path = workflow.source_path("../scripts/figure_theme.R"),
        figure_utils_path = workflow.source_path("../scripts/figure_utils.R"),
        # thresholds for EP/TA categorization
        adjp_th = 0.05,
        fdr_threshold = 0.05,
        lfc_th = 1,
        ave_expr_th = 0,
        max_genes_tf_plot = 25
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","CorcesINT_figures.log"),
    script:
        "../scripts/CorcesINT_figures.R"

#### Papalexi - Figure panels (custom rule) ####
rule Papalexi_figures:
    input:
        CORRECTED_umap_coords = os.path.join("results/Papalexi2021scCRISPR/unsupervised_analysis/merged_CORRECTED/UMAP/UMAP_correlation_10_0.1_2_data.csv"),
        CORRECTED_metadata = os.path.join("results/Papalexi2021scCRISPR/scrnaseq_processing_seurat/merged/CORRECTED/metadata.csv"),
        MIXSCAPE_umap_coords = os.path.join("results/Papalexi2021scCRISPR/unsupervised_analysis/merged_MIXSCAPE_LDA/UMAP/UMAP_correlation_10_0.1_2_data.csv"),
        MIXSCAPE_metadata = os.path.join("results/Papalexi2021scCRISPR/mixscape_seurat/merged/FILTERED_metadata.csv"),
        KO_crossprediction_adj_mtx = os.path.join("results/Papalexi2021scCRISPR/special_analyses/crossprediction/adjacency_matrix.csv"),
        SPI1_TA_results = os.path.join("results/Papalexi2021scCRISPR/enrichment_analysis/SPI1/preranked_GSEApy/Corces_TA_signatures/SPI1_Corces_TA_signatures.csv"),
        KO_DEA_results = os.path.join("results/Papalexi2021scCRISPR/dea_seurat/KO_mixscape/results.csv"),
        KO_enrichment_results_GOBP = os.path.join("results/Papalexi2021scCRISPR/enrichment_analysis/KO/preranked_GSEApy/GO_Biological_Process_2025/KO_GO_Biological_Process_2025_all.csv"),
        KO_enrichment_results_Reactome = os.path.join("results/Papalexi2021scCRISPR/enrichment_analysis/KO/preranked_GSEApy/ReactomePathways/KO_ReactomePathways_all.csv"),
    output:
        umap_corrected_KO_fig = os.path.join("paper/Papalexi/umap_CORRECTED_KO.pdf"),
        umap_corrected_phase_fig = os.path.join("paper/Papalexi/umap_CORRECTED_phase.pdf"),
        umap_corrected_fig = os.path.join("paper/Papalexi/umap_CORRECTED.pdf"),
        umap_lda_fig = os.path.join("paper/Papalexi/umap_LDA.pdf"),
        crossprediction_fig = os.path.join("paper/Papalexi/crossprediction.pdf"),
        spi1_ta_lollipop_fig = os.path.join("paper/Papalexi/SPI1_TA_lollipop.pdf"),
        ko_DEA_heatmap_fig = os.path.join("paper/Papalexi/KO_DEA_heatmap.pdf"),
        ko_enrichment_GOBP_fig = os.path.join("paper/Papalexi/KO_enrichment_GOBP_bubbleplot.pdf"),
    params:
        figure_theme_path = workflow.source_path("../scripts/figure_theme.R"),
        figure_utils_path = workflow.source_path("../scripts/figure_utils.R"),
        ko_column = "gene",
        phase_column = "Phase",
        fdr_threshold = 0.05
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","Papalexi_figures.log"),
    script:
        "../scripts/Papalexi_figures.R"
