#### libraries
# source this for libraries, MrBiomics theme and utility functions
source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# Root of the repository on the current machine
# Change if running locally: 
# repo_root <- "/Users/rbednarsky/projects/MrBiomics"
repo_root <- "/nobackup/lab_bock/projects/MrBiomics"

# Inputs
CORRECTED_umap_coords_path  <- file.path(repo_root, "results/Papalexi2021scCRISPR/unsupervised_analysis/merged_CORRECTED/UMAP/UMAP_correlation_10_0.1_2_data.csv")
CORRECTED_metadata_path     <- file.path(repo_root, "results/Papalexi2021scCRISPR/scrnaseq_processing_seurat/merged/CORRECTED/metadata.csv")
MIXSCAPE_umap_coords_path   <- file.path(repo_root, "results/Papalexi2021scCRISPR/unsupervised_analysis/merged_MIXSCAPE_LDA/UMAP/UMAP_correlation_10_0.1_2_data.csv")
MIXSCAPE_metadata_path      <- file.path(repo_root, "results/Papalexi2021scCRISPR/mixscape_seurat/merged/FILTERED_metadata.csv")
KO_crossprediction_adj_mtx_path <- file.path(repo_root, "results/Papalexi2021scCRISPR/special_analyses/crossprediction/adjacency_matrix.csv")
SPI1_TA_results_path   <- file.path(repo_root, "results/Papalexi2021scCRISPR/enrichment_analysis/SPI1/preranked_GSEApy/Corces_TA_signatures/SPI1_Corces_TA_signatures.csv")
KO_DEA_results_path <- file.path(repo_root, "results/Papalexi2021scCRISPR/dea_seurat/KO_mixscape/results.csv")
# KO enrichment (both DBs)
KO_enrichment_results_path_GOBP <- file.path(repo_root, "results/Papalexi2021scCRISPR/enrichment_analysis/KO/preranked_GSEApy/GO_Biological_Process_2025/KO_GO_Biological_Process_2025_all.csv")
KO_enrichment_results_path_Reactome <- file.path(repo_root, "results/Papalexi2021scCRISPR/enrichment_analysis/KO/preranked_GSEApy/ReactomePathways/KO_ReactomePathways_all.csv")

# Outputs
umap_corrected_KO_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED_KO.pdf")
umap_corrected_phase_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED_phase.pdf")
umap_corrected_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED.pdf")
umap_lda_fig_path        <- file.path(repo_root, "paper/Papalexi/umap_LDA.pdf")
crossprediction_fig_path <- file.path(repo_root, "paper/Papalexi/crossprediction.pdf")
spi1_ta_lollipop_fig_path <- file.path(repo_root, "paper/Papalexi/SPI1_TA_lollipop.pdf")
ko_DEA_heatmap_fig_path <- file.path(repo_root, "paper/Papalexi/KO_DEA_heatmap.pdf")
# KO enrichment outputs include DB substring
ko_enrichment_GOBP_fig_path <- file.path(repo_root, "paper/Papalexi/KO_enrichment_GOBP_bubbleplot.pdf")
# ko_enrichment_Reactome_fig_path <- file.path(repo_root, "paper/Papalexi/KO_enrichment_Reactome_bubbleplot.pdf")

dir.create(dirname(umap_corrected_KO_fig_path), recursive = TRUE, showWarnings = FALSE)

ko_column <- "gene"
phase_column <- "Phase"
fdr_threshold <- 0.05

umap_corrected_panels_plot <- umap_panels_ko_and_phase_highlights(
    data_path = CORRECTED_umap_coords_path,
    metadata_path = CORRECTED_metadata_path,
    fig_path = umap_corrected_fig_path,
    ko_column = ko_column,
    phase_column = phase_column,
    point_size_all = 0.000001,
    point_size_highlight = 0.000001,
)

umap_lda_plot <- umap_plot_with_metadata(
    data_path = MIXSCAPE_umap_coords_path,
    metadata_path = MIXSCAPE_metadata_path,
    fig_path = umap_lda_fig_path,
    category_col = ko_column,
    title = "LDA-transformed Mixscape perturbation signatures\nof 11 KOs with phenotype (colored by KO)",
    min_points_for_label = 1
)

crosspred_p <- plot_crossprediction_for_kos(
    adjacency_matrix_path = KO_crossprediction_adj_mtx_path,
    fig_path = crossprediction_fig_path,
    cut_off = 0.2,
    use_hierarchy_layout = FALSE,
    label = "Functional similarity graph of\nKO phenotypes derived by\ncross-prediction (colored by KO)"
)

# SPI1 KO enrichment in Corces cell type-specific TA signatures (lollipop)
spi1_ta_lollipop_plot <- plot_ko_ta_lollipop(
    results_csv_path = SPI1_TA_results_path,
    fig_path = spi1_ta_lollipop_fig_path,
    title = "SPI1 KO signature enrichment in\ntranscriptional abundance gene sets"
)

# KO logFC heatmap from Mixscape DEA results
ko_mixscape_heatmap_plot <- plot_ko_logfc_heatmap(
    dea_results_path = KO_DEA_results_path,
    fig_path = ko_DEA_heatmap_fig_path,
    adj.P.Val_col = "p_val_adj",
    logFC_col = "avg_log2FC",
    group_renaming_map = NULL,
    fdr_threshold = 0.05,
    log2FC_threshold = 1,
    title = "DEA of knockouts",
    label_box_size_factor = 0.9,
    q_mask = 0.025, 
    n_clusters = 100
)


create_papalexi_enrichment_df <- function(enrichment_results_path, fdr_threshold = 0.05) {
    # Load enrichment analysis result
    df <- data.frame(fread(file.path(enrichment_results_path), header=TRUE))
    df_formatted <- df %>%
        rename(statistic = FDR_q_val, score = NES)
    return(df_formatted)
}

# KO enrichment heatmaps: GOBP and Reactome

# GO Biological Process 2025
ko_df_formatted_gobp <- create_papalexi_enrichment_df(KO_enrichment_results_path_GOBP, fdr_threshold)

ko_heatmap_df_gobp <- prepare_for_heatmap(
    df_formatted = ko_df_formatted_gobp,
    fdr_threshold = fdr_threshold,
    top_n_per_name = 1,
    dict_to_sort_by = PAPALEXI_KO_COLORS, 
    pos_and_neg = TRUE
)

ko_heatmap_df_gobp %>% filter(grepl("Antigen", Term)) %>% as.data.frame() 

ko_enrichment_plot_gobp <- plot_clustered_enrichment_heatmap(
    heatmap_df = ko_heatmap_df_gobp,
    fig_path = ko_enrichment_GOBP_fig_path,
    fill_lab = "NES",
    size_lab = "-log10(q-adj.)",
    title = "",
    ylabel = "Enrichment term\n(preranked GSEA, GOBP 2025)",
    xlabel = "Knockout",
    ct_clst_dist = "euclidean",
    ct_clst_method = "ward.D2",
    n_clusters = 4,
    terms_to_drop = c(  # remove a few that are duplicates in the message they send
        "Cytoplasmic Translation",
        "Ribosomal Small Subunit Biogenesis",
        "Antigen Processing and Presentation of Exogenous Peptide Antigen via MHC Class II",
        "Antigen Processing and Presentation of Peptide Antigen via MHC Class II",
        "Regulation of Immune Response"
    ),
    mask_tiny_p_values = FALSE
)

# # Reactome Pathways
# ko_df_formatted_reactome <- create_papalexi_enrichment_df(KO_enrichment_results_path_Reactome, fdr_threshold)
# ko_heatmap_df_reactome <- prepare_for_heatmap(
#     df_formatted = ko_df_formatted_reactome,
#     fdr_threshold = fdr_threshold,
#     top_n_per_name = 1,
#     dict_to_sort_by = PAPALEXI_KO_COLORS,
#     pos_and_neg = TRUE
# )
# ko_enrichment_plot_reactome <- plot_clustered_enrichment_heatmap(
#     heatmap_df = ko_heatmap_df_reactome,
#     fig_path = ko_enrichment_Reactome_fig_path,
#     fill_lab = "NES",
#     size_lab = "-log10(q-adj.)",
#     title = "",
#     ylabel = "Enrichment term\n(preranked GSEA, Reactome)",
#     xlabel = "Knockout",
#     ct_clst_dist = "euclidean",
#     ct_clst_method = "ward.D2", 
#     n_clusters = 1,
#     mask_tiny_p_values = FALSE
# )