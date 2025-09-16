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
# KO enrichment in Corces TA signatures
# KO_TA_all_results_path <- file.path(repo_root, "results/Papalexi2021scCRISPR/enrichment_analysis/KO/preranked_GSEApy/Corces_TA_signatures/KO_Corces_TA_signatures_all.csv")
SPI1_TA_results_path   <- file.path(repo_root, "results/Papalexi2021scCRISPR/enrichment_analysis/SPI1/preranked_GSEApy/Corces_TA_signatures/SPI1_Corces_TA_signatures.csv")

# Outputs
umap_corrected_KO_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED_KO.pdf")
umap_corrected_phase_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED_phase.pdf")
umap_corrected_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED.pdf")
umap_lda_fig_path        <- file.path(repo_root, "paper/Papalexi/umap_LDA.pdf")
crossprediction_fig_path <- file.path(repo_root, "paper/Papalexi/crossprediction.pdf")
spi1_ta_lollipop_fig_path <- file.path(repo_root, "paper/Papalexi/SPI1_TA_lollipop.pdf")
# ko_ta_bubble_fig_path <- file.path(repo_root, "paper/Papalexi/KO_TA_TFs_bubble.pdf")

dir.create(dirname(umap_corrected_KO_fig_path), recursive = TRUE, showWarnings = FALSE)

ko_column <- "gene"
phase_column <- "Phase"

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

# # KO TA enrichment bubble plot for selected TFs (SPI1, STAT5A, IRF1, ETV7)  -> only significant for SPI1, not useful
# selected_tfs <- c("SPI1", "STAT5A", "IRF1", "ETV7")

# ko_ta_all_df <- data.frame(fread(file.path(KO_TA_all_results_path), header = TRUE))

# # Prepare dataframe compatible with plot_enrichment_heatmap
# ko_ta_bubble_df <- ko_ta_all_df %>%
#     filter(name %in% selected_tfs) %>%
#     transmute(
#         name = as.character(name),      # TF on x-axis
#         Term = as.character(Term),      # Cell type on y-axis
#         score = NES,                    # color by NES
#         statistic = FDR_q_val           # size by -log10(q-adj.)
#     ) %>%
#     mutate(
#         sig = (!is.na(statistic)) & is.finite(statistic) & (statistic < 0.05),
#         neg_log10_statistic = -log10(statistic)
#     )

# # Cap infinite sizes if any q-adj. are zero
# if (any(is.infinite(ko_ta_bubble_df$neg_log10_statistic), na.rm = TRUE)) {
#     max_finite <- suppressWarnings(max(ko_ta_bubble_df$neg_log10_statistic[is.finite(ko_ta_bubble_df$neg_log10_statistic)], na.rm = TRUE))
#     if (is.finite(max_finite)) {
#         ko_ta_bubble_df$neg_log10_statistic[!is.finite(ko_ta_bubble_df$neg_log10_statistic)] <- max_finite * 1.1
#     }
# }

# # Order axes
# ko_ta_bubble_df$name <- factor(ko_ta_bubble_df$name, levels = selected_tfs)
# if (all(unique(ko_ta_bubble_df$Term) %in% names(CELL_TYPE_COLORS))) {
#     ko_ta_bubble_df$Term <- factor(ko_ta_bubble_df$Term, levels = names(CELL_TYPE_COLORS))
# }

# ko_ta_bubble_plot <- plot_enrichment_heatmap(
#     heatmap_df = ko_ta_bubble_df,
#     fig_path = ko_ta_bubble_fig_path,
#     fill_lab = "NES",
#     size_lab = "-log10(q-adj.)",
#     title = "KO signatures in TA gene sets",
#     ylabel = "Cell type"
# )

# SPI1 KO enrichment in Corces cell type-specific TA signatures (lollipop)
spi1_ta_lollipop_plot <- plot_ko_ta_lollipop(
    results_csv_path = SPI1_TA_results_path,
    fig_path = spi1_ta_lollipop_fig_path,
    title = "SPI1 KO signature enrichment in\ntranscriptional abundance gene sets"
)
