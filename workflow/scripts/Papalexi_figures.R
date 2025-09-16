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
KO_mixscape_results_path <- file.path(repo_root, "results/Papalexi2021scCRISPR/dea_seurat/KO_mixscape/results.csv")

# Outputs
umap_corrected_KO_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED_KO.pdf")
umap_corrected_phase_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED_phase.pdf")
umap_corrected_fig_path  <- file.path(repo_root, "paper/Papalexi/umap_CORRECTED.pdf")
umap_lda_fig_path        <- file.path(repo_root, "paper/Papalexi/umap_LDA.pdf")
crossprediction_fig_path <- file.path(repo_root, "paper/Papalexi/crossprediction.pdf")
spi1_ta_lollipop_fig_path <- file.path(repo_root, "paper/Papalexi/SPI1_TA_lollipop.pdf")
ko_mixscape_heatmap_fig_path <- file.path(repo_root, "paper/Papalexi/KO_mixscape_heatmap.pdf")

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

# SPI1 KO enrichment in Corces cell type-specific TA signatures (lollipop)
spi1_ta_lollipop_plot <- plot_ko_ta_lollipop(
    results_csv_path = SPI1_TA_results_path,
    fig_path = spi1_ta_lollipop_fig_path,
    title = "SPI1 KO signature enrichment in\ntranscriptional abundance gene sets"
)

# KO logFC heatmap from Mixscape DEA results
ko_mixscape_heatmap_plot <- plot_ko_logfc_heatmap(
    dea_results_path = KO_mixscape_results_path,
    fig_path = ko_mixscape_heatmap_fig_path,
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
