#### libraries
# source this for libraries, MrBiomics theme and utility functions

source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# Inputs
CORRECTED_umap_coords_path  <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/unsupervised_analysis/merged_CORRECTED/UMAP/UMAP_correlation_10_0.1_2_data.csv"
CORRECTED_metadata_path     <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/scrnaseq_processing_seurat/merged/CORRECTED/metadata.csv"
MIXSCAPE_umap_coords_path   <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/unsupervised_analysis/merged_MIXSCAPE_LDA/UMAP/UMAP_correlation_10_0.1_2_data.csv"
MIXSCAPE_metadata_path      <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/mixscape_seurat/merged/FILTERED_metadata.csv"
KO_crossprediction_adj_mtx_path <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/special_analyses/crossprediction/adjacency_matrix.csv"
# KO enrichment in Corces TA signatures
KO_TA_all_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/enrichment_analysis/KO/preranked_GSEApy/Corces_TA_signatures/KO_Corces_TA_signatures_all.csv"
SPI1_TA_results_path   <- "/nobackup/lab_bock/projects/MrBiomics/results/Papalexi2021scCRISPR/enrichment_analysis/SPI1/preranked_GSEApy/Corces_TA_signatures/SPI1_Corces_TA_signatures.csv"

# Outputs
umap_corrected_fig_path  <- "/nobackup/lab_bock/projects/MrBiomics/paper/Papalexi/umap_CORRECTED.pdf"
umap_lda_fig_path        <- "/nobackup/lab_bock/projects/MrBiomics/paper/Papalexi/umap_LDA.pdf"
crossprediction_fig_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/Papalexi/crossprediction.pdf"
spi1_ta_lollipop_fig_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/Papalexi/SPI1_TA_lollipop.pdf"

dir.create(dirname(umap_corrected_fig_path), recursive = TRUE, showWarnings = FALSE)

ko_column <- "gene"

umap_corrected_plot <- umap_plot_with_metadata(
    data_path = CORRECTED_umap_coords_path,
    metadata_path = CORRECTED_metadata_path,
    fig_path = umap_corrected_fig_path,
    category_col = ko_column,
    title = "scCRISPR-seq of 25 gene knockouts (KO) in\nstimulated THP-1 cells (colored by KO)",
    min_points_for_label = 1,
    label_points = FALSE
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
