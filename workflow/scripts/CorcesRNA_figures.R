
# set correct working directory -> project/repo folder (only for dev, not prod)
# getwd()
# setwd('..')
# getwd()

#### libraries
# source this for libraries, MrBiomics theme and utility functions
source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# FIXME snakemakeify
# ## libraries are loaded in the source files
# source(snakemake@params[["figure_theme_path"]])
# figure_utils.R
# ## input
# enrichment_results_path <- snakemake@input[["enrichment_results"]]
# crossprediction_adj_mtx_path <- snakemake@input[["crossprediction_adj_mtx"]]
# ## output
# enrichment_plot_path <- snakemake@output[["enrichment_plot"]]
# crossprediction_plot_path <- snakemake@output[["crossprediction_plot"]]
# ## params
# # enrichment analysis
# fdr_threshold <- snakemake@params[["fdr_threshold"]]
# # lineage reconstructions
# cut_off <- snakemake@params[["cut_off"]]

# input
CorcesRNA_umap_coords_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/unsupervised_analysis/normCQN_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"
CorcesRNA_dea_OvA_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/dea_limma/normCQN_OvA_cell_type/results.csv"
CorcesRNA_enrichment_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_all.csv"
CorcesRNA_crossprediction_adj_mtx_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"
# output
rna_umap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/umap.pdf"
rna_dea_heatmap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/differential_heatmap.pdf"
rna_enrichment_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/enrichment.pdf"
rna_crossprediction_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/crossprediction.pdf"
rna_crossprediction_coordinates_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/crossprediction_coordinates.csv"
# params
# enrichment analysis
fdr_threshold <- 0.05
log2FC_threshold <- 3
lineage_tree_cut_off <- 0.07
hierarchy_coordinates <- TRUE

######### UMAPs (from unsupervised analysis) ############
# Create UMAP plots
rna_umap_plot <- umap_plot(CorcesRNA_umap_coords_path, rna_umap_path, title = "RNA")


######### DEA HEATMAP ############
rna_dea_heatmap_plot <- plot_differential_features_heatmap(
    dea_results_path = CorcesRNA_dea_OvA_path,
    fig_path = rna_dea_heatmap_path,
    fdr_threshold = fdr_threshold,
    log2FC_threshold = log2FC_threshold,
    title = "RNA",
    feature = 'Genes',
    ct_clst_dist = "euclidean",
    ct_clst_method = "ward.D2",
    feature_clst_dist = "maximum",  # maximum to focus on the most differentially expressed groups when sorting genes
    feature_clst_method = "ward.D2",
    q_mask = 0.025,
    label_box_size_factor = 1
)

######### ENRICHMENT HEATMAP ############
# Function to create enrichment heatmap
create_rna_enrichment_df <- function(enrichment_results_path, fdr_threshold = 0.05) {
    # Load enrichment analysis result
    df <- data.frame(fread(file.path(enrichment_results_path), header=TRUE))

    df_formatted <- df %>%
        rename(statistic = FDR_q_val, score = NES) %>%
        mutate(name = recode(name, !!!DATA_TO_CELL_TYPE_COLORS_MAPPING))
    
    return(df_formatted)
}

# Create enrichment heatmaps
rna_df_formatted <- create_rna_enrichment_df(CorcesRNA_enrichment_results_path, fdr_threshold)
rna_heatmap_df <- prepare_for_heatmap(df_formatted = rna_df_formatted, fdr_threshold = fdr_threshold)
rna_enrichment_plot <- plot_enrichment_heatmap(
    heatmap_df = rna_heatmap_df,
    fig_path = rna_enrichment_path,
    fill_lab = "NES",
    size_lab = "-log10(q-adj.)",
    title = "RNA",
    ylabel = "Enrichment term\n(preranked GSEA, Azimuth 2023)"
)

######### Lineage reconstruction using crossprediction ############
# plot adapted from: https://gist.github.com/dsparks/4331058

# set plot parameters
node_shape <- 19
stroke_max <- 5
point_size <- 12
fontsize <- 3
spacing_between_layers <- 5
layer_jitter <- 1
outcome_title <- "Compared to\nCorces et al. (2016)"
tp_name <- "Consistent"
fp_name <- "Additional"
fn_name <- "Missing"

crosspred_p <- plot_crossprediction_from_adjacency(
    adjacency_matrix_path = CorcesRNA_crossprediction_adj_mtx_path,
    fig_path = rna_crossprediction_path,
    coordinates_out_path = rna_crossprediction_coordinates_path,
    lineage_tree_cut_off = lineage_tree_cut_off,
    hierarchy_coordinates = hierarchy_coordinates,
    modality_label = "RNA",
    root_node = "HSC",
    spacing_between_layers = spacing_between_layers,
    jitter_step = layer_jitter,
    tp_name = tp_name,
    fp_name = fp_name,
    fn_name = fn_name,
    outcome_title = outcome_title,
    node_shape = node_shape,
    point_size = point_size,
    fontsize = fontsize,
    stroke_max = stroke_max,
    fig_width = 6
)
