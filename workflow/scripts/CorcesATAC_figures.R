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
# ## input
# enrichment_results_path <- snakemake@input[["enrichment_results"]]
# ## output
# enrichment_plot_path <- snakemake@output[["enrichment_plot"]]
# ## params
# fdr_threshold <- snakemake@params[["fdr_threshold"]]

# input
CorcesATAC_umap_coords_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesATAC/unsupervised_analysis/normCQN_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"
CorcesATAC_dea_OvA_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesATAC/dea_limma/normCQN_OvA_cell_type/results.csv"
CorcesATAC_enrichment_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesATAC/enrichment_analysis/cell_types_up/GREAT/Azimuth_2023/cell_types_up_Azimuth_2023_all.csv"
CorcesATAC_crossprediction_adj_mtx_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesATAC/special_analyses/crossprediction/adjacency_matrix.csv"
# output
atac_umap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesATAC/umap.pdf"
atac_dea_heatmap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesATAC/differential_heatmap.pdf"
atac_enrichment_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesATAC/enrichment.pdf"
atac_crossprediction_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesATAC/crossprediction.pdf"
atac_crossprediction_coordinates_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesATAC/crossprediction_coordinates.csv"
# params
fdr_threshold <- 0.05
log2FC_threshold <- 3
lineage_tree_cut_off <- 0.2
hierarchy_coordinates <- TRUE
######### UMAPs (from unsupervised analysis) ############
# Create UMAP plots
atac_umap_plot <- umap_plot(CorcesATAC_umap_coords_path, atac_umap_path, title = "ATAC")


######### DEA HEATMAP ############
atac_dea_heatmap_plot <- plot_differential_features_heatmap(
    dea_results_path = CorcesATAC_dea_OvA_path,
    fig_path = atac_dea_heatmap_path,
    fdr_threshold = fdr_threshold,
    log2FC_threshold = log2FC_threshold,
    title = "ATAC",
    feature = 'Regions',
    ct_clst_dist = "euclidean",
    ct_clst_method = "ward.D2",
    feature_clst_dist = "maximum",
    feature_clst_method = "ward.D2",
    q_mask = 0.025,
    label_box_size_factor = 1
)


######### ENRICHMENT HEATMAP ############
create_atac_enrichment_df <- function(enrichment_results_path, fdr_threshold = 0.05) {
    # Load enrichment analysis result
    df <- data.frame(fread(file.path(enrichment_results_path), header=TRUE))
    
    # Adapt for ATAC data
    df_formatted <- df %>%
        rename(Term = description, statistic = p_adjust, score = fold_enrichment) %>%
        mutate(name = sub("_up$", "", name),
                score = ifelse(is.infinite(log2(score)), NaN, log2(score))) %>%
        mutate(name = recode(name, !!!DATA_TO_CELL_TYPE_COLORS_MAPPING))
    
    return(df_formatted)
}

atac_df_formatted <- create_atac_enrichment_df(CorcesATAC_enrichment_results_path, fdr_threshold)
atac_heatmap_df <- prepare_for_heatmap(df_formatted = atac_df_formatted, fdr_threshold = fdr_threshold)
atac_enrichment_plot <- plot_enrichment_heatmap(
    heatmap_df = atac_heatmap_df, 
    fig_path = atac_enrichment_path,
    fill_lab = "log2(FE)",
    size_lab = "-log10(p-adj.)",
    title = "ATAC",
    ylabel = "Enrichment term\n(GREAT, Azimuth 2023)"
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

atac_crosspred_p <- plot_crossprediction_from_adjacency(
    adjacency_matrix_path = CorcesATAC_crossprediction_adj_mtx_path,
    fig_path = atac_crossprediction_path,
    coordinates_out_path = atac_crossprediction_coordinates_path,
    lineage_tree_cut_off = lineage_tree_cut_off,
    hierarchy_coordinates = hierarchy_coordinates,
    modality_label = "ATAC",
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