# set correct working directory -> project/repo folder (only for dev, not prod)
# getwd()
# setwd('..')
# getwd()

#### libraries
# source this for libraries, MrBiomics theme and utility functions
source(snakemake@params[["figure_theme_path"]])
source(snakemake@params[["figure_utils_path"]])

# FIXME snakemakeify
# ## libraries are loaded in the source files
# source(snakemake@params[["figure_theme_path"]])
# ## input
# enrichment_results_path <- snakemake@input[["enrichment_results"]]
# ## output
# enrichment_plot_path <- snakemake@output[["enrichment_plot"]]
# ## params
# fdr_threshold <- snakemake@params[["fdr_threshold"]]

# inputs
CorcesATAC_umap_coords_path <- snakemake@input[["umap_coords"]]
CorcesATAC_dea_OvA_path <- snakemake@input[["dea_ova"]]
CorcesATAC_enrichment_results_path <- snakemake@input[["enrichment_results"]]
CorcesATAC_crossprediction_adj_mtx_path <- snakemake@input[["crossprediction_adj_mtx"]]

# outputs
atac_umap_path <- snakemake@output[["umap_plot"]]
atac_dea_heatmap_path <- snakemake@output[["dea_heatmap_plot"]]
atac_enrichment_path <- snakemake@output[["enrichment_plot"]]
atac_crossprediction_path <- snakemake@output[["crossprediction_plot"]]
atac_crossprediction_coordinates_path <- snakemake@output[["crossprediction_coordinates"]]

# params
fdr_threshold <- snakemake@params[["fdr_threshold"]]
log2FC_threshold <- snakemake@params[["log2FC_threshold"]]
lineage_tree_cut_off <- snakemake@params[["lineage_tree_cut_off"]]
hierarchy_coordinates <- snakemake@params[["hierarchy_coordinates"]]
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
atac_heatmap_df <- prepare_for_heatmap(df_formatted = atac_df_formatted,
                                      fdr_threshold = fdr_threshold,
                                      tissues_to_keep = c("PBMC", "Bone Marrow"))
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


atac_crosspred_p <- plot_crossprediction_from_adjacency(
    adjacency_matrix_path = CorcesATAC_crossprediction_adj_mtx_path,
    fig_path = atac_crossprediction_path,
    coordinates_out_path = atac_crossprediction_coordinates_path,
    lineage_tree_cut_off = lineage_tree_cut_off,
    hierarchy_coordinates = hierarchy_coordinates,
    modality_label = "ATAC",
)