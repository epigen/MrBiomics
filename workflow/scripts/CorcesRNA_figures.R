
# set correct working directory -> project/repo folder (only for dev, not prod)
# getwd()
# setwd('..')
# getwd()

#### libraries
# source this for libraries, MrBiomics theme and utility functions
source(snakemake@params[["figure_theme_path"]])
source(snakemake@params[["figure_utils_path"]])

# inputs
CorcesRNA_umap_coords_path <- snakemake@input[["umap_coords"]]
CorcesRNA_dea_OvA_path <- snakemake@input[["dea_ova"]]
CorcesRNA_enrichment_results_azimuth_path <- snakemake@input[["enrichment_results_azimuth"]]
CorcesRNA_enrichment_results_reactome_path <- snakemake@input[["enrichment_results_reactome"]]
CorcesRNA_crossprediction_adj_mtx_path <- snakemake@input[["crossprediction_adj_mtx"]]

# outputs
rna_umap_path <- snakemake@output[["umap_plot"]]
rna_dea_heatmap_path <- snakemake@output[["dea_heatmap_plot"]]
rna_enrichment_azimuth_path <- snakemake@output[["enrichment_azimuth_plot"]]
rna_enrichment_reactome_path <- snakemake@output[["enrichment_reactome_plot"]]
rna_crossprediction_path <- snakemake@output[["crossprediction_plot"]]
rna_crossprediction_coordinates_path <- snakemake@output[["crossprediction_coordinates"]]

# params
fdr_threshold <- snakemake@params[["fdr_threshold"]]
log2FC_threshold <- snakemake@params[["log2FC_threshold"]]
lineage_tree_cut_off <- snakemake@params[["lineage_tree_cut_off"]]
hierarchy_coordinates <- snakemake@params[["hierarchy_coordinates"]]

######### UMAPs (from unsupervised analysis) ############
# Create UMAP plots
rna_umap_plot <- umap_plot(CorcesRNA_umap_coords_path, rna_umap_path, title = "Gene expression")

######### DEA HEATMAP ############
rna_dea_heatmap_plot <- plot_differential_features_heatmap(
    dea_results_path = CorcesRNA_dea_OvA_path,
    fig_path = rna_dea_heatmap_path,
    fdr_threshold = fdr_threshold,
    log2FC_threshold = log2FC_threshold,
    title = "Gene expression",
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
rna_df_formatted_azimuth <- create_rna_enrichment_df(CorcesRNA_enrichment_results_azimuth_path, fdr_threshold)
rna_heatmap_df_azimuth <- prepare_for_heatmap(df_formatted = rna_df_formatted_azimuth,
                                      fdr_threshold = fdr_threshold,
                                      tissues_to_keep = c("PBMC", "Bone Marrow"))
rna_enrichment_azimuth_plot <- plot_enrichment_heatmap(
    heatmap_df = rna_heatmap_df_azimuth,
    fig_path = rna_enrichment_azimuth_path,
    fill_lab = "NES",
    size_lab = "-log10(p-adj.)",
    title = "Gene expression",
    ylabel = "Enrichment term\n(preranked GSEA, Azimuth 2023)"
)

rna_df_formatted_reactome <- create_rna_enrichment_df(CorcesRNA_enrichment_results_reactome_path, fdr_threshold)
rna_heatmap_df_reactome <- prepare_for_heatmap(df_formatted = rna_df_formatted_reactome,
                                               fdr_threshold = fdr_threshold)
rna_enrichment_reactome_plot <- plot_enrichment_heatmap(
    heatmap_df = rna_heatmap_df_reactome,
    fig_path = rna_enrichment_reactome_path,
    fill_lab = "NES",
    size_lab = "-log10(p-adj.)",
    title = "Gene expression",
    ylabel = "Enrichment term\n(preranked GSEA, Reactome)"
)


######### Lineage reconstruction using crossprediction ############
# plot adapted from: https://gist.github.com/dsparks/4331058

crosspred_p <- plot_crossprediction_from_adjacency(
    adjacency_matrix_path = CorcesRNA_crossprediction_adj_mtx_path,
    fig_path = rna_crossprediction_path,
    coordinates_out_path = rna_crossprediction_coordinates_path,
    lineage_tree_cut_off = lineage_tree_cut_off,
    hierarchy_coordinates = hierarchy_coordinates,
    modality_label = "Gene expression"
)