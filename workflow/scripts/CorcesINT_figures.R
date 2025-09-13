source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# FIXME snakemakeify
#### configuration ####
# Root of the repository on the current machine
# Change if running locally: 
# repo_root <- "/Users/rbednarsky/projects/MrBiomics"
repo_root <- "/nobackup/lab_bock/projects/MrBiomics"
# input
unintegrated_umap_coords_path <- file.path(repo_root, "results/CorcesINT/unsupervised_analysis/normupperquartile/UMAP/UMAP_correlation_15_0.1_2_data.csv")
integrated_umap_coords_path <-   file.path(repo_root, "results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv")
unintegrated_cfa_path <- file.path(repo_root, "results/CorcesINT/spilterlize_integrate/all/normupperquartile_CFA.csv")
integrated_cfa_path <- file.path(repo_root, "results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated_CFA.csv")
norm_counts_path <- file.path(repo_root, "results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated.csv")
metadata_path <- file.path(repo_root, "results/CorcesINT/spilterlize_integrate/all/annotation.csv")
dea_results_path <- file.path(repo_root, "results/CorcesINT/dea_limma/normupperquartile_integrated/results.csv")
gene_annotation_path <- file.path(repo_root, "results/CorcesRNA/rnaseq_pipeline/counts/gene_annotation.csv")
GO_enrichment_results_path <- file.path(repo_root, "results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/GO_Biological_Process_2025/cell_types_GO_Biological_Process_2025_all.csv")
Reactome_enrichment_results_path <- file.path(repo_root, "results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/ReactomePathways/cell_types_ReactomePathways_all.csv")
Mono_TF_EP_path <- file.path(repo_root, "results/CorcesINT/enrichment_analysis/Mono_EP/RcisTarget/hg38_500bp_up_100bp_down_v10clust/Mono_EP_hg38_500bp_up_100bp_down_v10clust.csv")
Mono_TF_TA_path <- file.path(repo_root, "results/CorcesINT/enrichment_analysis/Mono_TA/RcisTarget/hg38_500bp_up_100bp_down_v10clust/Mono_TA_hg38_500bp_up_100bp_down_v10clust.csv")
HSC_TF_EP_path  <- file.path(repo_root, "results/CorcesINT/enrichment_analysis/HSC_EP/RcisTarget/hg38_500bp_up_100bp_down_v10clust/HSC_EP_hg38_500bp_up_100bp_down_v10clust.csv")
HSC_TF_TA_path  <- file.path(repo_root, "results/CorcesINT/enrichment_analysis/HSC_TA/RcisTarget/hg38_500bp_up_100bp_down_v10clust/HSC_TA_hg38_500bp_up_100bp_down_v10clust.csv")

# parameters
adjp_th <- 0.05
fdr_threshold <- 0.05
lfc_th <- 1
ave_expr_th <- 0
max_genes_tf_plot <- 25
TFs_in_papalexi <- c(
  "ATF2", "BRD4", "CAV1", "CD86", "CMTM6", "CUL3", "ETV7", "IFNGR1", "IFNGR2", "IRF1", "IRF7", "JAK2", "CMIR", "MARCH8",
  "MYC", "NFKB1A", "PDCD1LG2", "POU2F2", "SMAD4", "SPI1", "STAT1", "STAT2", "STAT3", "STAT5A", "TNFRSF14", "UBE2L6"
)
# output
unintegrated_cfa_plot_path <- file.path(repo_root, "paper/CorcesINT/unintegrated_cfa.pdf")
integrated_cfa_plot_path <- file.path(repo_root, "paper/CorcesINT/integrated_cfa.pdf")
integrated_umap_plot_path <- file.path(repo_root, "paper/CorcesINT/integrated_umap.pdf")
unintegrated_umap_plot_path <- file.path(repo_root, "paper/CorcesINT/unintegrated_umap.pdf")
epigenetic_scatter_dir <- file.path(repo_root, "paper/CorcesINT/correlation_plots")
GO_enrichment_path <- file.path(repo_root, "paper/CorcesINT/GO_enrichment.pdf")
Reactome_enrichment_path <- file.path(repo_root, "paper/CorcesINT/Reactome_enrichment.pdf")
TF_plot_path <- file.path(repo_root, "paper/CorcesINT/TF_lollipop.pdf")
########################################################################################################################
### LOAD DATA ##########################################################################################################
########################################################################################################################
V1_to_rowname <- function(df) {
    rownames(df) <- df$V1
    df <- df[,-1]
    return(df)
}

unintegrated_umap_coords <- data.frame(fread(file.path(unintegrated_umap_coords_path), header=TRUE))
integrated_umap_coords <- data.frame(fread(file.path(integrated_umap_coords_path), header=TRUE))
norm_counts <- data.frame(fread(file.path(norm_counts_path), header=TRUE), row.names=1)
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)
dea_results <- data.frame(fread(file.path(dea_results_path), header=TRUE))
gene_annotation <- data.frame(fread(file.path(gene_annotation_path), header=TRUE))
unintegrated_cfa_data <- V1_to_rowname(data.frame(fread(file.path(unintegrated_cfa_path), header=TRUE)))
integrated_cfa_data <- V1_to_rowname(data.frame(fread(file.path(integrated_cfa_path), header=TRUE)))

########################################################################################################################
### UMAP PLOT ##########################################################################################################
########################################################################################################################
unintegrated_umap_plot <- umap_plot(unintegrated_umap_coords_path, unintegrated_umap_plot_path, title = "Unintegrated",
                                    modality_by_shape = TRUE)
integrated_umap_plot <- umap_plot(integrated_umap_coords_path, integrated_umap_plot_path, title = "Integrated",
                                    modality_by_shape = TRUE)

########################################################################################################################
### CFA PLOT ##########################################################################################################
########################################################################################################################
plot_cfa_heatmap <- function(cfa_mat, title, path, var_max=NULL, nPCs=10, 
                             metadata_rows=c('cell_type', 'modality', 'donor')
                             ){
  var_explained_df <- t(cfa_mat['var_explained', 1:nPCs]) %>%
      as.data.frame() %>%
      rownames_to_column(var = 'PC') %>%
      mutate(PC = factor(PC, levels = rev(colnames(cfa_mat))))
  cfa_mat <- cfa_mat[metadata_rows, 1:nPCs]

  if (is.null(var_max)) {
      x_max <- max(var_explained_df$var_explained, na.rm = TRUE)
  } else {
      x_max <- var_max
  }

  # unstack cfa_mat, remembering rownames and colnames
  cfa_mat_long <- cfa_mat %>%
      rownames_to_column(var = "metadata_type") %>%
      pivot_longer(cols = -metadata_type, names_to = "PC", values_to = "stat") %>%
      mutate(PC = factor(PC, levels = rev(colnames(cfa_mat))),
              metadata_type = factor(metadata_type, levels = metadata_rows),
              text_color = ifelse(stat > (min(stat) + 0.75 * (max(stat) - min(stat))), "white", "black"),
              text_fontface = ifelse(stat > (min(stat) + 0.75 * (max(stat) - min(stat))), "bold", "plain")
              )

  # barplot of var_explained
  var_explained_plot <- ggplot(var_explained_df, aes(y = PC, x = var_explained)) +
      geom_bar(stat = "identity", fill = 'grey80') +
      MrBiomics_theme() +
      scale_x_continuous(limits = c(0, x_max), breaks = c(0, x_max/2, x_max),
                          labels = function(x) format(x, digits = 2)) +
      theme(
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(),
          panel.border = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      labs(title = NULL, x = "Variance explained", y = NULL)

  cfa_heatmap <- ggplot(cfa_mat_long, aes(x = metadata_type, y = PC, fill = stat)) +
      geom_tile(linewidth = 0) +
      geom_text(aes(label = round(stat, 1), color = text_color, fontface = text_fontface)) +
      scale_fill_gradient(low = "white", high = as.character(RdBu_extremes["up"]), name='-log10(p-adj.)\nfor association\nof PC &\nmetadata') +
      scale_color_identity(guide = "none") +
      scale_x_discrete(labels = function(x) tools::toTitleCase(gsub("_", " ", x)), expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      MrBiomics_theme() +
      theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      ) +
      coord_fixed() +
      labs(title = title, y = paste0("PC(1-", ncol(cfa_mat), ")"), x = NULL)

  cfa_plot <- cfa_heatmap + var_explained_plot + plot_spacer() + guide_area() +
      plot_layout(ncol = 4, widths = c(2, 2, 0.2, 1.4), guides = "collect") +
      plot_annotation(theme = theme(plot.margin = margin(10, 20, 10, 10))) &
      theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5), legend.box.margin = margin(5.5, 5.5, 5.5, 10),
              legend.margin = margin(5.5, 5.5, 5.5, 5.5))

  ggsave_all_formats(path = path,
                      plot = cfa_plot,
                      width = PLOT_SIZE_3_PER_ROW-0.5,
                      height = PLOT_SIZE_3_PER_ROW)
  return(cfa_plot)
}

var_max <- max(unintegrated_cfa_data['var_explained',], integrated_cfa_data['var_explained',])

unintegrated_cfa_plot <- plot_cfa_heatmap(cfa_mat=unintegrated_cfa_data, title='Unintegrated',
     path=unintegrated_cfa_plot_path, var_max=var_max, nPCs=10)
integrated_cfa_plot <- plot_cfa_heatmap(cfa_mat=integrated_cfa_data, title='Integrated',
     path=integrated_cfa_plot_path, var_max=var_max, nPCs=10)



########################################################################################################################
### EPIGENETIC POTENTIAL PLOT ##########################################################################################
########################################################################################################################
# transform norm_counts for plotting (for each group)
dea <- dea_results |> mutate(cell_type = sub("cell_type(.*?)__.*", "\\1", group))

for(ct in unique(metadata$cell_type)){

    ## ---------- summarise mean signal per gene ----------
    atac <- rownames(metadata)[metadata$cell_type == ct & metadata$modality == "ATAC"]
    rna  <- rownames(metadata)[metadata$cell_type == ct & metadata$modality == "RNA" ]
    
    df <- data.frame(
    feature = rownames(norm_counts),
    ATAC = rowMeans(norm_counts[ , atac, drop = FALSE]),
    RNA  = rowMeans(norm_counts[ , rna , drop = FALSE])
    )
    
    ## ---------- add category EP or TA ----------
    stat <- dea |>
    filter(cell_type == ct) |>
    select(feature, logFC, adj.P.Val, AveExpr)
    
    df <- df |>
    left_join(stat, by = "feature") |>
    mutate(category = case_when(
      adj.P.Val > adjp_th | abs(logFC) < lfc_th | AveExpr < ave_expr_th ~ "Convergent",
      adj.P.Val <= adjp_th & logFC >= lfc_th & AveExpr >= ave_expr_th ~ "Epigenetic potential",
      adj.P.Val <= adjp_th & logFC  <= -lfc_th & AveExpr >= ave_expr_th ~ "Transcriptional abundance"
    ))
    
    ## ---------- labels with counts and stats ----------
    counts <- table(df$category)
    
    # prepare statistics
    r_val <- cor(df$ATAC, df$RNA, method = "pearson", use = "complete.obs")

    df$log10_adjp <- -log10(df$adj.P.Val)
    df$gene_symbol <- gene_annotation$external_gene_name[match(df$feature, gene_annotation$ensembl_gene_id)]

    # prepare HAEMATOPOIESIS_MARKERS for annotation (divergent only)
    df$markers <- ifelse(
      !is.na(df$gene_symbol) & (df$gene_symbol %in% HAEMATOPOIESIS_MARKERS) & (df$category != "Convergent"),
      df$gene_symbol,
      ""
    )

    # split norm_counts for plotting
    df_conv <-  df %>% filter(category == "Convergent")
    df_div  <-  df %>% filter(category != "Convergent")

    # prepare annotations instead of legend
    cat_cols <- c("Convergent"             = "grey50",
              "Epigenetic potential"   = as.character(RdBu_extremes["up"]),
              "Transcriptional abundance"= as.character(RdBu_extremes["down"]))

    # annotation labels
    max_x <- max(df$ATAC)
    min_x <- min(df$ATAC)
    max_y <- max(df$RNA)
    min_y <- min(df$RNA)
    ann_df <- data.frame(
        category = c("Transcriptional abundance", "Convergent", "Epigenetic potential"),
        x = c(min_x,  min_x,  max_x),   # TL, BL, BR
        y = c( max_y, min_y, min_y),
        hjust = c(0, 0, 1),        # keep text inside panel
        vjust = c(1, 0, 0),
        lab = c(
            paste0("Transcriptional\nabundance\n(", counts["Transcriptional abundance"],")"),
            paste0("Convergent\n(", counts["Convergent"],")"),
            paste0("Epigenetic\npotential\n(", counts["Epigenetic potential"],")")
            )
        )
    
    # breaks for size legend based on original -log10(adjusted p)
    size_breaks <- pretty(range(df_div$log10_adjp, na.rm = TRUE), n = 4)

    ## ---------- plot ----------
    p <- ggplot(df, aes(x = ATAC, y = RNA, label = markers)) +
    # plot convergent first (small, faint – background)
    rasterise(geom_point(data = df_conv,
             aes(colour = category),
             size  = 0.5,
             alpha = 0.1,
             stroke = 0,
             show.legend = FALSE)) +
    # overlay divergent classes (larger, more opaque – foreground)
    rasterise(geom_point(data = df_div %>% filter(markers == ""),
             aes(colour = category, size = log10_adjp),
             alpha = 0.5,
             stroke = 0)) +
    # colour palette 
    scale_colour_manual(
    values = cat_cols,
    name   = NULL, 
    guide = "none"
    ) +
    # size legend using original -log10(adjusted p) values
    scale_size_continuous(
      name = "-log10(p-adj.)",
      range = c(0.1, 2),
      breaks = size_breaks,
      labels = function(x) format(x, digits = 2)
    ) +
    # annotation boxes with counts, also serving as legend
    geom_text(
    data = ann_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = lab, color = category, hjust = hjust, vjust = vjust),
    size = FONT_SIZE_NORMAL / .pt
    ) +
    # annotate divergent HAEMATOPOIESIS_MARKERS
    geom_text_repel(
        color = "black",
        box.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 42,
        show.legend = FALSE
    ) +
    # # redraw labeled points on top so labels don't occlude points
    # geom_point(
    #     data = df %>% filter(markers != ""),
    #     aes(colour = category, size = log10_adjp),
    #     alpha = 0.5,
    #     stroke = 0
    # ) +
    # axis titles
    labs(
    title = paste0(ct,"\n",sprintf("Pearson's R = %.2f", r_val)),
    x = "Chromatin accessibility\n(normalized & integrated)",
    y = "Gene expression\n(normalized & integrated)"
    ) +
    MrBiomics_theme() + 
    theme(aspect.ratio = 1)

    # print(p)
    
  ## ---------- save plot ----------
  ggsave_all_formats(file.path(epigenetic_scatter_dir,paste0(ct,"_correlation.png")), p,
                    width = PLOT_SIZE_3_PER_ROW+1,
                    height = PLOT_SIZE_3_PER_ROW
                    )
}

########################################################################################################################
### ENRICHMENT HEATMAP #################################################################################################
########################################################################################################################
create_int_enrichment_df <- function(enrichment_results_path, fdr_threshold = 0.05) {
    # Load enrichment analysis result
    df <- data.frame(fread(file.path(enrichment_results_path), header=TRUE))
    df_formatted <- df %>%
        rename(statistic = FDR_q_val, score = NES) %>%
        mutate(name = recode(name, !!!DATA_TO_CELL_TYPE_COLORS_MAPPING))
    return(df_formatted)
}

GO_df_formatted <- create_int_enrichment_df(GO_enrichment_results_path, fdr_threshold)
GO_heatmap_df <- prepare_for_heatmap(df_formatted = GO_df_formatted, fdr_threshold = fdr_threshold, top_n_per_name = 3)
GO_enrichment_plot <- plot_clustered_enrichment_heatmap(
    heatmap_df = GO_heatmap_df,
    fig_path = GO_enrichment_path,
    fill_lab = "NES",
    size_lab = "-log10(q-adj.)",
    title = "",
    ylabel = "Enrichment term\n(preranked GSEA,\nGOBP 2025)",
    ct_clst_dist = "euclidean",
    ct_clst_method = "ward.D2",
    term_clst_dist = "euclidean",
    term_clst_method = "ward.D2",
    n_clusters = 25,
    width = PLOT_SIZE_2_PER_ROW
)

Reactome_df_formatted <- create_int_enrichment_df(Reactome_enrichment_results_path, fdr_threshold)
Reactome_heatmap_df <- prepare_for_heatmap(df_formatted = Reactome_df_formatted, fdr_threshold = fdr_threshold,
                                           top_n_per_name = 2)
Reactome_enrichment_plot <- plot_clustered_enrichment_heatmap(
    heatmap_df = Reactome_heatmap_df,
    fig_path = Reactome_enrichment_path,
    fill_lab = "NES",
    size_lab = "-log10(q-adj.)",
    title = "",
    ylabel = "Enrichment term\n(preranked GSEA, Reactome)",
    ct_clst_dist = "euclidean",
    ct_clst_method = "ward.D2",
    term_clst_dist = "euclidean",
    term_clst_method = "ward.D2",
    n_clusters = 25,
    width = PLOT_SIZE_2_PER_ROW
)

dir.create(file.path(repo_root, "paper/CorcesINT/tmp"), showWarnings = FALSE)
Reactome_heatmap_df %>% filter(statistic < 0.05) %>% as.data.frame() %>% write.csv(file.path(repo_root, "paper/CorcesINT/tmp/Reactome_enrichment_sig.csv"), row.names = FALSE)
GO_heatmap_df %>% filter(statistic < 0.05) %>% as.data.frame() %>% write.csv(file.path(repo_root, "paper/CorcesINT/tmp/GO_enrichment_sig.csv"), row.names = FALSE)


########################################################################################################################
### TF PLOT ############################################################################################################
########################################################################################################################
# helper: load and preprocess TF enrichment CSV (Mono/HSC)
load_prepare_tf <- function(csv_path) {
  df <- data.frame(fread(file.path(csv_path), header=TRUE))
  df <- df %>%
    mutate(.row_id = row_number()) %>%
    select(.row_id, NES, nEnrGenes, TF_highConf) %>%
    mutate(
      TF_highConf = coalesce(TF_highConf, ""),
      TF_highConf = str_remove_all(TF_highConf, "\\([^)]*\\)"),   # drop (directAnnotation) etc.
      TF_highConf = str_replace_all(TF_highConf, ",", " "),       # remove stray commas
      TF_highConf = str_squish(TF_highConf)
    ) %>%
    # split TFs on semicolons or periods used as separators
    separate_rows(TF_highConf, sep = ";|\\.") %>%
    mutate(
      TF = str_squish(TF_highConf),
      TF = str_remove(TF, "\\.$")   # drop trailing dots if any remain
    ) %>%
    filter(TF != "") %>%
    distinct(.row_id, TF, .keep_all = TRUE) %>%
    select(TF, NES, nEnrGenes) %>%
    filter(nEnrGenes >= 30) %>%
    as.data.frame()
  return(df)
}

# helper: order TFs by max NES (ascending) with tiebreaker by sum(NES)
order_by_max_nes <- function(df) {
  tf_order <- df %>%
    group_by(TF) %>%
    summarise(
      max_NES = max(NES, na.rm = TRUE),
      sum_NES = sum(NES, na.rm = TRUE)
    ) %>%
    arrange(max_NES, sum_NES) %>%
    pull(TF)
  df$TF <- factor(df$TF, levels = tf_order)
  return(df)
}

# load & prepare Mono EP and TA
Mono_TF_EP_df <- load_prepare_tf(Mono_TF_EP_path)
Mono_TF_TA_df <- load_prepare_tf(Mono_TF_TA_path)

# load & prepare HSC EP and TA
HSC_TF_EP_df <- load_prepare_tf(HSC_TF_EP_path)
HSC_TF_TA_df <- load_prepare_tf(HSC_TF_TA_path)

# save all the dfs as csv right next to the plot path
Mono_TF_EP_df %>% write.csv(file.path(str_replace(TF_plot_path, "\\.pdf$", "__Mono_TF_EP_df.csv")), row.names = FALSE)
Mono_TF_TA_df %>% write.csv(file.path(str_replace(TF_plot_path, "\\.pdf$", "__Mono_TF_TA_df.csv")), row.names = FALSE)
HSC_TF_EP_df %>% write.csv(file.path(str_replace(TF_plot_path, "\\.pdf$", "__HSC_TF_EP_df.csv")), row.names = FALSE)
HSC_TF_TA_df %>% write.csv(file.path(str_replace(TF_plot_path, "\\.pdf$", "__HSC_TF_TA_df.csv")), row.names = FALSE)

# balance TF counts within each cell type by trimming the dataset with more TFs (lowest NES removed first)
# Mono
n_EP <- n_distinct(Mono_TF_EP_df$TF)
n_TA <- n_distinct(Mono_TF_TA_df$TF)
if (n_EP > n_TA) {
  keep_ep <- Mono_TF_EP_df %>%
    group_by(TF) %>%
    summarise(max_NES = max(NES, na.rm = TRUE)) %>%
    slice_max(order_by = max_NES, n = n_TA, with_ties = FALSE) %>%
    pull(TF)
  Mono_TF_EP_df <- Mono_TF_EP_df %>% filter(TF %in% keep_ep)
} else if (n_TA > n_EP) {
  keep_ta <- Mono_TF_TA_df %>%
    group_by(TF) %>%
    summarise(max_NES = max(NES, na.rm = TRUE)) %>%
    slice_max(order_by = max_NES, n = n_EP, with_ties = FALSE) %>%
    pull(TF)
  Mono_TF_TA_df <- Mono_TF_TA_df %>% filter(TF %in% keep_ta)
}

# HSC
n_EP_h <- n_distinct(HSC_TF_EP_df$TF)
n_TA_h <- n_distinct(HSC_TF_TA_df$TF)
if (n_EP_h > n_TA_h) {
  keep_ep_h <- HSC_TF_EP_df %>%
    group_by(TF) %>%
    summarise(max_NES = max(NES, na.rm = TRUE)) %>%
    slice_max(order_by = max_NES, n = n_TA_h, with_ties = FALSE) %>%
    pull(TF)
  HSC_TF_EP_df <- HSC_TF_EP_df %>% filter(TF %in% keep_ep_h)
} else if (n_TA_h > n_EP_h) {
  keep_ta_h <- HSC_TF_TA_df %>%
    group_by(TF) %>%
    summarise(max_NES = max(NES, na.rm = TRUE)) %>%
    slice_max(order_by = max_NES, n = n_EP_h, with_ties = FALSE) %>%
    pull(TF)
  HSC_TF_TA_df <- HSC_TF_TA_df %>% filter(TF %in% keep_ta_h)
}

# After balancing, select top N TFs by max NES within each cell type
N_matched_mono <- min(n_distinct(Mono_TF_EP_df$TF), n_distinct(Mono_TF_TA_df$TF), max_genes_tf_plot)
N_matched_hsc  <- min(n_distinct(HSC_TF_EP_df$TF),  n_distinct(HSC_TF_TA_df$TF),  max_genes_tf_plot)

select_top_by_max_nes <- function(df, N) {
  # Prioritize TFs from TFs_in_papalexi, then fill remaining slots by max NES
  tf_stats <- df %>%
    group_by(TF) %>%
    summarise(max_NES = max(NES, na.rm = TRUE))

  # 1) Take up to N TFs that are in TFs_in_papalexi, highest max NES first
  pri <- tf_stats %>%
    filter(TF %in% TFs_in_papalexi) %>%
    slice_max(order_by = max_NES, n = N, with_ties = FALSE) %>%
    pull(TF)

  remaining_n <- max(N - length(pri), 0)

  # 2) Fill the remaining with non-priority TFs by highest max NES
  if (remaining_n > 0) {
    non_pri <- tf_stats %>%
      filter(!(TF %in% pri)) %>%
      slice_max(order_by = max_NES, n = remaining_n, with_ties = FALSE) %>%
      pull(TF)
    keep <- c(pri, non_pri)
  } else {
    keep <- pri
  }

  df %>% filter(TF %in% keep)
}

Mono_TF_EP_df <- select_top_by_max_nes(Mono_TF_EP_df, N_matched_mono)
Mono_TF_TA_df <- select_top_by_max_nes(Mono_TF_TA_df, N_matched_mono)
HSC_TF_EP_df  <- select_top_by_max_nes(HSC_TF_EP_df,  N_matched_hsc)
HSC_TF_TA_df  <- select_top_by_max_nes(HSC_TF_TA_df,  N_matched_hsc)

# order y-axis for each dataset (with tiebreaker)
Mono_TF_EP_df <- order_by_max_nes(Mono_TF_EP_df)
Mono_TF_TA_df <- order_by_max_nes(Mono_TF_TA_df)
HSC_TF_EP_df  <- order_by_max_nes(HSC_TF_EP_df)
HSC_TF_TA_df  <- order_by_max_nes(HSC_TF_TA_df)

# unified scales across all four panels
max_nes_both <- max(c(Mono_TF_EP_df$NES, Mono_TF_TA_df$NES, HSC_TF_EP_df$NES, HSC_TF_TA_df$NES), na.rm = TRUE)
size_range <- c(1, 2.5)

segment_color <- "grey80"
segment_width <- 0.8

# helper: build rich y-axis labels (bold if in Papalexi; color by shared status)
red_col <- as.character(RdBu_extremes["up"])    # epigenetic potential color
blue_col <- as.character(RdBu_extremes["down"])  # transcriptional abundance color
violet_col <- "purple"

Mono_EP_TFs <- unique(as.character(Mono_TF_EP_df$TF))
Mono_TA_TFs <- unique(as.character(Mono_TF_TA_df$TF))
HSC_EP_TFs  <- unique(as.character(HSC_TF_EP_df$TF))
HSC_TA_TFs  <- unique(as.character(HSC_TF_TA_df$TF))

build_tf_axis_labels <- function(tfs, papalexi_set, red_set, blue_set) {
  vapply(tfs, function(tf) {
    is_pap <- tf %in% papalexi_set
    in_red <- tf %in% red_set
    in_blue <- tf %in% blue_set

    if (in_red && in_blue) {
      col <- violet_col
    } else if (in_red) {
      col <- red_col
    } else if (in_blue) {
      col <- blue_col
    } else {
      col <- NA_character_
    }

    lab <- tf
    if (is_pap) {
      lab <- paste0("<b>", lab, "</b>")
    }
    if (!is.na(col)) {
      lab <- paste0("<span style='color:", col, "'>", lab, "</span>")
    }
    lab
  }, character(1))
}

# EP plot (red scale) MONO
Mono_TF_EP_plot <- ggplot(Mono_TF_EP_df, aes(x = NES, y = TF, size = nEnrGenes, color = NES)) +
  geom_segment(aes(x = 0, xend = NES, y = TF, yend = TF), color = segment_color, linewidth = segment_width) +
  geom_point(alpha = 0.8, shape = 16) +
  MrBiomics_theme() +
  theme(axis.text.y = ggtext::element_markdown(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(limits = c(0, max_nes_both)) +
  scale_y_discrete(labels = function(x) build_tf_axis_labels(x, TFs_in_papalexi, HSC_EP_TFs, HSC_TA_TFs)) +
  scale_size_continuous(range = size_range, name = "Number of\nenriched\ngenes", guide = "none") +
  scale_color_gradient(limits = c(0, max_nes_both), low = "white", high = as.character(RdBu_extremes["up"]), name = "NES", guide = "none") +
  labs(x = "NES", y = NULL, title = "Mono\nEpigenetic\npotential")

# TA plot (blue scale) MONO
Mono_TF_TA_plot <- ggplot(Mono_TF_TA_df, aes(x = NES, y = TF, size = nEnrGenes, color = NES)) +
  geom_segment(aes(x = 0, xend = NES, y = TF, yend = TF), color = segment_color, linewidth = segment_width) +
  geom_point(alpha = 0.8, shape = 16) +
  MrBiomics_theme() +
  theme(axis.text.y = ggtext::element_markdown(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(limits = c(0, max_nes_both)) +
  scale_y_discrete(labels = function(x) build_tf_axis_labels(x, TFs_in_papalexi, HSC_EP_TFs, HSC_TA_TFs)) +
  scale_size_continuous(range = size_range, name = "Number of\nenriched\ngenes") +
  scale_color_gradient(limits = c(0, max_nes_both), low = "white", high = as.character(RdBu_extremes["down"]), name = "NES", guide = "none") +
  labs(x = "NES", y = NULL, title = "Mono\nTranscriptional\nabundance")

# EP plot (red scale) HSC
HSC_TF_EP_plot <- ggplot(HSC_TF_EP_df, aes(x = NES, y = TF, size = nEnrGenes, color = NES)) +
  geom_segment(aes(x = 0, xend = NES, y = TF, yend = TF), color = segment_color, linewidth = segment_width) +
  geom_point(alpha = 0.8, shape = 16) +
  MrBiomics_theme() +
  theme(axis.text.y = ggtext::element_markdown(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(limits = c(0, max_nes_both)) +
  scale_y_discrete(labels = function(x) build_tf_axis_labels(x, TFs_in_papalexi, Mono_EP_TFs, Mono_TA_TFs)) +
  scale_size_continuous(range = size_range, name = "Number of\nenriched\ngenes", guide = "none") +
  scale_color_gradient(limits = c(0, max_nes_both), low = "white", high = as.character(RdBu_extremes["up"]), name = "NES", guide = "none") +
  labs(x = "NES", y = NULL, title = "HSC\nEpigenetic\npotential")

# TA plot (blue scale) HSC
HSC_TF_TA_plot <- ggplot(HSC_TF_TA_df, aes(x = NES, y = TF, size = nEnrGenes, color = NES)) +
  geom_segment(aes(x = 0, xend = NES, y = TF, yend = TF), color = segment_color, linewidth = segment_width) +
  geom_point(alpha = 0.8, shape = 16) +
  MrBiomics_theme() +
  theme(axis.text.y = ggtext::element_markdown(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(limits = c(0, max_nes_both)) +
  scale_y_discrete(labels = function(x) build_tf_axis_labels(x, TFs_in_papalexi, Mono_EP_TFs, Mono_TA_TFs)) +
  scale_size_continuous(range = size_range, name = "Number of\nenriched\ngenes", guide = "none") +
  scale_color_gradient(limits = c(0, max_nes_both), low = "white", high = as.character(RdBu_extremes["down"]), name = "NES", guide = "none") +
  labs(x = "NES", y = NULL, title = "HSC\nTranscriptional\nabundance")

# combine with patchwork and save
# dummy plot to provide a single black-white NES colorbar legend
legend_dummy_df <- data.frame(x = 1:2, y = 1, NES = c(0, max_nes_both))
legend_dummy <- ggplot(legend_dummy_df, aes(x = x, y = y, color = NES)) +
  geom_point(alpha = 0, size = 0) +
  scale_color_gradient(
    limits = c(0, max_nes_both), low = "white", high = "black", name = "NES",
    guide = guide_colorbar(direction = "vertical", barheight = grid::unit(2, "cm"), barwidth = grid::unit(0.22, "cm"))
  ) +
  MrBiomics_void() +
  theme(legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0))

# combine with collected legends (only one colorbar from dummy, one size legend from Mono TA)
TF_combined_plot <- (HSC_TF_TA_plot | HSC_TF_EP_plot | Mono_TF_TA_plot | Mono_TF_EP_plot | legend_dummy)
TF_combined_plot <- TF_combined_plot +
  plot_layout(ncol = 5, widths = c(1, 1, 1, 1, 0.1), guides = "collect") &
  theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0))

ggsave_all_formats(path = TF_plot_path, plot = TF_combined_plot,
                   width = PLOT_SIZE_3_PER_ROW*2, height = PLOT_SIZE_3_PER_ROW)