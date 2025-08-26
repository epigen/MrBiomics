source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# FIXME snakemakeify
#### configuration ####
# input
unintegrated_umap_coords_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/unsupervised_analysis/normupperquartile/UMAP/UMAP_correlation_15_0.1_2_data.csv"
integrated_umap_coords_path <-   "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/unsupervised_analysis/normupperquartile_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"
unintegrated_cfa_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/spilterlize_integrate/all/normupperquartile_CFA.csv"
integrated_cfa_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated_CFA.csv"
norm_counts_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated.csv"
metadata_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/spilterlize_integrate/all/annotation.csv"
dea_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/dea_limma/normupperquartile_integrated/results.csv"
gene_annotation_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/rnaseq_pipeline/counts/gene_annotation.csv"
enrichment_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/GO_Biological_Process_2025/cell_types_GO_Biological_Process_2025_all.csv"
# enrichment_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/enrichment_analysis/cell_types/preranked_GSEApy/ReactomePathways/cell_types_ReactomePathways_all.csv"

# parameters
adjp_th <- 0.05
fdr_threshold <- 0.05
lfc_th <- 1
ave_expr_th <- 0
# output
unintegrated_cfa_plot_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/unintegrated_cfa.pdf"
integrated_cfa_plot_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/integrated_cfa.pdf"
integrated_umap_plot_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/integrated_umap.pdf"
unintegrated_umap_plot_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/unintegrated_umap.pdf"
epigenetic_scatter_dir <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/correlation_plots"
int_enrichment_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/enrichment.pdf"

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
plot_cfa_heatmap <- function(cfa_mat, title, path, var_max=NULL, nPCs=10, metadata_rows=c('cell_type', 'modality', 'donor')){
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
        axis.text.y = element_blank(),integrated_cfa.png
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
    geom_text(aes(label = round(stat, 1), color = text_color, fontface = text_fontface), size = 3) +
    scale_fill_gradient(low = "white", high = as.character(RdBu_extremes["up"]), name='-log10(p-adj.)') +
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

int_df_formatted <- create_int_enrichment_df(enrichment_results_path, fdr_threshold)
int_heatmap_df <- prepare_for_heatmap(df_formatted = int_df_formatted, fdr_threshold = fdr_threshold, top_n_per_name = 2)
int_enrichment_plot <- plot_clustered_enrichment_heatmap(
    heatmap_df = int_heatmap_df,
    fig_path = int_enrichment_path,
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