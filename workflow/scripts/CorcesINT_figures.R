source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# FIXME snakemakeify
#### configuration ####
# input
data_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated.csv"
metadata_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/spilterlize_integrate/all/annotation.csv"
dea_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesINT/dea_limma/normupperquartile_integrated/results.csv"
gene_annotation_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/rnaseq_pipeline/counts/gene_annotation.csv"
# parameters
adjp_th <- 0.05
lfc_th <- 1
ave_expr_th <- 0
# output
plot_dir <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesINT/correlation_plots"

# load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)
dea_results <- data.frame(fread(file.path(dea_results_path), header=TRUE))
gene_annotation <- data.frame(fread(file.path(gene_annotation_path), header=TRUE))
head(gene_annotation)

# transform data for plotting (for each group)
dea <- dea_results |> mutate(cell_type = sub("cell_type(.*?)__.*", "\\1", group))

#### iterate over groups ####
for(ct in unique(metadata$cell_type)){

    ## ---------- summarise mean signal per gene ----------
    atac <- rownames(metadata)[metadata$cell_type == ct & metadata$modality == "ATAC"]
    rna  <- rownames(metadata)[metadata$cell_type == ct & metadata$modality == "RNA" ]
    
    df <- data.frame(
    feature = rownames(data),
    ATAC = rowMeans(data[ , atac, drop = FALSE]),
    RNA  = rowMeans(data[ , rna , drop = FALSE])
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

    # split data for plotting
    df_conv <-  df %>% filter(category == "Convergent")
    df_div  <-  df %>% filter(category != "Convergent")

    # prepare annotations instead of legend
    cat_cols <- c("Convergent"             = "grey80",
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
    geom_point(data = df_conv,
             aes(colour = category),
             size  = 0.5,
             alpha = 0.1,
             stroke = 0,
             show.legend = FALSE) +
    # overlay divergent classes (larger, more opaque – foreground)
    geom_point(data = df_div %>% filter(markers == ""),
             aes(colour = category, size = log10_adjp),
             alpha = 0.5,
             stroke = 0) +
    # colour palette 
    scale_colour_manual(
    values = cat_cols,
    name   = NULL, 
    guide = "none"
    ) +
    # size legend using original -log10(adjusted p) values
    scale_size_continuous(
      name = "-log10(p-adj.)",
      range = c(0.5, 4),
      breaks = size_breaks,
      labels = function(x) format(x, digits = 2)
    ) +
    # annotation boxes with counts, also serving as legend
    geom_text(
    data = ann_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = lab, color = category, hjust = hjust, vjust = vjust),
    size = 5,
    ) +
    # annotate divergent HAEMATOPOIESIS_MARKERS
    geom_text_repel(
        color = "black",
        size = 4,
        box.padding = 0.5,
        max.overlaps = Inf,
        seed = 42,
        show.legend = FALSE
    ) +
    # redraw labeled points on top so labels don't occlude points
    geom_point(
        data = df %>% filter(markers != ""),
        aes(colour = category, size = log10_adjp),
        alpha = 0.5,
        stroke = 0
    ) +
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
  ggsave(file.path(plot_dir,paste0(ct,"_correlation.png")), p, width = PLOT_HEIGHT+1, height = PLOT_HEIGHT,
         dpi = 300, bg="white")
}
