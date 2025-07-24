
#### libraries ####
library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")

#### configuration ####

# input
data_path <- snakemake@input[["data"]]
metadata_path <- snakemake@input[["metadata"]]
dea_results_path <- snakemake@input[["dea_results"]]

# output
plot_paths <- snakemake@output[["correlation_plots"]]

# parameters
adjp_th <- snakemake@params[["adjp_th"]] #0.05
lfc_th <- snakemake@params[["lfc_th"]] #2
ave_expr_th <- snakemake@params[["ave_expr_th"]] #1

# plotting
width <- 3
height <- 3
# options(repr.plot.width=width, repr.plot.height=height)

# load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)
dea_results <- data.frame(fread(file.path(dea_results_path), header=TRUE))

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

    # split data for plotting
    df_conv <-  df %>% filter(category == "Convergent")
    df_div  <-  df %>% filter(category != "Convergent")   

    # prepare annotations instead of legend
    cat_cols <- c("Convergent"             = "grey30",
              "Epigenetic potential"   = "red",
              "Transcriptional abundance"= "deepskyblue")

    # annotation labels
    ann_df <- data.frame(
        category = c("Transcriptional abundance", "Convergent", "Epigenetic potential"),
        x = c(-Inf,  -Inf,  Inf),   # TL, BL, BR
        y = c( Inf,  -Inf, -Inf),
        hjust = c(0, 0, 1),        # keep text inside panel
        vjust = c(1, 0, 0),
        lab = c(
            paste0("Transcriptional\nabundance\n (", counts["Transcriptional abundance"],")"),
            paste0("Convergent\n (", counts["Convergent"],")"),
            paste0("Epigenetic\npotential\n (", counts["Epigenetic potential"],")")
            )
        )
    
    ## ---------- plot ----------
    p <- ggplot(df, aes(x = ATAC, y = RNA)) +
    # plot convergent first (small, faint – background)
    geom_point(data = df_conv,
             aes(colour = category),
             size  = 0.25,
             alpha = 0.05) +
    # overlay divergent classes (larger, more opaque – foreground)
    geom_point(data = df_div,
             aes(colour = category),
             size  = 0.25,
             alpha = 0.25) +
    # colour palette 
    scale_colour_manual(
    values = cat_cols,
    name   = NULL
    ) +
    # thin density contour on top
    stat_density_2d(
    data  = df,
    aes(x = ATAC, y = RNA,
        colour = "grey", group = category),
    contour_var = "ndensity",
    bins  = 5,
    linewidth  = 0.3, alpha = 1, show.legend = FALSE
  )+
    # annotation boxes with counts, also serving as legend
    geom_label(
    data = ann_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = lab, fill = category, hjust = hjust, vjust = vjust),
    colour = "black",
    size = 3,
    label.size = 0.3,
    lineheight = 0.85,
    label.padding = grid::unit(0.4, "lines"),
    label.r = grid::unit(0.1, "lines")
  ) +
  scale_fill_manual(values = scales::alpha(cat_cols, 0.25), guide = "none")+
    # axis titles
    labs(
    title = paste0(ct,"\n",sprintf("Pearson's R = %.2f", r_val)),
    x = "Chromatin accessibility\n(signal intensity, normalized & integrated)",
    y = "Gene expression\n(signal intensity, normalized & integrated)"
    ) +
    theme_minimal() +
    theme(plot.title      = element_text(size = 8, hjust = 0.5),
          axis.title.x    = element_text(size = 8),
    axis.title.y    = element_text(size = 8),
        legend.position = "none",  # remove legend
          )

    # print(p)
    
  ## ---------- save plot ----------
  ggsave(file.path(dirname(plot_paths[[1]]),paste0(ct,"_correlation.png")), p, width = width, height = height, dpi = 300, bg="white")
}
