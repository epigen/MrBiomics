#### UMAP plotting ####
umap_plot <- function(data_path, fig_path) {
    # Load UMAP data
    umap_data <- data.frame(fread(file.path(data_path), header=TRUE))
    
    # Parse sample names to extract cell type
    umap_data$cell_type <- sapply(strsplit(umap_data$sample_name, "_"), function(x) x[1])
    
    # Translate cell type names to match celltype_colors
    umap_data$cell_type_colored <- data_to_colors_mapping[umap_data$cell_type]
    
    # Order cell types according to celltype_colors
    umap_data$cell_type_colored <- factor(umap_data$cell_type_colored, 
                                         levels = names(celltype_colors))
    
    # Create scatterplot
    umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type_colored)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = celltype_colors, name = "Cell type") +
        labs(x = "UMAP 1", y = "UMAP 2") +
        theme_minimal() +
        MrBiomics_theme() +
        theme(
            legend.position = "right",
            panel.grid.minor = element_blank()
        )
    
    # Save plot
    width <- 6
    height <- 5
    ggsave_all_formats(path = fig_path,
                       plot = umap_plot,
                       width = width,
                       height = height)
    
    return(umap_plot)
}

#### ENRICHMENT HEATMAP ####
filter_top_terms <- function(df, fdr_threshold) {
    df_sig <- df %>%
      filter(score > 0, statistic < fdr_threshold)

    top_terms <- df_sig %>%
      group_by(name) %>%
      slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
      pull(Term)

    df_top <- df %>%
      filter(Term %in% top_terms)

    return(df_top)
}

reshape_for_heatmap <- function(df_top, df_formatted, fdr_threshold) {
    mat_df <- df_top %>%
      select(name, Term, score) %>%
      pivot_wider(names_from = name, values_from = score, values_fill = 0)

    mat <- mat_df %>% column_to_rownames("Term") %>% as.matrix()

    row_dendro <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
    row_order <- order.dendrogram(row_dendro)

    heatmap_df <- mat_df %>%
      pivot_longer(-Term, names_to = "name", values_to = "score") %>%
      left_join(df_formatted %>% select(name, Term, statistic), by = c("name", "Term")) %>%
      mutate(
        Term = factor(Term, levels = rownames(mat)[row_order]),
        name = factor(name, levels = names(celltype_colors)),
        sig = (statistic < fdr_threshold) & (!is.na(statistic)) & (!is.infinite(statistic)),
        neg_log10_statistic = ifelse(is.infinite(-log10(statistic)),
                                     max(
                                      -log10(statistic[!is.infinite(-log10(statistic)) & !is.na(statistic)]), na.rm = TRUE
                                      )*1.1,
                                     -log10(statistic))
      )

    return(heatmap_df)
}

plot_enrichment_heatmap <- function(heatmap_df, fig_path, fill_lab, size_lab) {
    # make plot
    enrichment_plot <- ggplot(heatmap_df, aes(x = name, y = Term, size=neg_log10_statistic, fill = score)) +
      geom_point(shape=21, stroke=0.25) +
      # add star for significance
      geom_text(aes(label = ifelse(sig, "✳︎", "")), vjust = 0.5, size=3, color = "white") +
      scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)*max(abs(heatmap_df$score)), name = fill_lab) +
      scale_size_continuous(name = size_lab) +
      # ensure square tiles
      coord_fixed() +
      MrBiomics_theme() + 
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_blank()
      )

    # Save plot
    width <- 10
    height <- 5
    ggsave_all_formats(path = fig_path,
                       plot = enrichment_plot,
                       width = width,
                       height = height)
    
    return(enrichment_plot)
}