#### UMAP plotting ####
umap_plot <- function(data_path, fig_path, title = NULL) {
    # Load UMAP data
    umap_data <- data.frame(fread(file.path(data_path), header=TRUE))
    
    # Parse sample names to extract cell type
    umap_data$cell_type <- sapply(strsplit(umap_data$sample_name, "_"), function(x) x[1])
    
    # Translate cell type names to match CELL_TYPE_COLORS
    umap_data$cell_type_colored <- DATA_TO_CELL_TYPE_COLORS_MAPPING[umap_data$cell_type]
    
    # Order cell types according to CELL_TYPE_COLORS
    umap_data$cell_type_colored <- factor(umap_data$cell_type_colored, 
                                         levels = names(CELL_TYPE_COLORS))
    
    # Create scatterplot
    umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type_colored)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = CELL_TYPE_COLORS, name = "Cell type") +
        labs(x = "UMAP 1", y = "UMAP 2", title = title) +
        theme_minimal() +
        MrBiomics_theme() +
        theme(
            legend.position = "right",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()
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


#### DEA HEATMAP ####
get_top_differential_features <- function(dea_results_path, fdr_threshold, log2FC_threshold) {
    # read data
    dea_df <- data.frame(fread(file.path(dea_results_path), header=TRUE))

    # get the features that are significant and have a logFC above the threshold
    significant_features <- dea_df %>%
        filter(adj.P.Val < fdr_threshold, abs(logFC) > log2FC_threshold) %>%
        pull(feature) %>%
        unique()
    
    # filter the original df for these features
    dea_top_features_df <- dea_df %>%
        filter(feature %in% significant_features)
    
    # rename groups to nicer names
    dea_top_features_df <- dea_top_features_df %>%
        mutate(group = DATA_TO_CELL_TYPE_COLORS_MAPPING[group])

    return(dea_top_features_df)
}

plot_differential_features_heatmap <- function(dea_results_path, fig_path, fdr_threshold, log2FC_threshold, title = NULL, 
                                               feature, ct_clst_method, ct_clst_dist, q_mask=0) {
    if (feature == 'Genes') {
        feature_col <- 'feature_name'
        y_label <- 'Differentially expressed genes'
    } else if (feature == 'Regions') {
        feature_col <- 'feature'
        y_label <- 'Differentially accessible regions'
    }
    
    # Get data for both up and down regulated features
    heatmap_df <- get_top_differential_features(dea_results_path, fdr_threshold, log2FC_threshold)

    # use hclust and dendsort to order the rows
    mat <- heatmap_df %>%
        select(group, all_of(feature_col), logFC) %>%
        pivot_wider(names_from = group, values_from = logFC, values_fill = 0) %>%
        column_to_rownames(feature_col) %>%
        as.matrix()
    order_dendro <- order.dendrogram(dendsort(as.dendrogram(hclust(dist(mat)))))
    mat <- mat[order_dendro, ]
    
    # Quantile masking
    if(q_mask > 0) {
        upper_limit <- quantile(heatmap_df$logFC, probs = 1 - q_mask, na.rm=TRUE)
        lower_limit <- quantile(heatmap_df$logFC, probs = q_mask, na.rm=TRUE)
        heatmap_df$logFC <- ifelse(heatmap_df$logFC < lower_limit, lower_limit, heatmap_df$logFC)
        heatmap_df$logFC <- ifelse(heatmap_df$logFC > upper_limit, upper_limit, heatmap_df$logFC)
    }

    plot_limits <- c(-1, 1) * max(abs(heatmap_df$logFC), na.rm=TRUE)

    dendsort_wrapper <- function(hc, ...) {
        dendsort::dendsort(as.dendrogram(hc))
    }

    heatmap_plot <- ggheatmap(mat) + 
        scale_fill_distiller(palette = "RdBu", limits = plot_limits) +
        # labs(x = NULL, y = y_label, title = title) +
        MrBiomics_theme() + 
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.text.y = element_blank(),
        ) + 
        anno_top(size=0.1) + 
        align_dendro(reorder_dendrogram=dendsort_wrapper, 
                     method = ct_clst_method, 
                     distance = ct_clst_dist) + 
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank()
        )


    width <- 4
    height <- 5 # Adjusted for two plots stacked.
    ggsave_all_formats(path = fig_path,
                       plot = heatmap_plot,
                       width = width,
                       height = height)
    
    return(heatmap_plot)
}


#### ENRICHMENT HEATMAP ####
filter_top_terms <- function(df, fdr_threshold) {
    df_sig <- df %>%
      filter(score > 0, statistic < fdr_threshold) 

    tissues_to_keep <- c("PBMC", "Bone Marrow")
    df_sig <- df_sig %>%
      filter(grepl(paste(tissues_to_keep, collapse = "|"), Term))

    # Identify top term for each cell type
    top_hits <- df_sig %>%
      group_by(name) %>%
      slice_max(order_by = score, n = 2, with_ties = FALSE) %>%
      ungroup()

    # Get the list of top terms
    top_term_strings <- top_hits %>%
      pull(Term) %>%
      unique()

    # Create a mapping from term to the cell type for which it's a top hit
    term_to_top_name_map <- top_hits %>%
      select(Term, top_for_name = name)

    # Filter original dataframe for top terms and add the 'top_for_name' info
    df_top <- df %>%
      filter(Term %in% top_term_strings)

    # make the terms a factor that is ordered by the order of the matching cell types in CELL_TYPE_COLORS
    term_to_top_name_map$top_for_name <- factor(term_to_top_name_map$top_for_name, levels = names(CELL_TYPE_COLORS))
    term_to_top_name_map <- term_to_top_name_map %>% arrange(top_for_name)
    df_top$Term <- factor(df_top$Term, levels = unique(term_to_top_name_map$Term))
    
    return(df_top)
}

prepare_for_heatmap <- function(df_formatted, fdr_threshold) {
    df_top <- filter_top_terms(df_formatted, fdr_threshold)

    mat_df <- df_top %>%
      select(name, Term, score) %>%
      pivot_wider(names_from = name, values_from = score, values_fill = 0)

    mat <- mat_df %>% column_to_rownames("Term") %>% as.matrix()

    # Shorten term names for plot readability
    original_term_levels <- levels(df_top$Term)
    short_term_levels <- original_term_levels

    # Abbreviate tissue and remove level information
    short_term_levels <- gsub("Bone Marrow", "BM", short_term_levels)
    # short_term_levels <- gsub("-L[1-9]-", "-", short_term_levels)

    # Define and apply cell type abbreviations
    cell_type_replacements <- c(
        "Hematopoeitic Stem Cell" = "HSC",
        "Lymphoid Primed Multipotent Progenitor" = "LMPP",
        "Hematopoietic Stem And Progenitor Cell" = "HSPC",
        "Basophil Eosinophil Mast Progenitor" = "BEMP",
        "Early Erythroid" = "Pro Ery",
        "Common Lymphoid Progenitor" = "CLP",
        "Granulocyte Monocyte Progenitor" = "GMP",
        "Erythroid Cell" = "Ery",
        "CD14 Monocyte" = "Mono",
        "Intermediate B Cell, Kappa Light Chain" = "Int B κ+",
        "CD8 Memory" = "CD8 mem",
        "Mucosal Associated Invariant T" = "MAIT",
        "Intermediate B Cell" = "Int B",
        "CD56-dim Natural Killer" = "CD56dim NK",
        "Natural Killer" = "NK",
        "Monocyte" = "Mono",
        "CD4 T" = "CD4",
        "Naive B" = "Naive B",
        "B" = "B",
        "Progenitor B" = "Pro B",
        "Naive B Cell, Kappa Light Chain" = "Naive B κ+",
        "Precursor Plasmacytoid Dendritic Cell" = "pDC",
        "Conventional Dendritic Cell 1" = "cDC1",
        "Transitional B" = "Transitional B",
        "Platelet" = "Platelet",
        "gamma-delta T 2" = "Vδ2 T cell",
        "CD16 Monocyte" = "CD16 Mono",
        "Innate Lymphoid Cell" = "ILC",
        "Hematopoeitic Stem And Progenitor Cell" = "HSPC",
        "Erythroid Megakaryocyte Progenitor" = "MEP"
    )
    
    # Sort by length to avoid partial matches on shorter strings
    cell_type_replacements <- cell_type_replacements[order(nchar(names(cell_type_replacements)), decreasing = TRUE)]
    for (i in seq_along(cell_type_replacements)) {
        short_term_levels <- gsub(names(cell_type_replacements)[i], cell_type_replacements[i], short_term_levels)
    }

    # Print the mapping of original to shortened names
    term_map_df <- data.frame(Original = original_term_levels, Shortened = short_term_levels)
    message("Shortening term names for heatmap:")
    print(term_map_df)
    
    # Create a named vector for applying the new names
    term_map <- setNames(short_term_levels, original_term_levels)

    heatmap_df <- mat_df %>%
      pivot_longer(-Term, names_to = "name", values_to = "score") %>%
      left_join(df_formatted %>% select(name, Term, statistic), by = c("name", "Term")) %>%
      mutate(
        Term = factor(term_map[Term], levels = rev(unique(short_term_levels))),
        name = factor(name, levels = names(CELL_TYPE_COLORS)),
        sig = (statistic < fdr_threshold) & (!is.na(statistic)) & (!is.infinite(statistic)),
        neg_log10_statistic = ifelse(is.infinite(-log10(statistic)),
                                     max(
                                      -log10(statistic[!is.infinite(-log10(statistic)) & !is.na(statistic)]),
                                      na.rm = TRUE
                                      )*1.1,
                                     -log10(statistic))
      )

    return(heatmap_df)
}

plot_enrichment_heatmap <- function(heatmap_df, fig_path, fill_lab, size_lab, title = NULL, ylabel = NULL) {
    # mask those values where the log statistic is tiny with NaN, to drop them in the heatmap
    heatmap_df$neg_log10_statistic <- ifelse(heatmap_df$neg_log10_statistic < 1, NaN, heatmap_df$neg_log10_statistic)

    enrichment_plot <- ggplot(heatmap_df, aes(x = name, y = Term, size=neg_log10_statistic, fill = score)) +
      geom_point(shape=21, stroke=0.25) +
      # add star for significance
      geom_text(aes(label = ifelse(sig, "✳︎", "")), vjust = 0.5, size=3, color = "white") +
      scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)*max(abs(heatmap_df$score)), name = fill_lab) +
      scale_size_continuous(name = size_lab) +
      labs(title = title, x = "Cell type", y = ylabel) +
      # ensure square tiles
      coord_fixed() +
      MrBiomics_theme() + 
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    # Save plot
    width <- 7
    height <- 5
    ggsave_all_formats(path = fig_path,
                       plot = enrichment_plot,
                       width = width,
                       height = height)
    
    return(enrichment_plot)
}