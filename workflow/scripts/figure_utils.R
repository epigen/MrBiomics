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

# Handcrafted non-overlapping label placement
# all_feature_levels: character vector of all row levels in heatmap order (same as factor levels in heatmap_df[[feature_col]])
# selected_features: character vector of features to label
# box_size: minimal vertical separation between labels (in row-index units)
compute_nonoverlapping_label_positions <- function(all_feature_levels, selected_features, box_size = 1) {
    # Number of rows in heatmap
    n_rows <- length(all_feature_levels)
    if (length(selected_features) == 0) {
        return(data.frame(feature = character(0), orig_y = numeric(0), label_y = numeric(0)))
    }

    # Convert to top-down index (top row = 1)
    levels_top_to_bottom <- rev(all_feature_levels)
    orig_idx <- match(selected_features, levels_top_to_bottom)

    # Remove any NA (features not found in heatmap levels)
    valid <- !is.na(orig_idx)
    selected_features <- selected_features[valid]
    orig_idx <- orig_idx[valid]

    if (length(orig_idx) <= 1) {
        return(data.frame(feature = selected_features, orig_y = orig_idx, label_y = orig_idx))
    }

    # Sort by original position (top to bottom)
    ord <- order(orig_idx)
    orig_sorted <- orig_idx[ord]
    feat_sorted <- selected_features[ord]

    # Initialize label positions with originals
    pos <- as.numeric(orig_sorted)

    # Anchor extremes if they sit exactly at boundaries
    # Forward pass (top -> bottom): enforce minimal spacing
    for (i in seq(2, length(pos))) {
        pos[i] <- max(pos[i], pos[i - 1] + box_size)
    }
    # Clamp to bottom boundary
    pos <- pmin(pos, n_rows)

    # Backward pass (bottom -> top): pull up if needed while keeping spacing
    if (length(pos) >= 2) {
        for (i in seq(length(pos) - 1, 1)) {
            pos[i] <- min(pos[i], pos[i + 1] - box_size)
        }
    }

    # Enforce boundaries
    pos <- pmax(1, pmin(n_rows, pos))

    # If original is at exact top or bottom, keep it fixed
    # Adjust neighbors accordingly if necessary
    if (orig_sorted[1] == 1) {
        pos[1] <- 1
        if (length(pos) >= 2) {
            for (i in seq(2, length(pos))) {
                pos[i] <- max(pos[i], pos[i - 1] + box_size)
            }
            pos <- pmin(pos, n_rows)
        }
    }
    if (orig_sorted[length(orig_sorted)] == n_rows) {
        pos[length(pos)] <- n_rows
        if (length(pos) >= 2) {
            for (i in seq(length(pos) - 1, 1)) {
                pos[i] <- min(pos[i], pos[i + 1] - box_size)
            }
            pos <- pmax(1, pos)
        }
    }

    # Return in original selected_features order
    result_sorted <- data.frame(feature = feat_sorted, orig_y = orig_sorted, label_y = pos, stringsAsFactors = FALSE)
    # Reorder back to the original selected_features order
    reord <- match(selected_features, result_sorted$feature)
    result <- result_sorted[reord, , drop = FALSE]

    return(result)
}

plot_differential_features_heatmap <- function(dea_results_path, fig_path, fdr_threshold, log2FC_threshold, title = NULL, 
                                               feature, ct_clst_method, ct_clst_dist, feature_clst_method,
                                               feature_clst_dist, q_mask=0, label_box_size = 50, max_marker_labels = 25) {
    if (feature == 'Genes') {
        feature_col <- 'feature_name'
        y_label <- 'Differentially expressed genes'
    } else if (feature == 'Regions') {
        feature_col <- 'feature'
        y_label <- 'Differentially accessible regions'
    }
    
    # Get data for both up and down regulated features
    heatmap_df <- get_top_differential_features(dea_results_path, fdr_threshold, log2FC_threshold)

    # For Regions, avoid adding display suffixes; instead, use a simple unique key for pivoting
    # based on the (feature, feature_name) pair so that duplicate regions remain duplicated rows.
    if (feature == 'Regions') {
        heatmap_df <- heatmap_df %>%
            mutate(feature_id = paste0(feature, "||", feature_name))
        feature_col_to_plot <- 'feature_id'
    } else {
        feature_col_to_plot <- feature_col
    }

    # use hclust and dendsort to order the rows
    mat <- heatmap_df %>%
        select(group, all_of(feature_col_to_plot), logFC) %>%
        pivot_wider(names_from = group, values_from = logFC, values_fill = 0) %>%
        column_to_rownames(feature_col_to_plot) %>%
        as.matrix()
    row_order <- order.dendrogram(dendsort(as.dendrogram(
            hclust(
                dist(mat, method = feature_clst_dist),
                method = feature_clst_method)
            )))
    col_dendro <- dendsort(as.dendrogram(
            hclust(
                dist(t(mat), method = ct_clst_dist),
                method = ct_clst_method)
            ))
    col_order <- order.dendrogram(col_dendro)

    feature_vector <- heatmap_df[[feature_col_to_plot]]
    heatmap_df[[feature_col_to_plot]] <- factor(feature_vector, levels = row.names(mat)[row_order])
    heatmap_df$group <- factor(heatmap_df$group, levels = colnames(mat)[col_order])
    
    # Annotation bars for cell type and lineage
    annotation_df <- data.frame(group = colnames(mat)[col_order])
    annotation_df$group <- factor(annotation_df$group, levels = colnames(mat)[col_order])
    annotation_df$lineage <- CELL_TYPE_TO_LINEAGE_MAPPING[as.character(annotation_df$group)]
    
    cell_type_plot <- ggplot(annotation_df, aes(x = group, y = 1, fill = group)) +
        geom_tile() +
        scale_fill_manual(values = CELL_TYPE_COLORS, name = "Cell type", guide = guide_legend(nrow = 4)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_void() +
        theme(axis.title.y = element_blank(), plot.margin = margin(t = -12, r = 0, b = -8, l = 0, unit = "pt"))

    lineage_plot <- ggplot(annotation_df, aes(x = group, y = 1, fill = lineage)) +
        geom_tile() +
        scale_fill_manual(values = LINEAGE_COLORS, name = "Lineage", guide = guide_legend(nrow = 2)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_void() +
        theme(axis.title.y = element_blank(), plot.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))

    # Quantile masking to remove outliers from color scale of heatmap
    if(q_mask > 0) {
        upper_limit <- quantile(heatmap_df$logFC, probs = 1 - q_mask, na.rm=TRUE)
        lower_limit <- quantile(heatmap_df$logFC, probs = q_mask, na.rm=TRUE)
        heatmap_df$logFC <- ifelse(heatmap_df$logFC < lower_limit, lower_limit, heatmap_df$logFC)
        heatmap_df$logFC <- ifelse(heatmap_df$logFC > upper_limit, upper_limit, heatmap_df$logFC)
    }

    plot_limits <- c(-1, 1) * max(abs(heatmap_df$logFC), na.rm=TRUE)

    dendro_data <- dendro_data(col_dendro, type = "rectangle")
    dendro_plot <- ggdendrogram(dendro_data, labels=FALSE) + 
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_void() +
        theme(plot.margin = margin(t = 0, r = 0, b = -14, l = 0, unit = "pt"))

    heatmap_plot <- ggplot(heatmap_df, aes_string(x = "group", y = feature_col_to_plot, fill = "logFC")) +
        geom_tile(linewidth = 0) +
        scale_fill_distiller(palette = "RdBu", limits = plot_limits) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        # labs(x = NULL, y = y_label, title = title) +
        MrBiomics_theme() +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(t = -14, r = 0, b = -8, l = 0, unit = "pt")
        )

    # Prepare data for marker labels
    if (feature == 'Regions') {
        # Select region rows whose mapped gene is a known marker; keep duplicates
        marker_labels_df <- heatmap_df %>%
            filter(feature_name %in% HAEMATOPOIESIS_MARKERS)
        label_text_map <- setNames(marker_labels_df$feature_name, marker_labels_df[[feature_col_to_plot]])
    } else {
        marker_labels_df <- heatmap_df %>%
            filter(as.character(.data[[feature_col_to_plot]]) %in% HAEMATOPOIESIS_MARKERS) %>%
            distinct(.data[[feature_col_to_plot]], .keep_all = TRUE)
        # For genes, label text is the gene itself
        label_text_map <- setNames(as.character(marker_labels_df[[feature_col_to_plot]]), as.character(marker_labels_df[[feature_col_to_plot]]))
    }

    # Create a separate plot for the marker labels using handcrafted layout
    if (nrow(marker_labels_df) > 0) {
        all_levels <- levels(heatmap_df[[feature_col_to_plot]])
        selected_features <- as.character(marker_labels_df[[feature_col_to_plot]])

        # Limit the number of marker labels to at most max_marker_labels.
        # If more exist, choose a subset that is spaced evenly along the full heatmap height
        # by selecting the nearest available marker to evenly spaced target rows.
        if (length(selected_features) > max_marker_labels) {
            feature_indices <- match(selected_features, all_levels)
            # Remove any NA indices defensively
            valid <- !is.na(feature_indices)
            selected_features <- selected_features[valid]
            feature_indices <- feature_indices[valid]

            # Targets are evenly spaced along all rows, not just among existing markers
            n_rows <- length(all_levels)
            target_rows <- unique(round(seq(1, n_rows, length.out = max_marker_labels)))

            # Greedy nearest-neighbor selection without replacement
            keep <- logical(length(selected_features))
            taken <- rep(FALSE, length(selected_features))
            for (t in target_rows) {
                # compute distances to target for all available markers not yet taken
                available <- which(!taken)
                if (length(available) == 0) break
                dists <- abs(feature_indices[available] - t)
                pick_local <- available[which.min(dists)]
                keep[pick_local] <- TRUE
                taken[pick_local] <- TRUE
            }
            # Fallback in case unique targets < max (e.g., very few markers spread),
            # fill remaining slots by picking closest remaining markers in order top-to-bottom
            if (sum(keep) < min(max_marker_labels, length(selected_features))) {
                remaining <- which(!taken)
                # order remaining by feature index top-to-bottom
                ord_rem <- order(feature_indices[remaining])
                fill_n <- min(min(max_marker_labels, length(selected_features)) - sum(keep), length(remaining))
                keep[remaining[ord_rem][seq_len(fill_n)]] <- TRUE
            }
            selected_features <- selected_features[keep]
        }

        pos_df <- compute_nonoverlapping_label_positions(
            all_feature_levels = all_levels,
            selected_features = selected_features,
            box_size = label_box_size
        )

        n_rows <- length(all_levels)
        # Precompute plotting-space y coordinates (bottom = 1, top = n_rows)
        pos_df$label_y_plot <- n_rows + 1 - pos_df$label_y
        pos_df$orig_y_plot  <- n_rows + 1 - pos_df$orig_y
        # Attach label text (gene symbols for regions, feature itself for genes)
        pos_df$label <- as.character(label_text_map[pos_df$feature])

        marker_label_plot <- ggplot(pos_df) +
            # Connect calculated label positions to original row positions
            geom_segment(aes(x = 0.85, xend = 1.0, y = label_y_plot, yend = orig_y_plot),
                         size = 0.25, color = "grey40") +
            # Draw labels
            geom_text(aes(x = 0.8, y = label_y_plot, label = label),
                      size = 2, hjust = 1) +
            # Fixed limits matching number of rows; flush to edges
            scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = c(0, 0)) +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            theme_void() +
            coord_cartesian(clip = "off") +
            theme(plot.margin = margin(0, 0, 0, 0))
    } else {
        marker_label_plot <- plot_spacer()
    }
    
    # Assemble the final plot using patchwork's wrap_plots for robust layout
    plot_list <- list(
        plot_spacer(),       dendro_plot,
        marker_label_plot,   heatmap_plot,
        plot_spacer(),       cell_type_plot,
        plot_spacer(),       lineage_plot
    )

    gp <- wrap_plots(plot_list, ncol = 2,
                     widths = c(1, 4),
                     heights = c(0.6, 10, 0.35, 0.35)) +
        plot_layout(guides = "collect") +
        plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) &
        theme(legend.position = "bottom", legend.box="vertical", plot.margin = margin(0, 0, 0, 0))

    width <- 5 # Increased width to accommodate labels
    height <- 7
    ggsave_all_formats(path = fig_path,
                       plot = gp,
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