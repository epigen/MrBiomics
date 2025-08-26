########################################################################################################################
#### UMAP plotting #####################################################################################################
########################################################################################################################
umap_plot <- function(data_path, fig_path, title = NULL, modality_by_shape = FALSE) {
    # Load UMAP data
    umap_data <- data.frame(fread(file.path(data_path), header=TRUE))
    
    # Parse sample names to extract cell type
    umap_data$cell_type <- sapply(strsplit(umap_data$sample_name, "_"), function(x) x[1])
    
    # Translate cell type names to match CELL_TYPE_COLORS
    umap_data$cell_type_colored <- DATA_TO_CELL_TYPE_COLORS_MAPPING[umap_data$cell_type]
    
    # Order cell types according to CELL_TYPE_COLORS
    umap_data$cell_type_colored <- factor(umap_data$cell_type_colored, 
                                         levels = names(CELL_TYPE_COLORS))

    if (modality_by_shape) {
        umap_data$modality <- sapply(strsplit(umap_data$sample_name, "_"), function(x) x[3])
        umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type_colored, shape = modality)) +
            scale_shape_manual(values = c("ATAC" = 17, "RNA" = 16), name = NULL)
    } else {
        umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type_colored))
    }

    umap_plot <- umap_plot +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(values = CELL_TYPE_COLORS, name = "Cell type") +
        labs(x = "UMAP 1", y = "UMAP 2", title = title) +
        MrBiomics_theme() +
        theme(
            legend.position = "right",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            aspect.ratio = 1
        )
    
    # Save plot
    ggsave_all_formats(path = fig_path,
                       plot = umap_plot,
                       width = PLOT_SIZE_3_PER_ROW+1,
                       height = PLOT_SIZE_3_PER_ROW)
    
    return(umap_plot)
}


########################################################################################################################
#### DEA HEATMAP #######################################################################################################
########################################################################################################################
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
# heatmap_rows: character vector of all row levels in heatmap order (same as factor levels in heatmap_df[[feature_col]])
# selected_features: character vector of features to label
# box_size: minimal vertical separation between labels (in row-index units)
compute_nonoverlapping_label_positions <- function(heatmap_rows, selected_features, box_size = 1) {
	# Number of rows in heatmap
	n_rows <- length(heatmap_rows)
	if (length(selected_features) == 0) {
		return(data.frame(feature = character(0), orig_y = numeric(0), label_y = numeric(0)))
	}

	# Convert to top-down index (top row = 1)
	levels_top_to_bottom <- rev(heatmap_rows)
	orig_idx <- match(selected_features, levels_top_to_bottom)

	# Remove any NA (features not found in heatmap levels)
	valid <- !is.na(orig_idx)
	selected_features <- selected_features[valid]
	orig_idx <- orig_idx[valid]

	# Define safe plotting bounds to avoid touching the frame
	lower_bound <- 1
	upper_bound <- n_rows - box_size/2
	if (upper_bound < lower_bound) {
		# Degenerate case: too large margin relative to number of rows
		lower_bound <- 1
		upper_bound <- n_rows
	}

	if (length(orig_idx) <= 1) {
		pos_single <- pmin(pmax(orig_idx, lower_bound), upper_bound)
		return(data.frame(feature = selected_features, orig_y = orig_idx, label_y = pos_single))
	}

	# Sort by original position (top to bottom)
	ord <- order(orig_idx)
	orig_sorted <- orig_idx[ord]
	feat_sorted <- selected_features[ord]

	# Initialize label positions with originals
	pos <- as.numeric(orig_sorted)

	# Forward pass (top -> bottom): enforce minimal spacing
	for (i in seq(2, length(pos))) {
		pos[i] <- max(pos[i], pos[i - 1] + box_size)
	}
	# Clamp to safe bottom boundary
	pos <- pmin(pos, upper_bound)

	# Backward pass (bottom -> top): pull up if needed while keeping spacing
	if (length(pos) >= 2) {
		for (i in seq(length(pos) - 1, 1)) {
			pos[i] <- min(pos[i], pos[i + 1] - box_size)
		}
	}

	# Enforce safe boundaries
	pos <- pmax(lower_bound, pmin(upper_bound, pos))

	# If original is at or beyond top safe bound, anchor there and push others
	if (orig_sorted[1] <= lower_bound) {
		pos[1] <- lower_bound
		if (length(pos) >= 2) {
			for (i in seq(2, length(pos))) {
				pos[i] <- max(pos[i], pos[i - 1] + box_size)
			}
			pos <- pmin(pos, upper_bound)
		}
	}
	# If original is at or beyond bottom safe bound, anchor there and pull others
	if (orig_sorted[length(orig_sorted)] >= upper_bound) {
		pos[length(pos)] <- upper_bound
		if (length(pos) >= 2) {
			for (i in seq(length(pos) - 1, 1)) {
				pos[i] <- min(pos[i], pos[i + 1] - box_size)
			}
			pos <- pmax(lower_bound, pos)
		}
	}

	# Return in original selected_features order
	result_sorted <- data.frame(feature = feat_sorted, orig_y = orig_sorted, label_y = pos, stringsAsFactors = FALSE)
	# Reorder back to the original selected_features order
	reord <- match(selected_features, result_sorted$feature)
	result <- result_sorted[reord, , drop = FALSE]

	return(result)
}

compute_diagonal_cluster_order <- function(mat,
                                           row_clst_dist = "euclidean",
                                           row_clst_method = "ward.D2",
                                           col_clst_dist = "euclidean",
                                           col_clst_method = "ward.D2",
                                           n_clusters = 25) {

    # Row clustering (features)
    row_hclust <- hclust(dist(mat, method = row_clst_dist), method = row_clst_method)
    row_dendro <- dendsort(as.dendrogram(row_hclust))
    row_order <- order.dendrogram(row_dendro)

    # Column clustering (cell types)
    col_hclust <- hclust(dist(t(mat), method = col_clst_dist), method = col_clst_method)
    col_dendro <- dendsort(as.dendrogram(col_hclust))
    col_order <- order.dendrogram(col_dendro)

    # Derive diagonal-friendly row order by clustering rows, assigning each cluster to
    # the column with the highest mean, then grouping clusters by column order
    rownames_in_hclust_order <- rownames(mat)[row_order]
    k_effective <- max(1, min(n_clusters, nrow(mat)))
    cluster_assignments <- cutree(row_hclust, k = k_effective)
    cluster_ids <- sort(unique(cluster_assignments))
    column_names_in_order <- colnames(mat)[col_order]

    cluster_summary <- lapply(cluster_ids, function(cl_id) {
        rows_in_cluster <- names(cluster_assignments)[cluster_assignments == cl_id]
        col_means <- colMeans(mat[rows_in_cluster, , drop = FALSE])
        best_col <- names(which.max(col_means))
        first_pos <- min(match(rows_in_cluster, rownames_in_hclust_order), na.rm = TRUE)
        list(
            cluster_id = cl_id,
            rows = rows_in_cluster,
            best_column = best_col,
            best_column_rank = match(best_col, column_names_in_order),
            first_pos = first_pos
        )
    })

    cluster_order_idx <- order(
        vapply(cluster_summary, function(x) x$best_column_rank, numeric(1)),
        vapply(cluster_summary, function(x) x$first_pos, numeric(1))
    )
    cluster_summary <- cluster_summary[cluster_order_idx]

    final_row_order <- unlist(lapply(cluster_summary, function(cs) {
        rows <- cs$rows
        rows[order(match(rows, rownames_in_hclust_order))]
    }), use.names = FALSE)

    return(list(
        row_order = final_row_order,
        col_order = column_names_in_order,
        col_dendro = col_dendro
    ))
}

build_column_annotations <- function(col_order, col_dendro, title = NULL) {
    annotation_df <- data.frame(group = col_order)
    annotation_df$group <- factor(annotation_df$group, levels = col_order)
    annotation_df$lineage <- CELL_TYPE_TO_LINEAGE_MAPPING[as.character(annotation_df$group)]

    cell_type_plot <- ggplot(annotation_df, aes(x = group, y = 1, fill = group)) +
        geom_tile(color = 'white', linewidth = 0.5) +
        scale_fill_manual(values = CELL_TYPE_COLORS, name = "Cell type", guide = "none") +
        scale_x_discrete(expand = expansion(add = 0.6)) +
        scale_y_continuous(expand = c(0, 0)) +
        MrBiomics_void() +
        coord_cartesian(clip = "off") +
        theme(axis.title.y = element_blank(), plot.margin = margin(t = -12, r = 0, b = -8, l = 0, unit = "pt"))

    lineage_plot <- ggplot(annotation_df, aes(x = group, y = 1, fill = lineage)) +
        geom_tile(color = 'white', linewidth = 0.5) +
        scale_fill_manual(
            values = LINEAGE_COLORS,
            name = "Lineage",
            guide = guide_legend(ncol = 1, keyheight = grid::unit(0.4, "lines"), keywidth = grid::unit(0.6, "lines"))
        ) +
        scale_x_discrete(expand = expansion(add = 0.6)) +
        scale_y_continuous(expand = c(0, 0)) +
        MrBiomics_void() +
        coord_cartesian(clip = "off") +
        theme(axis.title.y = element_blank(), plot.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"))

    cell_type_label_plot <- ggplot() +
        annotate("text", x = 1, y = 0.5, label = "Cell type", hjust = 1, vjust = 0.5) +
        xlim(0, 1) + ylim(0, 1) +
        MrBiomics_void() +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(t = -12, r = 0, b = -8, l = 4, unit = "pt"))

    lineage_label_plot <- ggplot() +
        annotate("text", x = 1, y = 0.5, label = "Lineage", hjust = 1, vjust = 0.5) +
        xlim(0, 1) + ylim(0, 1) +
        MrBiomics_void() +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(t = -10, r = 0, b = 0, l = 4, unit = "pt"))

    n_cols <- length(col_order)
    dendro_data_obj <- dendro_data(col_dendro, type = "rectangle")
    dendro_plot <- ggdendrogram(dendro_data_obj, labels = FALSE) +
        scale_x_continuous(limits = c(0.5, n_cols + 0.5), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = title) +
        MrBiomics_void() +
        theme(plot.margin = margin(t = 0, r = 0, b = -14, l = 0, unit = "pt"))

    return(list(
        cell_type_plot = cell_type_plot,
        lineage_plot = lineage_plot,
        cell_type_label_plot = cell_type_label_plot,
        lineage_label_plot = lineage_label_plot,
        dendro_plot = dendro_plot
    ))
}

plot_differential_features_heatmap <- function(dea_results_path, fig_path, fdr_threshold, log2FC_threshold, title = NULL, 
                                               feature, ct_clst_method, ct_clst_dist, feature_clst_method,
                                               feature_clst_dist, q_mask=0, label_box_size_factor = 1, n_clusters = 25,
                                               max_marker_labels = 30, test_marker_annot = FALSE) {
    if (feature == 'Genes') {
        feature_col <- 'feature_name'
        y_label <- 'Differentially expressed genes'
    } else if (feature == 'Regions') {
        feature_col <- 'feature'
        y_label <- 'Differentially accessible regions\n(mapped to closest gene)'
    }
    
    # Get data for both up and down regulated features
    heatmap_df <- get_top_differential_features(dea_results_path, fdr_threshold, log2FC_threshold)

    if (test_marker_annot) {
        stopifnot(feature == 'Regions')
        # add some made up regions that have very high or low logFC for all cell types
        cell_types <- unique(heatmap_df$group)
        new_regions <- data.frame(
            feature = rep(paste0("test_region", 1:2), each = length(cell_types)),  # regions different names
            feature_name = rep("test_region", length(cell_types)),  # all map to the same gene to test duplicates
            logFC = c(rep(20, length(cell_types)), rep(-20, length(cell_types))),
            group = cell_types
        )
        # Add NA values for any missing columns when binding the dataframes
        new_regions[setdiff(names(heatmap_df), names(new_regions))] <- NA
        heatmap_df <- rbind(heatmap_df, new_regions)

        # add to HAEMATOPOIESIS_MARKERS
        HAEMATOPOIESIS_MARKERS <- c(HAEMATOPOIESIS_MARKERS, "test_region")

    }

    # For Regions, avoid adding display suffixes; instead, use a simple unique key for pivoting
    # based on the (feature, feature_name) pair so that duplicate regions remain duplicated rows.
    if (feature == 'Regions') {
        heatmap_df <- heatmap_df %>%
            mutate(feature_id = paste0(feature, "||", feature_name))
        feature_col_to_plot <- 'feature_id'
    } else {
        feature_col_to_plot <- feature_col
    }

    # use hclust and dendsort to order rows and columns
    mat <- heatmap_df %>%
        select(group, all_of(feature_col_to_plot), logFC) %>%
        pivot_wider(names_from = group, values_from = logFC, values_fill = 0) %>%
        column_to_rownames(feature_col_to_plot) %>%
        as.matrix()

    ord <- compute_diagonal_cluster_order(
        mat,
        row_clst_dist = feature_clst_dist,
        row_clst_method = feature_clst_method,
        col_clst_dist = ct_clst_dist,
        col_clst_method = ct_clst_method,
        n_clusters = n_clusters
    )

    feature_vector <- heatmap_df[[feature_col_to_plot]]
    heatmap_df[[feature_col_to_plot]] <- factor(feature_vector, levels = rev(ord$row_order))
    heatmap_df$group <- factor(heatmap_df$group, levels = ord$col_order)
    
    # Annotation bars and dendrogram for cell type and lineage
    ann_plots <- build_column_annotations(ord$col_order, ord$col_dendro, title = title)

    # Quantile masking to remove outliers from color scale of heatmap
    if(q_mask > 0) {
        upper_limit <- quantile(heatmap_df$logFC, probs = 1 - q_mask, na.rm=TRUE)
        lower_limit <- quantile(heatmap_df$logFC, probs = q_mask, na.rm=TRUE)
        heatmap_df$logFC <- ifelse(heatmap_df$logFC < lower_limit, lower_limit, heatmap_df$logFC)
        heatmap_df$logFC <- ifelse(heatmap_df$logFC > upper_limit, upper_limit, heatmap_df$logFC)
    }

    plot_limits <- c(-1, 1) * max(abs(heatmap_df$logFC), na.rm=TRUE)

    dendro_plot <- ann_plots$dendro_plot

    heatmap_plot <- ggplot(heatmap_df, aes_string(x = "group", y = feature_col_to_plot, fill = "logFC")) +
        geom_tile(linewidth = 0) +
        scale_fill_distiller(
            palette = "RdBu",
            limits = plot_limits,
            guide = guide_colorbar(direction = "vertical"),  # , barheight = grid::unit(2, "cm"), barwidth = grid::unit(0.25, "cm")
            name = "log2(FC)"
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        labs(x = 'Cell type') +
        MrBiomics_theme() +
        theme(
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
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
        all_features <- levels(heatmap_df[[feature_col_to_plot]])
        selected_features <- as.character(marker_labels_df[[feature_col_to_plot]])

        # If there are more candidate markers than allowed, choose top significant markers per cell type
        if (length(selected_features) > max_marker_labels) {
            # backup original candidates in case filtering removes all
            marker_labels_df_backup <- marker_labels_df
            label_text_map_backup <- label_text_map

            num_cell_types <- length(levels(heatmap_df$group))
            n_markers_per_ct <- floor(max_marker_labels / max(1, num_cell_types))
            if (n_markers_per_ct < 1) n_markers_per_ct <- 1

            selected_features_new <- c()
            marker_labels_no_dupl <- marker_labels_df
            for (ct in levels(heatmap_df$group)) {
                ct_df <- marker_labels_no_dupl %>%
                    filter(group == ct)
                if (nrow(ct_df) > n_markers_per_ct) {
                    ct_df <- ct_df %>%
                        slice_min(order_by = abs(logFC), n = n_markers_per_ct)
                }
                selected_features_new <- c(selected_features_new, as.character(ct_df[[feature_col_to_plot]]))
                # drop already used markers to avoid duplicates between cell types
                marker_labels_no_dupl <- marker_labels_no_dupl %>%
                    filter(!.data[[feature_col_to_plot]] %in% ct_df[[feature_col_to_plot]])
            }
            selected_features <- selected_features_new
        }


        # adapt based on number of features in the heatmap and how many markers there can be at max
        label_box_size <- nrow(mat) / max_marker_labels * label_box_size_factor
        
        pos_df <- compute_nonoverlapping_label_positions(
            heatmap_rows = all_features,
            selected_features = selected_features,
            box_size = label_box_size
        )

        n_rows <- length(all_features)
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
            geom_text(aes(x = 0.8, y = label_y_plot, label = label), hjust = 1) +
            # Fixed limits matching number of rows; flush to edges
            scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = c(0, 0)) +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            labs(y = y_label) +
            MrBiomics_theme() +
            theme(
                # remove all grid and ticks and box and everything except the y-axis label
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_text(margin = margin(l = 10)),
                panel.border = element_blank(),
                # make transparent so that other plots are not cut off
                panel.background = element_rect(fill = "transparent", color = NA),
                plot.background = element_rect(fill = "transparent", color = NA)
            ) +
            coord_cartesian(clip = "off") +
            # margin only on the left so that the labels are not cut off
            theme(plot.margin = margin(0, 0, 0, 10))
    } else {
        marker_label_plot <- plot_spacer()
    }
    
    # Assemble the final plot using patchwork's wrap_plots for robust layout
    plot_list <- list(
        plot_spacer(),          dendro_plot,
        ann_plots$cell_type_label_plot,   ann_plots$cell_type_plot,
        ann_plots$lineage_label_plot,     ann_plots$lineage_plot,
        marker_label_plot,      heatmap_plot
    )
    gp <- wrap_plots(plot_list, ncol = 2,
                     widths = c(1.2, 3),
                     heights = c(0.6, 0.35, 0.35, 10)) +
        plot_layout(guides = "collect") +
        plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) &
        theme(
            legend.position = "right",
            legend.box = "vertical",
            plot.margin = margin(0, 0, 0, 0),
            panel.spacing = grid::unit(0, "pt")
        )

    ggsave_all_formats(path = fig_path,
                       plot = gp,
                       width = PLOT_SIZE_3_PER_ROW,
                       height = PLOT_SIZE_3_PER_ROW, 
                       dpi = 1000)
    
    return(heatmap_plot)
}

########################################################################################################################
#### ENRICHMENT HEATMAP #################################################################################################
########################################################################################################################
filter_top_terms <- function(df, fdr_threshold, tissues_to_keep = NULL, top_n_per_name = 2) {
    df_sig <- df %>%
      filter(score > 0, statistic < fdr_threshold) 

    if (!is.null(tissues_to_keep)) {
        df_sig <- df_sig %>%
            filter(grepl(paste(tissues_to_keep, collapse = "|"), Term))
    }

    # Identify top term(s) for each cell type
    top_hits <- df_sig %>%
      group_by(name) %>%
      slice_max(order_by = score, n = top_n_per_name, with_ties = FALSE) %>%
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

prepare_for_heatmap <- function(df_formatted, fdr_threshold, tissues_to_keep = NULL, top_n_per_name = 2) {
    df_top <- filter_top_terms(df = df_formatted, fdr_threshold = fdr_threshold, tissues_to_keep = tissues_to_keep, top_n_per_name = top_n_per_name)

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

    # for the GOBP terms, shorten the names by removing the text in the brackets at the end (incl the brackets)
    short_term_levels <- gsub("\\s*\\(.*\\)", "", short_term_levels)

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

plot_enrichment_heatmap <- function(heatmap_df, fig_path, fill_lab, size_lab, title = NULL, ylabel = NULL, width = NULL) {
    # mask those values where the log statistic is tiny with NaN, to drop them in the heatmap
    heatmap_df$neg_log10_statistic <- ifelse(heatmap_df$neg_log10_statistic < 1, NaN, heatmap_df$neg_log10_statistic)
    lim <- c(-1, 1) * max(abs(heatmap_df$score), na.rm = TRUE)

    enrichment_plot <- ggplot(heatmap_df, aes(x = name, y = Term, size=neg_log10_statistic, fill = score)) +
      geom_point(shape=21, stroke=0.25) +
      # add star for significance (overlay centered star points)
      geom_point(
        data = subset(heatmap_df, sig),
        aes(x = name, y = Term),
        inherit.aes = FALSE,
        shape = 8,
        size = 0.5,
        alpha = 0.8,
        stroke = 0.25,
        color = "white",
        show.legend = FALSE
      ) +
      scale_fill_distiller(
        palette = "RdBu",
        limits = lim,
        name = fill_lab,
        guide = guide_colorbar(direction = "vertical")
      ) +
      scale_size_continuous(name = size_lab, range = c(0.5, 3)) +
    #   guides(size = guide_legend(keyheight = grid::unit(0.4, "lines"), keywidth = grid::unit(0.6, "lines"))) +
      labs(title = title, x = "Cell type", y = ylabel) +
      # ensure square tiles
      coord_fixed() +
      MrBiomics_theme() + 
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = FONT_SIZE_SMALL),
      )

    if (is.null(width)) {
        width <- PLOT_SIZE_3_PER_ROW + 2
    }

    # Save plot
    ggsave_all_formats(path = fig_path,
                       plot = enrichment_plot,
                       width = width,
                       height = PLOT_SIZE_3_PER_ROW)
    
    return(enrichment_plot)
}

plot_clustered_enrichment_heatmap <- function(heatmap_df,
                                              fig_path,
                                              fill_lab,
                                              size_lab,
                                              title = NULL,
                                              ylabel = NULL,
                                              ct_clst_dist = "euclidean",
                                              ct_clst_method = "ward.D2",
                                              term_clst_dist = "euclidean",
                                              term_clst_method = "ward.D2",
                                              n_clusters = 25,
                                              width = NULL) {

    # Build term-by-celltype matrix of enrichment scores for clustering
    mat_df <- heatmap_df %>%
        select(name, Term, score) %>%
        tidyr::pivot_wider(names_from = name, values_from = score, values_fill = 0)
    mat <- mat_df %>% column_to_rownames("Term") %>% as.matrix()

    ord <- compute_diagonal_cluster_order(
        mat,
        row_clst_dist = term_clst_dist,
        row_clst_method = term_clst_method,
        col_clst_dist = ct_clst_dist,
        col_clst_method = ct_clst_method,
        n_clusters = n_clusters
    )

    # Apply orders
    heatmap_df$Term <- factor(as.character(heatmap_df$Term), levels = rev(ord$row_order))
    heatmap_df$name <- factor(as.character(heatmap_df$name), levels = ord$col_order)

    # Mask tiny p-values from size scale
    heatmap_df$neg_log10_statistic <- ifelse(heatmap_df$neg_log10_statistic < 1, NaN, heatmap_df$neg_log10_statistic)

    # Color scale limits based on scores
    lim <- c(-1, 1) * max(abs(heatmap_df$score), na.rm = TRUE)

    # FIXME break long term names into multiple lines, breaking at the first space after 20 characters

    # Bubble heatmap with significance stars (y-axis labels moved to separate plot)
    bubble <- ggplot(heatmap_df, aes(x = name, y = Term, size = neg_log10_statistic, fill = score)) +
        geom_point(shape = 21, stroke = 0.25) +
        # add asterisks for significance (overlay centered star points)
        geom_point(
            data = subset(heatmap_df, sig),
            aes(x = name, y = Term),
            inherit.aes = FALSE,
            shape = 8,
            size = 0.5,
            alpha = 0.8,
            stroke = 0.25,
            color = "white",
            show.legend = FALSE
        ) +
        scale_fill_distiller(
            palette = "RdBu",
            limits = lim,
            name = fill_lab,
            guide = guide_colorbar(direction = "vertical")
        ) +
        scale_size_continuous(name = size_lab, range = c(0.5, 2.3)) +
        scale_x_discrete(expand = expansion(add = 0.6)) +
        scale_y_discrete(expand = expansion(add = 0.6)) +
        labs(title = NULL, x = "Cell type", y = NULL) +
        coord_cartesian(clip = "off") +
        MrBiomics_theme() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(t = -18, r = 4, b = -8, l = 4, unit = "pt")
        )

    # Top annotations + dendrogram
    ann_plots <- build_column_annotations(ord$col_order, ord$col_dendro, title = title)

    # Left-side term labels (y-axis moved out of bubble plot)
    term_levels <- levels(heatmap_df$Term)
    if (is.null(term_levels)) {
        term_levels <- rev(ord$row_order)
    }
    n_rows <- length(term_levels)
    term_labels_df <- data.frame(Term = term_levels, y = seq_len(n_rows))
    term_label_plot <- ggplot(term_labels_df, aes(y = y)) +
        geom_text(aes(x = 0.98, label = Term), hjust = 1) +
        scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = c(0, 0)) +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        labs(y = ylabel) +
        MrBiomics_theme() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_text(margin = margin(l = 20)),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA),
            plot.margin = margin(0, 0, 0, 20)
        ) +
        coord_cartesian(clip = "off")

    # Assemble    
    plot_list <- list(
        plot_spacer(),              ann_plots$dendro_plot,
        ann_plots$cell_type_label_plot, ann_plots$cell_type_plot,
        ann_plots$lineage_label_plot,   ann_plots$lineage_plot,
        term_label_plot,            bubble
    )

    gp <- wrap_plots(plot_list, ncol = 2,
                     widths = c(1.5, 1),
                     heights = c(0.6, 0.35, 0.35, 10)) +
        plot_layout(guides = "collect") +
        plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) &
        theme(legend.position = "right", legend.box = "vertical", plot.margin = margin(0, 0, 0, 0),
               panel.spacing = grid::unit(0, "pt"), legend.margin = margin(l = -10))

    if (is.null(width)) {
        width <- PLOT_SIZE_3_PER_ROW + 2
    }

    ggsave_all_formats(path = fig_path,
                       plot = gp,
                       width = width,
                       height = PLOT_SIZE_3_PER_ROW)

    return(gp)
}

########################################################################################################################
#### Lineage reconstruction via crossprediction #######################################################################
########################################################################################################################
# Generate hierarchical layout coordinates based on a predefined tree
generate_hierarchy_layout_coordinates <- function(nodes, children_list, root_node,
                                                 spacing_between_layers , jitter_step) {
    nodes_in_current_layer <- c(root_node)
    y <- 0
    x <- 0
    coordinates <- list()
    coordinates[[root_node]] <- c(x, y)

    while (length(nodes_in_current_layer) > 0) {
        y <- y - spacing_between_layers
        children_in_next_layer <- c()

        for (current_node in nodes_in_current_layer) {
            if (current_node %in% names(children_list)) {
                children <- children_list[[current_node]]
                children_in_next_layer <- c(children_in_next_layer, children)
            }
        }

        children_in_next_layer <- unique(children_in_next_layer)

        if (length(children_in_next_layer) > 0) {
            x_positions <- 1:length(children_in_next_layer) - 1
            x_positions <- x_positions - ((length(children_in_next_layer) - 1) / 2)
            x_positions <- x_positions * (spacing_between_layers / (length(children_in_next_layer)))
            for (i in seq_along(children_in_next_layer)) {
                y <- y + jitter_step
                jitter_step <- jitter_step * -1
                coordinates[[children_in_next_layer[i]]] <- c(x_positions[i], y)
            }
        }

        nodes_in_current_layer <- children_in_next_layer
    }

    layoutCoordinates <- do.call(rbind, coordinates)
    layoutCoordinates <- layoutCoordinates[nodes, , drop = FALSE]
    rownames(layoutCoordinates) <- NULL
    colnames(layoutCoordinates) <- c("x", "y")
    return(layoutCoordinates)
}

# Build and save the crossprediction plot from an adjacency matrix CSV
plot_crossprediction_from_adjacency <- function(adjacency_matrix_path,
                                               fig_path,
                                               coordinates_out_path,
                                               modality_label,
                                               lineage_tree_cut_off = 0.05,
                                               hierarchy_coordinates = TRUE,
                                               root_node = "HSC",
                                               spacing_between_layers = 5,
                                               jitter_step = 1.5,
                                               outcome_title = "Compared to\nliterature",
                                               node_shape = 19,
                                               point_size = 8, 
                                               tp_name = "Consistent",
                                               fp_name = "Additional",
                                               fn_name = "Missing") {



    # Load adjacency matrix and harmonize labels
    adj_mtx <- data.frame(fread(file.path(adjacency_matrix_path), header = TRUE), row.names = 1)
    rownames(adj_mtx) <- recode(rownames(adj_mtx), !!!DATA_TO_CELL_TYPE_COLORS_MAPPING)
    colnames(adj_mtx) <- recode(colnames(adj_mtx), !!!DATA_TO_CELL_TYPE_COLORS_MAPPING)

    adjacencyMatrix <- as.matrix(adj_mtx)
    adjacencyMatrix[adjacencyMatrix < lineage_tree_cut_off] <- 0

    nodes <- rownames(adj_mtx)
    node_colors <- CELL_TYPE_COLORS[nodes]

    # Create layout coordinates
    if (hierarchy_coordinates) {
        layoutCoordinates <- generate_hierarchy_layout_coordinates(
            nodes = nodes,
            children_list = children_list,
            root_node = root_node,
            spacing_between_layers = spacing_between_layers,
            jitter_step = jitter_step
        )
    } else {
        # Use Fruchterman-Reingold layout
        if (.Platform$OS.type == "windows") {
            pdf(file = "NUL")
        } else {
            pdf(file = "/dev/null")
        }
        layoutCoordinates <- gplot(adjacencyMatrix, mode = "fruchtermanreingold")
        dev.off()
    }

    # Build symmetric ground-truth adjacency from lineage tree
    ground_truth_adj_mat <- matrix(0, nrow = length(nodes), ncol = length(nodes))
    rownames(ground_truth_adj_mat) <- nodes
    colnames(ground_truth_adj_mat) <- nodes
    for (parent in names(children_list)) {
        children <- children_list[[parent]]
        children <- intersect(children, nodes)
        if (length(children) > 0 && parent %in% nodes) {
            ground_truth_adj_mat[parent, children] <- 1
            ground_truth_adj_mat[children, parent] <- 1
        }
    }

    # Capture full set of undirected ground-truth edges before subtracting predicted edges
    ground_truth_pairs_all <- ground_truth_adj_mat %>%
        reshape2::melt() %>%
        dplyr::filter(value == 1) %>%
        dplyr::mutate(
            a = pmin(as.character(Var1), as.character(Var2)),
            b = pmax(as.character(Var1), as.character(Var2))
        ) %>%
        dplyr::select(a, b) %>%
        dplyr::distinct()
    gt_keys <- paste(ground_truth_pairs_all$a, ground_truth_pairs_all$b, sep = ">")

    # Convert to list of ties only (>= cutoff) and remove ground-truth edges that are predicted
    adjacencyList <- reshape2::melt(adjacencyMatrix)
    adjacencyList <- adjacencyList[adjacencyList$value >= lineage_tree_cut_off, , drop = FALSE]
    if (nrow(adjacencyList) > 0) {
        for (i in seq_len(nrow(adjacencyList))) {
            v1 <- adjacencyList[i, 'Var1']
            v2 <- adjacencyList[i, 'Var2']
            ground_truth_adj_mat[v1, v2] <- 0
            ground_truth_adj_mat[v2, v1] <- 0
        }
    }

    # Unstack remaining ground-truth edges and add coordinates
    layoutCoordinates_df <- data.frame(layoutCoordinates)
    layoutCoordinates_df['node_name'] <- nodes
    ground_truth_adj_mat_long <- ground_truth_adj_mat %>%
        reshape2::melt() %>%
        dplyr::filter(value == 1) %>%
        dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
        dplyr::filter(Var1 < Var2) %>%
        dplyr::select(-value) %>%
        dplyr::left_join(layoutCoordinates_df, by = c("Var1" = "node_name"), copy = TRUE) %>%
        dplyr::rename(x_start = x, y_start = y) %>%
        dplyr::left_join(layoutCoordinates_df, by = c("Var2" = "node_name"), copy = TRUE) %>%
        dplyr::rename(x_end = x, y_end = y)

    # Edge path generator with type and probability
    edgeMaker <- function(whichRow, len = 100) {
        fromC <- layoutCoordinates[adjacencyList[whichRow, 1], ]
        toC <- layoutCoordinates[adjacencyList[whichRow, 2], ]

        edge <- data.frame(bezier(c(fromC[1], toC[1]), c(fromC[2], toC[2]), evaluation = len))
        edge$Sequence <- 1:len
        edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")

        # Determine undirected pair key and type based on ground truth membership
        a_val <- pmin(as.character(adjacencyList[whichRow, 'Var1']), as.character(adjacencyList[whichRow, 'Var2']))
        b_val <- pmax(as.character(adjacencyList[whichRow, 'Var1']), as.character(adjacencyList[whichRow, 'Var2']))
        pair_key <- paste(a_val, b_val, sep = ">")
        edge$edge_type <- ifelse(pair_key %in% gt_keys, tp_name, fp_name)

        # Interpolate between the weights of the current edge and its counterpart
        x <- c(1, len)
        y <- c(adjacencyMatrix[adjacencyList[whichRow, 'Var1'], adjacencyList[whichRow, 'Var2']],
               adjacencyMatrix[adjacencyList[whichRow, 'Var2'], adjacencyList[whichRow, 'Var1']])
        edge$Probability <- approx(x, y, xout = 1:len)$y

        return(edge)
    }

    # Generate edge paths
    allEdges <- lapply(seq_len(nrow(adjacencyList)), edgeMaker, len = 500)
    if (length(allEdges) > 0) {
        allEdges <- do.call(rbind, allEdges)
    } else {
        allEdges <- data.frame()
    }

    # Build plot
    p <- ggplot()
    if (nrow(ground_truth_adj_mat_long) > 0) {
        p <- p + geom_segment(
            data = ground_truth_adj_mat_long,
            aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                colour = fn_name, linetype = fn_name),
            size = 1
        )
    }
    if (nrow(allEdges) > 0) {
        p <- p + geom_path(
            data = allEdges,
            aes(x = x, y = y, group = Group, size = Probability, colour = edge_type, linetype = edge_type),
            na.rm = TRUE
        )
    }
    p <- p + geom_point(data = data.frame(layoutCoordinates),
                        aes(x = x, y = y), shape = node_shape, size = point_size, color = node_colors)
    p <- p + geom_text(data = data.frame(layoutCoordinates),
                       aes(x = x, y = y, label = nodes), color = 'white', hjust = 0.5, vjust = 0.5, fontface = "bold")

    p <- p +
        scale_colour_manual(
            name = outcome_title,
            values = setNames(c("black", "grey", "grey"), c(tp_name, fp_name, fn_name)),
            breaks = c(tp_name, fp_name, fn_name),
            labels = c(tp_name, fp_name, fn_name)
        ) +
        scale_linetype_manual(
            name = outcome_title,
            values = setNames(c("solid", "solid", "11"), c(tp_name, fp_name, fn_name)),
            breaks = c(tp_name, fp_name, fn_name),
            labels = c(tp_name, fp_name, fn_name)
        )

    p <- p + scale_size(range = c(0.1, 3))

    # Add padding so nodes/edges aren't clipped at plot boundaries
    p <- p +
        scale_x_continuous(expand = expansion(mult = 0.12)) +
        scale_y_continuous(expand = expansion(mult = 0.12))

    # Add modality label as title text
    p <- p + geom_text(aes(x = min(data.frame(layoutCoordinates)$x), y = 1, label = modality_label),
                       hjust = 0, vjust = 1, fontface = "bold", size = FONT_SIZE_NORMAL / .pt)

    crosspred_plot <- p +
        guides(
            size = guide_legend(title = 'Average\ncross-prediction\nprobability'),
            colour = guide_legend(),
            linetype = guide_legend()
        ) +
        MrBiomics_void()

    # Save plot and coordinates
    ggsave_all_formats(path = fig_path,
                       plot = crosspred_plot,
                       width = PLOT_SIZE_3_PER_ROW+2,
                       height = PLOT_SIZE_3_PER_ROW)

    write.csv(layoutCoordinates, file.path(coordinates_out_path))

    return(crosspred_plot)
}