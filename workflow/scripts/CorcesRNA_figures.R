
# set correct working directory -> project/repo folder (only for dev, not prod)
# getwd()
# setwd('..')
# getwd()

#### libraries
# source this for libraries, MrBiomics theme and utility functions
source("workflow/scripts/figure_theme.R")
source("workflow/scripts/figure_utils.R")

# FIXME snakemakeify
# ## libraries are loaded in the source files
# source(snakemake@params[["figure_theme_path"]])
# figure_utils.R
# ## input
# enrichment_results_path <- snakemake@input[["enrichment_results"]]
# crossprediction_adj_mtx_path <- snakemake@input[["crossprediction_adj_mtx"]]
# ## output
# enrichment_plot_path <- snakemake@output[["enrichment_plot"]]
# crossprediction_plot_path <- snakemake@output[["crossprediction_plot"]]
# ## params
# # enrichment analysis
# fdr_threshold <- snakemake@params[["fdr_threshold"]]
# # lineage reconstructions
# cut_off <- snakemake@params[["cut_off"]]

# input
CorcesRNA_umap_coords_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/unsupervised_analysis/normCQN_integrated/UMAP/UMAP_correlation_15_0.1_2_data.csv"
CorcesRNA_dea_OvA_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/dea_limma/normCQN_OvA_cell_type/results.csv"
CorcesRNA_enrichment_results_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_all.csv"
CorcesRNA_crossprediction_adj_mtx_path <- "/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"
# output
rna_umap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/umap.pdf"
rna_dea_heatmap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/differential_heatmap.pdf"
rna_enrichment_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/enrichment.pdf"
rna_crossprediction_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/crossprediction.pdf"
rna_crossprediction_coordinates_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/crossprediction_coordinates.csv"
# params
# enrichment analysis
fdr_threshold <- 0.05
log2FC_threshold <- 3
lineage_tree_cut_off <- 0.05
hierarchy_coordinates <- TRUE

# ######### UMAPs (from unsupervised analysis) ############
# # Create UMAP plots
# rna_umap_plot <- umap_plot(CorcesRNA_umap_coords_path, rna_umap_path, title = "RNA")


# ######### DEA HEATMAP ############
# rna_dea_heatmap_plot <- plot_differential_features_heatmap(
#     dea_results_path = CorcesRNA_dea_OvA_path,
#     fig_path = rna_dea_heatmap_path,
#     fdr_threshold = fdr_threshold,
#     log2FC_threshold = log2FC_threshold,
#     title = "RNA",
#     feature = 'Genes',
#     ct_clst_dist = "euclidean",
#     ct_clst_method = "ward.D2",
#     feature_clst_dist = "maximum",  # maximum to focus on the most differentially expressed groups when sorting genes
#     feature_clst_method = "ward.D2",
#     q_mask = 0.025,
#     label_box_size_factor = 1
# )

# ######### ENRICHMENT HEATMAP ############
# # Function to create enrichment heatmap
# create_rna_enrichment_df <- function(enrichment_results_path, fdr_threshold = 0.05) {
#     # Load enrichment analysis result
#     df <- data.frame(fread(file.path(enrichment_results_path), header=TRUE))

#     df_formatted <- df %>%
#         rename(statistic = FDR_q_val, score = NES) %>%
#         mutate(name = recode(name, !!!DATA_TO_CELL_TYPE_COLORS_MAPPING))
    
#     return(df_formatted)
# }

# # Create enrichment heatmaps
# rna_df_formatted <- create_rna_enrichment_df(CorcesRNA_enrichment_results_path, fdr_threshold)
# rna_heatmap_df <- prepare_for_heatmap(df_formatted = rna_df_formatted, fdr_threshold = fdr_threshold)
# rna_enrichment_plot <- plot_enrichment_heatmap(
#     heatmap_df = rna_heatmap_df,
#     fig_path = rna_enrichment_path,
#     fill_lab = "NES",
#     size_lab = "-log10(FDR q-value)",
#     title = "RNA",
#     ylabel = "Enrichment term\n(preranked GSEA, Azimuth 2023)"
# )

######### Lineage reconstruction using crossprediction ############
# plot adapted from: https://gist.github.com/dsparks/4331058

# load adjacency matrix
adj_mtx <- data.frame(fread(file.path(CorcesRNA_crossprediction_adj_mtx_path), header=TRUE), row.names=1)
dim(adj_mtx)
# head(adj_mtx)

# rename colnames and rownames with nice names
rownames(adj_mtx) <- recode(rownames(adj_mtx), !!!DATA_TO_CELL_TYPE_COLORS_MAPPING)
colnames(adj_mtx) <- recode(colnames(adj_mtx), !!!DATA_TO_CELL_TYPE_COLORS_MAPPING)

children_list <- list(
    HSC = c("MPP"),
    MPP = c("LMPP", "CMP"),
    LMPP = c("CLP", "GMP"),
    CMP = c("GMP", "MEP"),
    GMP = c("Mono"),
    MEP = c("Ery"),
    CLP = c("CD4", "CD8", "NK", "B")
)

# set plot parameters
nodes <- rownames(adj_mtx)
node_colors <- CELL_TYPE_COLORS[nodes]
node_shape <- 19
stroke_max <- 5
point_size <- 12
fontsize <- 3
spacing_between_layers <- 5
layer_jitter <- 1
outcome_title <- "Compared to\nCorces et al. (2016)"
tp_name <- "Consistent"
fp_name <- "Additional"
fn_name <- "Missing"

adjacencyMatrix <- as.matrix(adj_mtx)

adjacencyMatrix[adjacencyMatrix < lineage_tree_cut_off] <- 0

new_theme_empty <- theme_bw() 
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, 0, 0), unit = "lines",
                                         valid.unit = 3L, class = "unit")

nodes_in_current_layer <- c('HSC')  # root node at top of the hierarchy
y <- 0
x <- 0
coordinates <- list('HSC' = c(x, y))
ground_truth_adj_mat <- matrix(0, nrow = length(nodes), ncol = length(nodes))
rownames(ground_truth_adj_mat) <- nodes
colnames(ground_truth_adj_mat) <- nodes
if (hierarchy_coordinates){
    while (length(nodes_in_current_layer) > 0){
        y <- y - spacing_between_layers
        children_in_next_layer <- c()
        for (current_node in nodes_in_current_layer){
            if (current_node %in% names(children_list)){
                children <- children_list[[current_node]]
                
                # symmetric adjacency matrix
                ground_truth_adj_mat[current_node, children] <- 1
                ground_truth_adj_mat[children, current_node] <- 1

                # remember children for next layer
                children_in_next_layer <- c(children_in_next_layer, children)
            }
        }
        
        # for duplicates, keep the first occurrence
        children_in_next_layer <- unique(children_in_next_layer)
        
        if (length(children_in_next_layer) > 0){
            # positions equally spaced
            x_positions <- 1:length(children_in_next_layer) - 1
            # center the positions
            x_positions <- x_positions - ((length(children_in_next_layer)-1) / 2)
            # scale the positions to be wider if there are more children in the layer
            x_positions <- x_positions * (spacing_between_layers/(length(children_in_next_layer)))
            # add jitter so that edges between nodes of the same layers don't overlap
            for (i in 1:length(children_in_next_layer)){
                y <- y + layer_jitter
                layer_jitter <- layer_jitter * -1
                coordinates[[children_in_next_layer[i]]] <- c(x_positions[i], y)
            }
        }
        # move to next layer
        nodes_in_current_layer <- children_in_next_layer
    }
    # make same data shape as one would get from gplot
    layoutCoordinates <- do.call(rbind, coordinates)
    layoutCoordinates <- layoutCoordinates[nodes, ]
    rownames(layoutCoordinates) <- NULL
    colnames(layoutCoordinates) <- c("x", "y")

} else {
    pdf(file = if (.Platform$OS.type == "windows") "NUL" else "/dev/null")  # route to null device (ie trash)
    layoutCoordinates <- gplot(adjacencyMatrix, mode="fruchtermanreingold")    # Get graph layout coordinates
    dev.off()
}

# capture full set of undirected ground-truth edges before subtracting predicted edges
ground_truth_pairs_all <- ground_truth_adj_mat %>% 
    reshape2::melt() %>%
    filter(value == 1) %>%
    mutate(
        a = pmin(as.character(Var1), as.character(Var2)),
        b = pmax(as.character(Var1), as.character(Var2))
    ) %>%
    select(a, b) %>%
    distinct()


adjacencyList <- reshape2::melt(adjacencyMatrix)  # Convert to list of ties only
adjacencyList <- adjacencyList[adjacencyList$value >= lineage_tree_cut_off, ] # prune weak edges

# keys for ground-truth undirected edge membership
gt_keys <- paste(ground_truth_pairs_all$a, ground_truth_pairs_all$b, sep=">")

# only keep ground truth edges that are not in the data driven adjacency matrix
for (i in 1:nrow(adjacencyList)){
    ground_truth_adj_mat[adjacencyList[i,'Var1'], adjacencyList[i,'Var2']] <- 0
    ground_truth_adj_mat[adjacencyList[i,'Var2'], adjacencyList[i,'Var1']] <- 0
}
# unstack the matrix to long format and only keep the edges that are in the adjacency matrix
# since edges are symmetric, drop the ones that are not in alphabetical order to avoid duplicates
ground_truth_adj_mat_long <- ground_truth_adj_mat %>% 
    reshape2::melt() %>%
    filter(value == 1) %>%
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
    filter(Var1 < Var2) %>%
    select(-value)

# add the coordinates for plotting
layoutCoordinates_incl_names <- data.frame(layoutCoordinates)
layoutCoordinates_incl_names['node_name'] <- nodes
ground_truth_adj_mat_long <- ground_truth_adj_mat_long %>%
    left_join(layoutCoordinates_incl_names, by = c("Var1" = "node_name"), copy = TRUE) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(layoutCoordinates_incl_names, by = c("Var2" = "node_name"), copy = TRUE) %>%
    rename(x_end = x, y_end = y)


# Function to generate paths between each connected node
edgeMaker <- function(whichRow, len = 100){
    fromC <- layoutCoordinates[adjacencyList[whichRow, 1], ]  # Origin
    toC <- layoutCoordinates[adjacencyList[whichRow, 2], ]  # Terminus

    edge <- data.frame(bezier(c(fromC[1], toC[1]),  # Generate
                            c(fromC[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
    edge$Sequence <- 1:len  # For size and colour weighting in plot
    edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")

    # compute undirected pair key and type based on ground truth membership
    a_val <- pmin(as.character(adjacencyList[whichRow, 'Var1']), as.character(adjacencyList[whichRow, 'Var2']))
    b_val <- pmax(as.character(adjacencyList[whichRow, 'Var1']), as.character(adjacencyList[whichRow, 'Var2']))
    pair_key <- paste(a_val, b_val, sep=">")
    edge$edge_type <- ifelse(pair_key %in% gt_keys, tp_name, fp_name)

    # add a vector 'Probability' that linearly interpolates between the weights of the current edge and its counterpart
    x <- c(1,len)
    y <- c(adjacencyMatrix[adjacencyList[whichRow,'Var1'],adjacencyList[whichRow,'Var2']],
           adjacencyMatrix[adjacencyList[whichRow,'Var2'],adjacencyList[whichRow,'Var1']])
    edge$Probability <- approx(x,y,xout=1:len)$y

  return(edge)
}

# Generate a edge path for each pair of connected nodes
allEdges <- lapply(1:nrow(adjacencyList), edgeMaker, len = 500)
allEdges <- do.call(rbind, allEdges)  # a fine-grained path ^, with bend ^
head(allEdges)

zp1 <- ggplot(allEdges) 
zp1 <- zp1 + geom_segment(
                          data = ground_truth_adj_mat_long,
                          aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                              colour = fn_name, linetype = fn_name),
                          size = 2)  # linetype: draw 2, skip 1

zp1 <- zp1 + geom_path(aes(x = x, y = y, group = Group, size = Probability, colour = edge_type, linetype = edge_type),
                       na.rm = TRUE)  # taper with explicit type

zp1 <- zp1 + geom_point(data = data.frame(layoutCoordinates),  # Add nodes
                        aes(x = x, y = y), shape=node_shape, size = point_size, color = node_colors)

zp1 <- zp1 + geom_text(data = data.frame(layoutCoordinates), aes(x = x, y = y, label=nodes), color='white', hjust=0.5,
                       vjust=0.5, size = fontsize, fontface = "bold")

zp1 <- zp1 +
    scale_colour_manual(
        name = outcome_title,
        values = setNames(
            c("black", "grey", "grey"),
            c(tp_name, fp_name, fn_name)
        ),
        breaks = c(tp_name, fp_name, fn_name),
        labels = c(tp_name, fp_name, fn_name)
    ) +
    scale_linetype_manual(
        name = outcome_title,
        values = setNames(
            c("solid", "solid", "11"),
            c(tp_name, fp_name, fn_name)
        ),
        breaks = c(tp_name, fp_name, fn_name),
        labels = c(tp_name, fp_name, fn_name)
    )

zp1 <- zp1 + scale_size(range = c(0.3, stroke_max))#, guide = "none")  # Customize taper

# add padding so nodes/edges aren't clipped at plot boundaries
zp1 <- zp1 +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    scale_y_continuous(expand = expansion(mult = 0.12))

# add modality as text to act as title
zp1 <- zp1 + geom_text(aes(x = min(data.frame(layoutCoordinates)$x), y = 1, label = "RNA"),
                       hjust = 0, vjust = 1, size = 8, fontface = "bold")

crosspred_p <- zp1 + 
        new_theme_empty + 
        #theme(legend.position = c(0.5, 0.5))+ 
        guides(size = guide_legend(title='Average\ncross-prediction\nprobability')) # Clean up plot



width <- 6
# save plot
ggsave_all_formats(path=rna_crossprediction_path,
                   plot=crosspred_p,
                   width=width,
                   height=PLOT_HEIGHT
)

# save layout coordinates
write.csv(layoutCoordinates, file.path(rna_crossprediction_coordinates_path))

