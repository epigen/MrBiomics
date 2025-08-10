
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
# output
rna_umap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/umap.pdf"
rna_dea_heatmap_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/dea_heatmap.pdf"
rna_enrichment_path <- "/nobackup/lab_bock/projects/MrBiomics/paper/CorcesRNA/enrichment.pdf"
# params
# enrichment analysis
fdr_threshold <- 0.05
TOP_N_GENES <- 10

######### UMAPs (from unsupervised analysis) ############
# Create UMAP plots
rna_umap_plot <- umap_plot(CorcesRNA_umap_coords_path, rna_umap_path, title = "RNA")



######### DEA HEATMAP ############
rna_dea_heatmap_plot <- plot_dea_heatmap(
    dea_results_path = CorcesRNA_dea_OvA_path,
    fig_path = rna_dea_heatmap_path,
    top_n_genes = TOP_N_GENES,
    fdr_threshold = fdr_threshold,
    title = "RNA one-versus-all DEA"
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
rna_df_formatted <- create_rna_enrichment_df(CorcesRNA_enrichment_results_path, fdr_threshold)
rna_heatmap_df <- prepare_for_heatmap(df_formatted = rna_df_formatted, fdr_threshold = fdr_threshold)
rna_enrichment_plot <- plot_enrichment_heatmap(
    heatmap_df = rna_heatmap_df,
    fig_path = rna_enrichment_path,
    fill_lab = "NES",
    size_lab = "-log10(FDR q-value)",
    title = "RNA"
)


# ######### Lineage reconstruction using crossprediction ############
# # plot adapted from: https://gist.github.com/dsparks/4331058

# # load adjacency matrix
# adj_mtx <- data.frame(fread(file.path(crossprediction_adj_mtx_path), header=TRUE), row.names=1)
# dim(adj_mtx)
# # head(adj_mtx)

# # load previous best layout coordinates (TODO)
# # coordinates <- as.matrix(data.frame(fread(file.path(crossprediction_adj_mtx_path), header=TRUE), row.names=1))
# # dim(coordinates)
# # head(coordinates)

# # make cross prediction graph plot

# # set colors & shapes
# nodes <- rownames(adj_mtx)
# # node_cond <- rep(cond, times = dim(adj_mtx)[1])
# node_colors <- "#F8766D" #KO_col[nodes]
# node_shape <- 19 # cond_shapes[node_cond]
# border_color <- "#000000" #cond_colors[node_cond]
# stroke_size <- 1.5
# point_size <- 5

# # data
# adjacencyMatrix <- as.matrix(adj_mtx)
# # parameters
# curved=FALSE

# adjacencyMatrix[adjacencyMatrix<cut_off] <- 0

# # plot graph

# # Empty ggplot2 theme
# new_theme_empty <- theme_bw() 
# new_theme_empty$line <- element_blank()
# new_theme_empty$rect <- element_blank()
# new_theme_empty$strip.text <- element_blank()
# new_theme_empty$axis.text <- element_blank()
# new_theme_empty$plot.title <- element_blank()
# new_theme_empty$axis.title <- element_blank()
# new_theme_empty$plot.margin <- structure(c(0, 0, 0, 0), unit = "lines",
#                                          valid.unit = 3L, class = "unit")

# if (exists("coordinates")){
#     layoutCoordinates <- coordinates
# } else{
#     pdf(file = if (.Platform$OS.type == "windows") "NUL" else "/dev/null")  # route to null device (ie trash)
#     layoutCoordinates <- gplot(adjacencyMatrix, mode="fruchtermanreingold")    # Get graph layout coordinates
#     dev.off()
# }


# adjacencyList <- reshape2::melt(adjacencyMatrix)  # Convert to list of ties only
# adjacencyList <- adjacencyList[adjacencyList$value >= cut_off, ] # prune weak edges


# # Function to generate paths between each connected node
# edgeMaker <- function(whichRow, len = 100, curved = TRUE){
#     fromC <- layoutCoordinates[adjacencyList[whichRow, 1], ]  # Origin
#     toC <- layoutCoordinates[adjacencyList[whichRow, 2], ]  # Terminus

#     # Add curve:
#     graphCenter <- colMeans(layoutCoordinates)  # Center of the overall graph
#     bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
#     distance1 <- sum((graphCenter - bezierMid)^2)
#     if(distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)){
#     bezierMid <- c(toC[1], fromC[2])
#     }  # To select the best Bezier midpoint
#     bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
#     if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve

#     edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
#                             c(fromC[2], bezierMid[2], toC[2]),  # X & y
#                             evaluation = len))  # Bezier path coordinates
#     edge$Sequence <- 1:len  # For size and colour weighting in plot
#     edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")

#     # add a vector 'Probability' that linearly interpolates between the weights of the current edge and its counterpart
#     x <- c(1,len)
#     y <- c(adjacencyMatrix[adjacencyList[whichRow,'Var1'],adjacencyList[whichRow,'Var2']],adjacencyMatrix[adjacencyList[whichRow,'Var2'],adjacencyList[whichRow,'Var1']])
#     edge$Probability <- approx(x,y,xout=1:len)$y

#   return(edge)
# }

# # Generate a (curved) edge path for each pair of connected nodes
# allEdges <- lapply(1:nrow(adjacencyList), edgeMaker, len = 500, curved = curved)
# allEdges <- do.call(rbind, allEdges)  # a fine-grained path ^, with bend ^

# zp1 <- ggplot(allEdges) 
# zp1 <- zp1 + geom_path(aes(x = x, y = y, group = Group,  # Edges with gradient
#                            colour = Sequence, size = Probability))  # and taper

# zp1 <- zp1 + geom_point(data = data.frame(layoutCoordinates),  # Add nodes
# #                         aes(x = x, y = y),  pch = 21, size = point_size, fill = node_colors, colour = border_color, stroke = stroke_size) # fill & borders are informative
#                         aes(x = x, y = y), shape=node_shape, size = point_size, color = node_colors)#+scale_shape_manual(name = "condition", labels = names(cond_shapes),values = cond_shapes) # fill and shape are informative

# # zp1 <- zp1 + geom_text(data = data.frame(layoutCoordinates), aes(x = x, y = y, label=nodes),#rownames(adjacencyMatrix)),
# #                        hjust=0.5, vjust=-1)

# zp1 <- zp1 + geom_label(data = data.frame(layoutCoordinates), aes(x = x, y = y, label=nodes), fill=node_colors,#rownames(adjacencyMatrix)),
#                         hjust=0.5, vjust=-0.5)

# zp1 <- zp1 + scale_colour_gradient(low = gray(0), high = gray(0), guide = "none") # Customize gradient

# zp1 <- zp1 + scale_size(range = c(1/10, point_size-1))#, guide = "none")  # Customize taper

# crosspred_p <- zp1 + 
#         new_theme_empty + 
#         #theme(legend.position = c(0.5, 0.5))+ 
#         guides(size = guide_legend(title='Average\ncross-prediction\nprobability')) # Clean up plot


# # plot the graph
# width <- 6
# height <- 6
# options(repr.plot.width=width, repr.plot.height=height)

# # crosspred_p


# # # save layout coordinates
# # write.csv(layoutCoordinates, file.path("results/KO150/KO_classifier",cond,'layoutCoordinates.csv'))

# # save plot
# ggsave_all_formats(path=crossprediction_plot_path,
#                    plot=crosspred_p,
#                    width=width,
#                    height=height
# )
