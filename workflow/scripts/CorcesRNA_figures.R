
# set correct working directory -> project/repo folder (only for dev, not prod)
# getwd()
# setwd('..')
# getwd()

#### libraries
# source this for libraries, MrBiomics theme and utility functions
# source("workflow/scripts/figure_theme.R")
source(snakemake@params[["figure_theme_path"]])

#### configs

# input
enrichment_results_path <- snakemake@input[["enrichment_results"]]
#"/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_all.csv" 

crossprediction_adj_mtx_path <- snakemake@input[["crossprediction_adj_mtx"]]
#"/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"

# output
enrichment_plot_path <- snakemake@output[["enrichment_plot"]]

crossprediction_plot_path <- snakemake@output[["crossprediction_graph"]]
#"/nobackup/lab_bock/projects/MrBiomics/results/CorcesRNA/special_analyses/crossprediction/normCQN_integrated_graph.png"

# params (TODO: consider putting in rule)
# enrichment analysis
fdr_threshold <- 0.05
# lineage reconstructions
cut_off <- 0.2 # when pruning at 0.2 both HVF and integrated produce the same graph


######### Unsupervised Analysis (PCA & UMAP) ############

######### ... ############

######### Enrichment analysis ############

# load enrichment analysis result
enrichment_results <- data.frame(fread(file.path(enrichment_results_path), header=TRUE))
dim(enrichment_results)
# head(enrichment_results)

df <- enrichment_results

df_sig <- df %>%
  filter(NES > 0, FDR_q_val < fdr_threshold)  # keep only positive, significant terms

top_terms <- df_sig %>%
  group_by(name) %>%
  slice_max(order_by = NES, n = 1, with_ties = FALSE) %>%  # one row per group 
  pull(Term)

df_top <- df %>%
  filter(Term %in% top_terms) 

# wide matrix: rows = Term, cols = name, values = NES
mat_df <- df_top %>%
  select(name, Term, NES) %>%
  pivot_wider(names_from = name, values_from = NES, values_fill = 0)

mat <- mat_df %>% column_to_rownames("Term") %>% as.matrix()

# hierarchical clustering
row_dendro <- as.dendrogram(hclust(dist(mat), method = "ward.D2"))
col_dendro <- as.dendrogram(hclust(dist(t(mat)), method = "ward.D2"))

row_order <- order.dendrogram(row_dendro)  # get row order 
col_order <- order.dendrogram(col_dendro)  # get col order 

# prepare heatmap
heatmap_df <- mat_df %>%
  pivot_longer(-Term, names_to = "name", values_to = "NES") %>%
  left_join(df %>% select(name, Term, FDR_q_val), by = c("name", "Term")) %>%
  mutate(
    Term = factor(Term, levels = rownames(mat)[row_order]),
    name = factor(name, levels = colnames(mat)[col_order]),
    sig = FDR_q_val < fdr_threshold
  )

# make plot
enrichment_plot <- ggplot(heatmap_df, aes(x = name, y = Term, fill = NES)) +
  geom_tile() +
  # add star for significance
  geom_text(aes(label = ifelse(sig, "*", "")), vjust = 0.5) +
  # red-white-blue gradient centered at 0
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  # ensure square tiles
  coord_fixed() +
  # wrap long y-axis labels after 25 chars
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  # rotate x-axis labels
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank()
  ) + MrBiomics_theme()

# show plot
width <- 4
height <- 4
options(repr.plot.width=width, repr.plot.height=height)

# enrichment_plot

# save plot
ggsave_new(path=enrichment_plot_path,
           plot=enrichment_plot,
           width=width,
           height=height
          )


######### Lineage reconstruction using crossprediction ############
# plot adapted from: https://gist.github.com/dsparks/4331058

# load adjacency matrix
adj_mtx <- data.frame(fread(file.path(crossprediction_adj_mtx_path), header=TRUE), row.names=1)
dim(adj_mtx)
# head(adj_mtx)

# load previous best layout coordinates (TODO)
# coordinates <- as.matrix(data.frame(fread(file.path(crossprediction_adj_mtx_path), header=TRUE), row.names=1))
# dim(coordinates)
# head(coordinates)

# make cross prediction graph plot

# set colors & shapes
nodes <- rownames(adj_mtx)
# node_cond <- rep(cond, times = dim(adj_mtx)[1])
node_colors <- "#F8766D" #KO_col[nodes]
node_shape <- 19 # cond_shapes[node_cond]
border_color <- "#000000" #cond_colors[node_cond]
stroke_size <- 1.5
point_size <- 5

# data
adjacencyMatrix <- as.matrix(adj_mtx)
# parameters
curved=FALSE

adjacencyMatrix[adjacencyMatrix<cut_off] <- 0

# plot graph

# Empty ggplot2 theme
new_theme_empty <- theme_bw() 
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, 0, 0), unit = "lines",
                                         valid.unit = 3L, class = "unit")

if (exists("coordinates")){
    layoutCoordinates <- coordinates
} else{
    pdf(file = if (.Platform$OS.type == "windows") "NUL" else "/dev/null")  # route to null device (ie trash)
    layoutCoordinates <- gplot(adjacencyMatrix, mode="fruchtermanreingold")    # Get graph layout coordinates
    dev.off()
}


adjacencyList <- reshape2::melt(adjacencyMatrix)  # Convert to list of ties only
adjacencyList <- adjacencyList[adjacencyList$value >= cut_off, ] # prune weak edges


# Function to generate paths between each connected node
edgeMaker <- function(whichRow, len = 100, curved = TRUE){
    fromC <- layoutCoordinates[adjacencyList[whichRow, 1], ]  # Origin
    toC <- layoutCoordinates[adjacencyList[whichRow, 2], ]  # Terminus

    # Add curve:
    graphCenter <- colMeans(layoutCoordinates)  # Center of the overall graph
    bezierMid <- c(fromC[1], toC[2])  # A midpoint, for bended edges
    distance1 <- sum((graphCenter - bezierMid)^2)
    if(distance1 < sum((graphCenter - c(toC[1], fromC[2]))^2)){
    bezierMid <- c(toC[1], fromC[2])
    }  # To select the best Bezier midpoint
    bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
    if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve

    edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
    edge$Sequence <- 1:len  # For size and colour weighting in plot
    edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")

    # add a vector 'Probability' that linearly interpolates between the weights of the current edge and its counterpart
    x <- c(1,len)
    y <- c(adjacencyMatrix[adjacencyList[whichRow,'Var1'],adjacencyList[whichRow,'Var2']],adjacencyMatrix[adjacencyList[whichRow,'Var2'],adjacencyList[whichRow,'Var1']])
    edge$Probability <- approx(x,y,xout=1:len)$y

  return(edge)
}

# Generate a (curved) edge path for each pair of connected nodes
allEdges <- lapply(1:nrow(adjacencyList), edgeMaker, len = 500, curved = curved)
allEdges <- do.call(rbind, allEdges)  # a fine-grained path ^, with bend ^

zp1 <- ggplot(allEdges) 
zp1 <- zp1 + geom_path(aes(x = x, y = y, group = Group,  # Edges with gradient
                           colour = Sequence, size = Probability))  # and taper

zp1 <- zp1 + geom_point(data = data.frame(layoutCoordinates),  # Add nodes
#                         aes(x = x, y = y),  pch = 21, size = point_size, fill = node_colors, colour = border_color, stroke = stroke_size) # fill & borders are informative
                        aes(x = x, y = y), shape=node_shape, size = point_size, color = node_colors)#+scale_shape_manual(name = "condition", labels = names(cond_shapes),values = cond_shapes) # fill and shape are informative

# zp1 <- zp1 + geom_text(data = data.frame(layoutCoordinates), aes(x = x, y = y, label=nodes),#rownames(adjacencyMatrix)),
#                        hjust=0.5, vjust=-1)

zp1 <- zp1 + geom_label(data = data.frame(layoutCoordinates), aes(x = x, y = y, label=nodes), fill=node_colors,#rownames(adjacencyMatrix)),
                        hjust=0.5, vjust=-0.5)

zp1 <- zp1 + scale_colour_gradient(low = gray(0), high = gray(0), guide = "none") # Customize gradient

zp1 <- zp1 + scale_size(range = c(1/10, point_size-1))#, guide = "none")  # Customize taper

crosspred_p <- zp1 + 
        new_theme_empty + 
        #theme(legend.position = c(0.5, 0.5))+ 
        guides(size = guide_legend(title='Average\ncross-prediction\nprobability')) # Clean up plot


# plot the graph
width <- 6
height <- 6
options(repr.plot.width=width, repr.plot.height=height)

# crosspred_p


# # save layout coordinates
# write.csv(layoutCoordinates, file.path("results/KO150/KO_classifier",cond,'layoutCoordinates.csv'))

# save plot
ggsave_new(path=crossprediction_plot_path,
           plot=crosspred_p,
           width=width,
           height=height
          )


######### patchwork of sub figure ############
