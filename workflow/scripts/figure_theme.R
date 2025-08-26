### LIBRARIES
# select libraries
required_libs <- c(
    "ggplot2",
    "ggdendro",
    "dendsort",
    "tidyr",
    "tibble",
    "patchwork",
    "scales",
    "dplyr",
    "reshape2",
    "ggrepel",
    "ggrastr",
    "stringr",
    "gridGraphics",
    "ggplotify",
    "directlabels",
    "data.table",
    "svglite",
    "sna",
    "RColorBrewer",
    "Hmisc"
)

set.seed(42)

# load libraries
loaded_libs <- lapply(required_libs, function(x) suppressWarnings(suppressMessages(library(x, character.only = TRUE))))
options(stringsAsFactors=F, ggrastr.default.dpi = 500)

# in cm, for A4 paper
PAGE_WIDTH <- 21
PLOT_SIZE_5_PER_ROW <- PAGE_WIDTH/5
PLOT_SIZE_4_PER_ROW <- PAGE_WIDTH/4
PLOT_SIZE_3_PER_ROW <- PAGE_WIDTH/3
PLOT_SIZE_2_PER_ROW <- PAGE_WIDTH/2

# as per Nature and Science guidelines
FONT <- "Arial"
FONT_SIZE_NORMAL <- 7
FONT_SIZE_SMALL <- 5

# MrBiomics plotting theme (TODO)
MrBiomics_theme <- function(){
    
    # settings
    font <- FONT
    size <- FONT_SIZE_NORMAL
    
    theme_bw(
        base_size=size,
        base_family = font
    ) %+replace% 
    
    theme(
        ###  grid elements
        axis.ticks = element_blank(),          #remove axis ticks
        ### text elements
        text = element_text(              
                    family = font,           
                    size = size),
            
        plot.title = element_text(             #title
                    family = font,            #set font family
                    size = size,                #set font size
                    face = 'bold',            #bold typeface
                    hjust = 0,                #center align
                    vjust = 2),               #raise slightly
        
        plot.subtitle = element_text(          #subtitle
                    family = font,            #font family
                    size = size),               #font size
        
        plot.caption = element_text(           #caption
                    family = font,            #font family
                    size = size,                 #font size
                    hjust = 0.5),               #center align
        
        axis.title = element_text(             #axis titles
                    family = font,            #font family
                    size = size),               #font size
        
        axis.text = element_text(              #axis text
                    family = font,            #axis famuly
                    size = size),                #font size
            
        legend.text = element_text(              #axis text
                    family = font,            #axis famuly
                    size = size), 
        
        legend.title = element_text(              #axis text
                    family = font,            #axis famuly
                    size = size),
        
        # compact legends globally
        legend.key.height = grid::unit(0.4, "lines"),
        legend.key.width  = grid::unit(0.6, "lines"),
        legend.spacing.y  = grid::unit(0.1, "lines"),
        legend.spacing    = grid::unit(0.2, "lines"),
        legend.box.spacing= grid::unit(0.2, "lines")
      
#       axis.text.x = element_text(            #margin for axis text
#                     margin=margin(5, b = 10))
    )
}

# text in axis is in pt, but geom_text etc use other units --> fix here
# note: .pt is a conversion factor from mm to pt (setting size.units = "pt" did not work)
update_geom_defaults("text", list(size = FONT_SIZE_SMALL / .pt, family = FONT))
update_geom_defaults("text_repel", list(size = FONT_SIZE_SMALL / .pt, segment.size = 0.2, family = FONT))
update_geom_defaults("label", list(size = FONT_SIZE_SMALL / .pt, family = FONT))

# same theme as MrBiomics_theme but everything removed except for the title
MrBiomics_void <- function(){
    font <- FONT
    size <- FONT_SIZE_NORMAL
    
    theme_void() +
        theme(
            plot.title = element_text(hjust = 0, size = size, family = font, face = "bold", vjust = 2),
            legend.text = element_text(              #axis text
                        family = font,            #axis famuly
                        size = size), 
            
            legend.title = element_text(              #axis text
                        family = font,            #axis famuly
                        size = size),
            
            # compact legends globally
            legend.key.height = grid::unit(0.4, "lines"),
            legend.key.width  = grid::unit(0.6, "lines"),
            legend.spacing.y  = grid::unit(0.1, "lines"),
            legend.spacing    = grid::unit(0.2, "lines"),
            legend.box.spacing= grid::unit(0.2, "lines")
        )
}

### FUNCTIONS

# extended ggsave
ggsave_all_formats <- function(path, plot, width=PLOT_SIZE_3_PER_ROW, height=PLOT_SIZE_3_PER_ROW, dpi=300){

    # close any lingering devices before opening a new one
    while (!is.null(dev.list())) dev.off()
    
    # get paths
    dir_path <- file.path(dirname(path))
    filename <- sub("\\.[^.]+$", "", basename(path))
    
    # ensure directory exists
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    
    # SVG, PNG & PDF
    for (format in c('svg','png','pdf')){
        ggsave(
          paste0(filename,'.',format),
          plot = plot,
          device = if (format == "pdf") cairo_pdf else format,
          path = dir_path,
          scale = 1,
          dpi = dpi,
            width = width,
            height = height,
          limitsize = FALSE,
          units = "cm"
        )
    }
    print(paste0("Saved ", path))
}
                      
remove_term_suffix <- function(db, terms){
    
    if(grepl('GO', db, fixed = TRUE)){
        # remove "(GO:.......)" from terms & abbreviate if necessary
        return(gsub("GO:.......", "", gsub("\\(GO:.......)", "", terms)))
    }
    
    if(grepl('WikiPathways', db, fixed = TRUE)){
        # remove WikiPathway IDs WP+numbers
        return(gsub("WP.*", "", terms))
    }
    
    return(terms)
}

                      
### DEFINITIONS (e.g., shapes, colors,...)
IRONMAN_COLORS <- c(
    "red"="#AA0505",
    "darkred"="#6A0C0B",
    "gold"="#B97D10",
    "yellow"="#FBCA03",
    "lightblue"="#67C7EB"
)

## CorcesRNA & CorcesATAC
# blue green violet: https://coolors.co/0e536c-3b8dac-51aacc-67c7eb-977cba-c73188-7fc493-96c03a-86ac35-566e21
# red yellow orange: https://coolors.co/e95e30-ca321b-aa0505-b97d10-daa40a-fbca03-e1de3d-96c03a

# cell type colors
CELL_TYPE_COLORS <- c(
    'HSC'='#566E21',
    'MPP'='#96C03A',
    'LMPP'='#7FC493',
    'CMP'='#E1DE3D',
    'GMP'='#FBCA03',
    'Mono'='#B97D10',
    'CLP'='#67C7EB',
    'CD4'='#3B8DAC',
    'CD8'='#0E536C',
    'B'='#977CBA',
    'NK'='#C73188',
    'MEP'='#E95E30',
    'Ery'='#AA0505'
)

# Map from data names to CELL_TYPE_COLORS names
DATA_TO_CELL_TYPE_COLORS_MAPPING <- c(
    "Bcell" = "B",
    "CD4Tcell" = "CD4", 
    "CD8Tcell" = "CD8",
    "CLP" = "CLP",
    "CMP" = "CMP",
    "Ery" = "Ery",
    "GMP" = "GMP",
    "HSC" = "HSC",
    "LMPP" = "LMPP",
    "MEP" = "MEP",
    "MPP" = "MPP",
    "Mono" = "Mono",
    "NKcell" = "NK"
)

CELL_TYPE_TO_LINEAGE_MAPPING <- c(
    "HSC" = "Progenitor",
    "MPP" = "Progenitor",
    "LMPP" = "Progenitor",
    "CMP" = "Progenitor",
    "GMP" = "Myeloid",
    "Mono" = "Myeloid",
    "MEP" = "Erythroid",
    "Ery" = "Erythroid",
    "CLP" = "Lymphoid",
    "CD4" = "Lymphoid",
    "CD8" = "Lymphoid",
    "NK" = "Lymphoid",
    "B" = "Lymphoid"
)

LINEAGE_COLORS <- c(
    "Progenitor" = "#86AC35",
    "Myeloid" = "#DAA40A",
    "Erythroid" = "#CA321B",
    "Lymphoid" = "#51AACC"
)

CELL_TYPE_AND_LINEAGE_COLORS <- c(CELL_TYPE_COLORS, LINEAGE_COLORS)

HAEMATOPOIESIS_MARKERS <- c(
    "CD34",
    "CD38",
    "THY1",
    "KIT",
    "FLT3",
    "IL3RA",
    "PTPRC",
    "CSF1R",
    "IL7R",
    "MME",
    "CD19",
    "MS4A1",
    "CD3E",
    "CD4",
    "CD8A",
    "CD8B",
    "NCAM1",
    "FCGR3A",
    "CD14",
    "HLA-DRA",
    "GYPA",
    "TFRC"
)

children_list <- list(
    HSC = c("MPP"),
    MPP = c("LMPP", "CMP"),
    LMPP = c("CLP", "GMP"),
    CMP = c("GMP", "MEP"),
    GMP = c("Mono"),
    MEP = c("Ery"),
    CLP = c("CD4", "CD8", "NK", "B")
)

# picked from colorpalette RdBu from ggplot2
RdBu_extremes <- c(
    "up" = "#B6242F",  # red
    "down" = "#2569AD"  # blue
)