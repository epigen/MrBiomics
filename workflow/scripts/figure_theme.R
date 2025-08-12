### LIBRARIES
# select libraries
required_libs <- c(
    "ggplot2",
    "ggalign",
    "dendsort",
    "tidyr",
    "tibble",
    "patchwork",
    "scales",
    "dplyr",
    "reshape2",
    "ggrepel",
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
options(stringsAsFactors=F)
                      
# MrBiomics plotting theme (TODO)
MrBiomics_theme <- function(){
    
    # settings
    font <- "Arial"
    size <- 12
    
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
      
#       axis.text.x = element_text(            #margin for axis text
#                     margin=margin(5, b = 10))
    )
}

### FUNCTIONS

# extended ggsave
ggsave_all_formats <- function(path, plot, width=5, height=5, dpi=300){

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
# cell type colors
CELL_TYPE_COLORS <- c(
    'HSC'='#566E21',   # HSC
    'MPP'='#96C03A',  # MPP
    'LMPP'='#7FC493',  # LMPP
    'CMP'='#FBCA03',  # CMP
    'GMP'='#B97D10',  # GMP
    'MEP'='#E95E30',  # MEP
    'Mono'='#AA0505',  # Mono
    'Ery'='#6A0C0B',  # Ery
    'CLP'='#67C7EB',  # CLP
    'B'='#977CBA',  # B
    'CD4'='#3B8DAC',  # CD4
    'CD8'='#0E536C',  # CD8
    'NK'='#C73188'  # NK
)

# blue green violet: https://coolors.co/0e536c-3b8dac-67c7eb-7fc493-96c03a-566e21-977cba-c73188
# red yellow orange: https://coolors.co/e95e30-aa0505-6a0c0b-b97d10-fbca03

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