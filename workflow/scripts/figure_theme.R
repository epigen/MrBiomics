

########## TODOs
# - this is semi-cleaned and adapted. keep what is useful and remove the rest. this is just a starting point
# - optimize for downstream Adobe Illustrator/Inkscape usage ideal (e.g., always produce PDF for figures and PNG for slides)
# - ideally this can become a plotting template for future papers

### LIBRARIES
# select libraries
required_libs <- c(
    "ggplot2",
    "tidyr",
    "tibble",
    "patchwork",
    "scales",
    "dplyr",
    "reshape2",
    "pheatmap",
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
    size <- 6
    
    theme_bw(
        base_size=size,
        base_family = font
    ) %+replace% 
    
    theme(
      #grid elements
#       panel.grid.major = element_blank(),    #strip major gridlines
#       panel.grid.minor = element_blank(),    #strip minor gridlines
#       axis.ticks = element_blank(),          #strip axis ticks
      
#       strips axis lines ?
      
      #text elements
        text = element_text(              
                   family = font,           
                   size = size),
        
      plot.title = element_text(             #title
                   family = font,            #set font family
                   size = size,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0.5,                #center align
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
ggsave_new <- function(path, plot, width=5, height=5, dpi=300){

    # close any lingering devices before opening a new one
    while (!is.null(dev.list())) dev.off()
    
    # get paths
    dir_path <- file.path(dirname(path))
    filename <- sub("\\.[^.]+$", "", basename(path))
    
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
                      
## CorcesRNA & CorcesATAC
# cell type colors
celltype_colors <- c('HSC'='#707070')

# cell type shapes
celltype_shapes <- c('HSC' = 16)
