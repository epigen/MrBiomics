# Quickstart: Enrichment analysis applied to a sub-selection of genomic region files resulting 
# from a One-vs-All differential accessibility analysis of ATAC-seq data of healthy hematopoietic 
# cells from Corces et al. 2016 Nature Genetics (https://www.nature.com/articles/ng.3646).
# It focuses on regions that are significantly more open in four specific cell types 
# (Bcell, Ery, Mono, CD8Tcell) compared to the rest of the hematopoietic cell types.
# For details on the dataset and the analysis, see the Quickstart and the ATAC-seq recipe on the wiki.

### Quickstart - Enrichment Analysis ####
module quickstart_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="main")
    config:
        config_wf["quickstart_enrichment_analysis"]

use rule * from quickstart_enrichment_analysis as quickstart_enrichment_analysis_*