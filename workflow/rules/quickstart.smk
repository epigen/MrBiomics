# TODO: Add short description of this analysis here, see other rule files. 

### quickstart - Download Enrichment Analysis Resources (custom rules) ####


### quickstart - Enrichment Analysis ####
module quickstart_enrichment_analysis:
    snakefile: 
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.1")
    config: 
        config_wf["quickstart_enrichment_analysis"]

use rule * from quickstart_enrichment_analysis as quickstart_enrichment_analysis_*