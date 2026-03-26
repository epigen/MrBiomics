# Quick Start example:
# Run the enrichment_analysis module on the grouped quickstart BED sets
# defined in config/quickstart/quickstart_enrichment_analysis_annotation.csv.

module quickstart_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="main")
    config:
        config_wf["quickstart_enrichment_analysis"]

use rule * from quickstart_enrichment_analysis as quickstart_enrichment_analysis_*


