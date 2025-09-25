# RNA-seq Analysis Recipe applied to healthy hematopoeitic RNA-seq samples 
# from Corces et al. 2016 Nature Genetics (https://www.nature.com/articles/ng.3646)
# Subset of GEO SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75384
# Subset of GEO SuperSeries for **unstranded** RNA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246

### CorcesRNA - Fetch NGS ####
module CorcesRNA_fetch_ngs:
    snakefile:
       github("epigen/fetch_ngs", path="workflow/Snakefile", tag="v1.0.5")
    config:
        config_wf["CorcesRNA_fetch_ngs"]

use rule * from CorcesRNA_fetch_ngs as CorcesRNA_fetch_ngs_*

### CorcesRNA - RNA-seq processing ####
module CorcesRNA_rnaseq_pipeline:
    snakefile:
        github("epigen/rnaseq_pipeline", path="workflow/Snakefile", tag="v1.2.1")
    config:
        config_wf["CorcesRNA_rnaseq_pipeline"]

use rule * from CorcesRNA_rnaseq_pipeline as CorcesRNA_rnaseq_pipeline_*

#### CorcesRNA - Genome Tracks #### 
module CorcesRNA_genome_tracks:
    snakefile:
        github("epigen/genome_tracks", path="workflow/Snakefile", tag="v2.0.5")
    config:
        config_wf["CorcesRNA_genome_tracks"]

use rule * from CorcesRNA_genome_tracks as CorcesRNA_genome_tracks_*

#### CorcesRNA - Spilterlize & Integrate #### 
module CorcesRNA_spilterlize_integrate:
    snakefile:
        github("epigen/spilterlize_integrate", path="workflow/Snakefile", tag="v3.0.4")
    config:
        config_wf["CorcesRNA_spilterlize_integrate"]

use rule * from CorcesRNA_spilterlize_integrate as CorcesRNA_spilterlize_integrate_*

### CorcesRNA - Unsupervised Analysis ####
module CorcesRNA_unsupervised_analysis:
    snakefile:
        github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.3")
    config:
        config_wf["CorcesRNA_unsupervised_analysis"]

use rule * from CorcesRNA_unsupervised_analysis as CorcesRNA_unsupervised_analysis_*

#### CorcesRNA - Differential Expression Analysis #### 
module CorcesRNA_dea_limma:
    snakefile:
        github("epigen/dea_limma", path="workflow/Snakefile", tag="v2.2.1")
    config:
        config_wf["CorcesRNA_dea_limma"]

use rule * from CorcesRNA_dea_limma as CorcesRNA_dea_limma_*

#### CorcesRNA - Enrichment Analysis ####
module CorcesRNA_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["CorcesRNA_enrichment_analysis"]

use rule * from CorcesRNA_enrichment_analysis as CorcesRNA_enrichment_analysis_*

#### CorcesRNA - Lineage Reconstruction (custom rule) ####
rule CorcesRNA_reconstruct_lineage:
    input:
        data = os.path.join("results/CorcesRNA/spilterlize_integrate/all/normCQN_integrated_HVF.csv"),
        metadata = os.path.join("results/CorcesRNA/spilterlize_integrate/all/annotation.csv"),
        feature_annotation = os.path.join("results/CorcesRNA/rnaseq_pipeline/counts/gene_annotation.csv"),
    output:
        adjacency_matrix = os.path.join("results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"),
        top_features = os.path.join("results/CorcesRNA/special_analyses/crossprediction/top_features.csv"),
        graph = os.path.join("results/CorcesRNA/special_analyses/crossprediction/graph.png"),
    params:
        group_var = "cell_type",
        group_rm = "",
        top_features_n = 5,
        prune_th = 0.2,
        feature_annotation_var = "external_gene_name",
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","CorcesRNA_reconstruct_lineage.log"),
    script:
        "../scripts/crossprediction.py"
