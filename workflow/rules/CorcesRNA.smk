# RNA-seq Analysis Recipe for the healthy hematopoeitic RNA-seq samples 
# from Corces et al. 2016 Nature Genetics (https://www.nature.com/articles/ng.3646)
# Subset of GEO SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75384
# Subset of GEO SuperSeries for **unstranded** RNA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246

### CorcesRNA - Fetch NGS (optional) ####
module CorcesRNA_fetch_ngs:
    snakefile:
       github("epigen/fetch_ngs", path="workflow/Snakefile", tag="v1.0.5")
    config:
        config_wf["CorcesRNA_fetch_ngs"]

use rule * from CorcesRNA_fetch_ngs as CorcesRNA_fetch_ngs_*

### CorcesRNA - RNA-seq processing ####
module CorcesRNA_rnaseq_pipeline:
    snakefile:
        # "/home/sreichl/projects/rnaseq_pipeline/workflow/Snakefile"
        github("epigen/rnaseq_pipeline", path="workflow/Snakefile", tag="v1.1.1")
    config:
        config_wf["CorcesRNA_rnaseq_pipeline"]

use rule * from CorcesRNA_rnaseq_pipeline as CorcesRNA_rnaseq_pipeline_*

#### CorcesRNA - Genome Tracks #### 
module CorcesRNA_genome_tracks:
    snakefile:
        github("epigen/genome_tracks", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["CorcesRNA_genome_tracks"]

use rule * from CorcesRNA_genome_tracks as CorcesRNA_genome_tracks_*

#### CorcesRNA - Spilterlize & Integrate #### 
module CorcesRNA_spilterlize_integrate:
    snakefile:
        # "/home/sreichl/projects/spilterlize_integrate/workflow/Snakefile"
        github("epigen/spilterlize_integrate", path="workflow/Snakefile", tag="v3.0.1")
    config:
        config_wf["CorcesRNA_spilterlize_integrate"]

use rule * from CorcesRNA_spilterlize_integrate as CorcesRNA_spilterlize_integrate_*

### CorcesRNA - Unsupervised Analysis ####
module CorcesRNA_unsupervised_analysis:
    snakefile:
        github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.1")
    config:
        config_wf["CorcesRNA_unsupervised_analysis"]

use rule * from CorcesRNA_unsupervised_analysis as CorcesRNA_unsupervised_analysis_*

#### CorcesRNA - Differential Expression Analysis #### 
module CorcesRNA_dea_limma:
    snakefile:
        # "/home/sreichl/projects/dea_limma/workflow/Snakefile"
        # github("epigen/dea_limma", path="workflow/Snakefile", tag="v2.1.1")
        github("epigen/dea_limma", path="workflow/Snakefile", branch="main")
    config:
        config_wf["CorcesRNA_dea_limma"]

use rule * from CorcesRNA_dea_limma as CorcesRNA_dea_limma_*

#### CorcesRNA - Enrichment Analysis ####
module CorcesRNA_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.2")
    config:
        config_wf["CorcesRNA_enrichment_analysis"]

use rule * from CorcesRNA_enrichment_analysis as CorcesRNA_enrichment_analysis_*

#### CorcesRNA - Lineage Reconstruction (custom rule) ####
rule CorcesRNA_reconstruct_lineage:
    input:
        data = os.path.join("results/CorcesRNA/spilterlize_integrate/all/normCQN_integrated_HVF.csv"),
        metadata = os.path.join("results/CorcesRNA/spilterlize_integrate/all/annotation.csv"),
    output:
        adjacency_matrix = os.path.join("results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"),
        top_features = os.path.join("results/CorcesRNA/special_analyses/crossprediction/top_features.csv"),
    params:
        group_var = "cell_type",
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","CorcesRNA_reconstruct_lineage.log"),
    script:
        "../scripts/crossprediction.py"

#### CorcesRNA - Copy selected plots for documentation and visualization in wiki (custom rule) ####
rule CorcesRNA_plots:
    input:
        genome_tracks_MS4A1  = os.path.join("results/CorcesRNA/genome_tracks/tracks/MS4A1.pdf"),
        spilterlize_integrate_filtered = os.path.join("results/CorcesRNA/spilterlize_integrate/all/plots/filtered.png"),
        spilterlize_integrate_normCQN = os.path.join("results/CorcesRNA/spilterlize_integrate/all/plots/normCQN.png"),
        spilterlize_integrate_normCQN_integrated = os.path.join("results/CorcesRNA/spilterlize_integrate/all/plots/normCQN_integrated.png"),
        enrichment_analysis_summary = os.path.join("results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_summary.png"),
    output:
        genome_tracks_MS4A1  = os.path.join("docs/CorcesRNA/MS4A1.pdf"),
        spilterlize_integrate_filtered = os.path.join("docs/CorcesRNA/filtered.png"),
        spilterlize_integrate_normCQN = os.path.join("docs/CorcesRNA/normCQN.png"),
        spilterlize_integrate_normCQN_integrated = os.path.join("docs/CorcesRNA/normCQN_integrated.png"),
        enrichment_analysis_summary = os.path.join("docs/CorcesRNA/cell_types_Azimuth_2023_summary.png"),
    resources:
        mem_mb="1000",
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","CorcesRNA_plots.log"),
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        cp {input[2]} {output[2]}
        cp {input[3]} {output[3]}
        cp {input[4]} {output[4]}
        """


