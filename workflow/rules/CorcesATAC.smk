# ATAC-seq Analysis Recipe applied to healthy hematopoeitic ATAC-seq samples 
# from Corces et al. 2016 Nature Genetics (https://www.nature.com/articles/ng.3646)
# Subset of GEO SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75384
# Subset of GEO SuperSeries for ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74912

#### CorcesATAC - Fetch NGS ####
module CorcesATAC_fetch_ngs:
    snakefile:
        github("epigen/fetch_ngs", path="workflow/Snakefile", tag="v1.0.5")
    config:
        config_wf["CorcesATAC_fetch_ngs"]

use rule * from CorcesATAC_fetch_ngs as CorcesATAC_fetch_ngs_*

#### CorcesATAC - Get resources from Zenodo (custom helper rule) ####
# Downloads Bowtie2 indices for hg38 from Zenodo record 6344173 and unpacks them.
rule CorcesATAC_get_resources:
    output:
        "resources/CorcesATAC/atacseq_pipeline/hg38/gencode.v38.basic.annotation.gtf",
        resource_dir = directory("resources/CorcesATAC/atacseq_pipeline/hg38/"),
    params:
        zenodo_record = "6344173",
        zip_filename = "indices_for_Bowtie2.zip"
    conda:
        "../envs/zenodo_get.yaml"
    shell:
        """
        # Download the specific record to the target directory
        zenodo_get --record {params.zenodo_record} --output-dir={output.resource_dir}

        # Change directory, unzip the specific file, and remove the zip archive
        # Using && ensures commands run sequentially and stop if one fails
        cd {output.resource_dir} && \
        unzip {params.zip_filename} && \
        rm {params.zip_filename}
        """

### CorcesATAC - ATAC-seq processing ####
module CorcesATAC_atacseq_pipeline:
    snakefile:
        github("epigen/atacseq_pipeline", path="workflow/Snakefile", tag="v2.1.0")
    config:
        config_wf["CorcesATAC_atacseq_pipeline"]

use rule * from CorcesATAC_atacseq_pipeline as CorcesATAC_atacseq_pipeline_*

#### CorcesATAC - Genome Tracks #### 
module CorcesATAC_genome_tracks:
    snakefile:
        github("epigen/genome_tracks", path="workflow/Snakefile", tag="v2.0.5")
    config:
        config_wf["CorcesATAC_genome_tracks"]

use rule * from CorcesATAC_genome_tracks as CorcesATAC_genome_tracks_*

#### CorcesATAC - Spilterlize & Integrate #### 
module CorcesATAC_spilterlize_integrate:
    snakefile:
        github("epigen/spilterlize_integrate", path="workflow/Snakefile", tag="v3.0.5")
    config:
        config_wf["CorcesATAC_spilterlize_integrate"]

use rule * from CorcesATAC_spilterlize_integrate as CorcesATAC_spilterlize_integrate_*

### CorcesATAC - Unsupervised Analysis ####
module CorcesATAC_unsupervised_analysis:
    snakefile:
        github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.3")
    config:
        config_wf["CorcesATAC_unsupervised_analysis"]

use rule * from CorcesATAC_unsupervised_analysis as CorcesATAC_unsupervised_analysis_*

#### CorcesATAC - Differential Accessibility Analysis #### 
module CorcesATAC_dea_limma:
    snakefile:
        github("epigen/dea_limma", path="workflow/Snakefile", tag="v2.2.1")
    config:
        config_wf["CorcesATAC_dea_limma"]

use rule * from CorcesATAC_dea_limma as CorcesATAC_dea_limma_*

#### CorcesATAC - Feature lists to BED files (custom helper rule) ####
rule CorcesATAC_convert_features2bed:
    input:
        consensus_annotation = os.path.join("results/CorcesATAC/atacseq_pipeline/counts/consensus_annotation.csv"),
        features_txt = os.path.join("results/CorcesATAC/dea_limma/{analysis}/feature_lists/{feature_set}_features.txt"),
    output:
        features_bed = os.path.join("results/CorcesATAC/dea_limma/{analysis}/feature_lists/{feature_set}_features.bed"),
    params:
        region_col = "peak_id",
        chr_col = "gencode_chr",
        start_col = "gencode_start",
        end_col = "gencode_end",
    resources:
        mem_mb="4000",
    log:
        os.path.join("logs","rules","CorcesATAC_convert_features2bed_{analysis}_{feature_set}.log"),
    run:
        # load files as pandas df
        consensus_df = pd.read_csv(input.consensus_annotation)
        features_df = pd.read_csv(input.features_txt, header=None, names=[params.region_col])
        # map using params
        merged_df = pd.merge(features_df, consensus_df, on=params.region_col, how="inner")
        # Select and order the columns required for the BED file format.
        bed_df = merged_df[[params.chr_col, params.start_col, params.end_col, params.region_col]]
        # save in BED format
        bed_df.to_csv(output.features_bed, sep="\t", header=False, index=False)

#### CorcesATAC - Enrichment Analysis ####
module CorcesATAC_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["CorcesATAC_enrichment_analysis"]

use rule * from CorcesATAC_enrichment_analysis as CorcesATAC_enrichment_analysis_*

#### CorcesATAC - Lineage Reconstruction (custom rule) ####
rule CorcesATAC_reconstruct_lineage:
    input:
        data = os.path.join("results/CorcesATAC/spilterlize_integrate/all/normCQN_integrated_HVF.csv"),
        metadata = os.path.join("results/CorcesATAC/spilterlize_integrate/all/annotation.csv"),
        feature_annotation = os.path.join("results/CorcesATAC/atacseq_pipeline/counts/consensus_annotation.csv"),
    output:
        adjacency_matrix = os.path.join("results/CorcesATAC/special_analyses/crossprediction/adjacency_matrix.csv"),
        top_features = os.path.join("results/CorcesATAC/special_analyses/crossprediction/top_features.csv"),
        graph = os.path.join("results/CorcesATAC/special_analyses/crossprediction/graph.png"),
    params:
        group_var = "cell_type",
        group_rm = "",
        top_features_n = 5,
        prune_th = 0.2,
        feature_annotation_var = "homer_Gene_Name",
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","CorcesATAC_reconstruct_lineage.log"),
    script:
        "../scripts/crossprediction.py"
