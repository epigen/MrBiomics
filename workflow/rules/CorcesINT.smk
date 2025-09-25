# Integrative Analysis Recipe applied to matched(!) healthy hematopoeitic RNA-seq & ATAC-seq samples
# from Corces et al. 2016 Nature Genetics (https://www.nature.com/articles/ng.3646)
# Subset of GEO SuperSeries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75384
# Subset of GEO SuperSeries for **unstranded** RNA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246
# Subset of GEO SuperSeries for ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74912

#### CorcesINT - Merge count matrices (custom rule) #### 
rule CorcesINT_merge_counts:
    input:
        RNA_counts = os.path.join("results/CorcesRNA/rnaseq_pipeline/counts/counts.csv"),
        # here you can choose between promoter or TSS counts (difference explained in ATAC-seq pipeline module docs)
        ATAC_counts = os.path.join("results/CorcesATAC/atacseq_pipeline/counts/promoter_counts.csv"),
        # ATAC_counts = os.path.join("results/CorcesATAC/atacseq_pipeline/counts/TSS_counts.csv"),
    output:
        INT_counts = os.path.join("results/CorcesINT/counts.csv"),
    resources:
        mem_mb="4000",
    threads: config.get("threads", 1)
    log:
        "logs/rules/CorcesINT_merge_counts.log",
    run:
        # load data
        rna = pd.read_csv(input.RNA_counts, index_col=0).add_suffix('_RNA')
        atac = pd.read_csv(input.ATAC_counts, index_col=0).add_suffix('_ATAC')
        # merge data
        merged = rna.merge(atac, left_index=True, right_index=True, how='inner')
        # save merged data
        merged.to_csv(output.INT_counts)

#### CorcesINT - Spilterlize & Integrate #### 
module CorcesINT_spilterlize_integrate:
    snakefile:
        github("epigen/spilterlize_integrate", path="workflow/Snakefile", tag="v3.0.4")
    config:
        config_wf["CorcesINT_spilterlize_integrate"]

use rule * from CorcesINT_spilterlize_integrate as CorcesINT_spilterlize_integrate_*

### CorcesINT - Unsupervised Analysis ####
module CorcesINT_unsupervised_analysis:
    snakefile:
        github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.3")
    config:
        config_wf["CorcesINT_unsupervised_analysis"]

use rule * from CorcesINT_unsupervised_analysis as CorcesINT_unsupervised_analysis_*

#### CorcesINT - Differential Expression Analysis #### 
module CorcesINT_dea_limma:
    snakefile:
        github("epigen/dea_limma", path="workflow/Snakefile", tag="v2.2.1")
    config:
        config_wf["CorcesINT_dea_limma"]

use rule * from CorcesINT_dea_limma as CorcesINT_dea_limma_*

# #### CorcesINT - Visualize correlation between RNA-seq & ATAC-seq (custom rule) ####
rule CorcesINT_plot_correlation:
    input:
        data = os.path.join("results/CorcesINT/spilterlize_integrate/all/normupperquartile_integrated.csv"),
        metadata = os.path.join("results/CorcesINT/spilterlize_integrate/all/annotation.csv"),
        dea_results = os.path.join("results/CorcesINT/dea_limma/normupperquartile_integrated/results.csv"),
    output:
        correlation_plots = expand(os.path.join("results/CorcesINT/special_analysis/correlation_plots/{group}_correlation.png"), group=['CMP','LMPP','MEP','HSC','CD4Tcell','CD8Tcell','GMP','Mono','CLP','MPP','Ery','NKcell','Bcell']),
    params:
        adjp_th = config_wf["CorcesINT_dea_limma"]["filters"]["adj_pval"],
        lfc_th = config_wf["CorcesINT_dea_limma"]["filters"]["lfc"],
        ave_expr_th = config_wf["CorcesINT_dea_limma"]["filters"]["ave_expr"],
    log:
        "logs/rules/CorcesINT_plot_correlation.log",
    resources:
        mem_mb="8000",
    threads: config.get("threads", 1)
    
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/plot_integration_correlation.R"

# #### CorcesINT - Download Enrichment Analysis Resources (custom rules) ####
rule CorcesINT_download_motif_annotation:
    output:
        "resources/CorcesINT/enrichment_analysis/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    params:
        url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    log:
        "logs/aria2c/CorcesINT_download_motif_annotation.log",
    wrapper:
        "v7.2.0/utils/aria2c"

rule CorcesINT_download_rankings_database:
    output:
        "resources/CorcesINT/enrichment_analysis/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    params:
        url="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    log:
        "logs/aria2c/CorcesINT_download_rankings_database.log",
    wrapper:
        "v7.2.0/utils/aria2c"

#### CorcesINT - Enrichment Analysis ####
module CorcesINT_enrichment_analysis:
    snakefile:
        github("epigen/enrichment_analysis", path="workflow/Snakefile", tag="v2.0.3")
    config:
        config_wf["CorcesINT_enrichment_analysis"]

use rule * from CorcesINT_enrichment_analysis as CorcesINT_enrichment_analysis_*
