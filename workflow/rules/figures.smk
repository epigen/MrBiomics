# Snakemake file for creating the figure for the manuscript
#### CorcesRNA - Figure panels (custom rule) ####
rule CorcesRNA_figures:
    input:
        enrichment_results = os.path.join("results/CorcesRNA/enrichment_analysis/cell_types/preranked_GSEApy/Azimuth_2023/cell_types_Azimuth_2023_all.csv"),
        crossprediction_adj_mtx = os.path.join("results/CorcesRNA/special_analyses/crossprediction/adjacency_matrix.csv"),
    output:
        enrichment_plot = os.path.join("docs/CorcesRNA/enrichment_analysis.pdf"),
        crossprediction_plot = os.path.join("docs/CorcesRNA/crossprediction_plot.pdf"),
    params:
        figure_theme_path = workflow.source_path("../scripts/figure_theme.R"),
        # enrichment analysis
        fdr_threshold = 0.05,
        # lineage reconstructions
        cut_off = 0.2 # when pruning at 0.2 both HVF and integrated produce the same graph
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","CorcesRNA_figures.log"),
    script:
        "../scripts/CorcesRNA_figures.R"
