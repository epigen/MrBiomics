# digits Analysis

### digits - Load data with custom rule ####
rule load_digits:
    output:
        data = os.path.join('data','digits','digits_data.csv'),
        labels = os.path.join('data','digits','digits_labels.csv'),
    resources:
        mem_mb=1000,
    threads: 1
    conda:
        "../envs/sklearn.yaml"
    log:
        os.path.join("logs","rules","load_digits.log"),
    script:
        "../scripts/digits/load_digits.py"

### digits - Unsupervised Analysis ####
module digits_unsupervised_analysis:
    snakefile:
        "/research/home/sreichl/projects/unsupervised_analysis/workflow/Snakefile" # absolute path to my current clone
        # github("epigen/unsupervised_analysis", path="workflow/Snakefile", tag="v3.0.0") # Snakemake issue: https://github.com/snakemake/snakemake/issues/3212
    config:
        config_wf["digits_unsupervised_analysis"]

use rule * from digits_unsupervised_analysis as digits_unsupervised_analysis_*
