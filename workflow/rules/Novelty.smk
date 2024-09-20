from workflow_helpers import get_novelty_null_iters, get_novelty_years

rule count_journal_cooccurence_observed:
    input: ancient(rules.dl_gbq_journal_refs.output)
    output: OBSERVED_JOURNAL_COOCCURENCE_COUNTS
    conda: "../envs/python-minimal.yaml"
    resources:
        mem_mb=16000,
        runtime=45,
        cpus=8
    params: shuffle = False
    script: "../scripts/processing/count_journal_cooccurence.py"

rule count_journal_cooccurence_null:
    input: rules.dl_gbq_journal_refs.output
    output: NULL_JOURNAL_COOCCURENCE_COUNTS
    conda: "../envs/python-minimal.yaml"
    resources:
    resources:
        mem_mb=16000,
        runtime=45,
        cpus=8
    params: shuffle = True
    script: "../scripts/processing/count_journal_cooccurence.py"

rule calculate_novelty_zscores:
    input:
        observed = ancient(rules.count_journal_cooccurence_observed.output),
        null = lambda wc: collect(
                    rules.count_journal_cooccurence_null.output,
                    iter = get_novelty_null_iters(config),
                    year = wc.year
        )
    output: NOVELTY_ZSCORES
    conda: "../envs/python-minimal.yaml"
    resources:
    resources:
        mem_mb=16000,
        runtime=45,
        cpus=8
    script: "../scripts/processing/calculate_novelty_zscores.py"

rule agg_novelty_zscores:
    input:
        expand(
            rules.calculate_novelty_zscores.output, 
            year = get_novelty_years(config)
        )
    output: AGG_NOVELTY_ZSCORES
    conda: "../envs/python-minimal.yaml"
    resources:
    resources:
        mem_mb=12000,
        runtime=30,
        cpus=8
    script: "../scripts/processing/agg_novelty_zscores.py"
