from os.path import join as j

# Use a config file
configfile: "config.yaml"

RAW_DIR = j(config["data_dir"], "raw")
DERIVED_DIR = j(config["data_dir"], "derived")

JOURNAL_REFS = j(RAW_DIR, "dl", "journal_refs_{year}.csv")

AGG_OBSERVED_JOURNAL_COOCCURENCE_COUNTS = j(DERIVED_DIR, "novelty", "observed", "agg_observed_journal_cooccurence_counts.pickle")
AGG_NULL_JOURNAL_COOCCURENCE_COUNTS = j(DERIVED_DIR, "novelty", "null", "agg_null_journal_cooccurence_counts_iter{iter}.pickle")
NOVELTY_ZSCORES = j(DERIVED_DIR, "novelty", "journal_zscores_{year}.pickle")
AGG_NOVELTY_ZSCORES = j(DERIVED_DIR, "novelty", "agg_journal_zscores.pickle")
PAPER_NOVELTY_SCORES = j(DERIVED_DIR, "novelty", "paper_zscores.csv")

YEARS_OF_STUDY = list(range(2000, 2021))
NULL_ITERS = list(range(0,10))

rule all:
    input:
        expand(JOURNAL_REFS, year=YEARS_OF_STUDY),
        expand(OBSERVED_JOURNAL_COOCCURENCE_COUNTS, year=YEARS_OF_STUDY),
        expand(NULL_JOURNAL_COOCCURENCE_COUNTS, year=YEARS_OF_STUDY, iter=NULL_ITERS),
        expand(NOVELTY_ZSCORES, year = YEARS_OF_STUDY),
        AGG_NOVELTY_ZSCORES,
        PAPER_NOVELTY_SCORES

rule dl_gbq_journal_refs:
    output: JOURNAL_REFS
    script: "../scripts/dl/dl_gbq_journal_refs.py"

rule count_journal_cooccurence_observed:
    input: rules.dl_gbq_journal_refs.output
    output: OBSERVED_JOURNAL_COOCCURENCE_COUNTS
    params: shuffle = False
    script: "../scripts/processing/count_journal_cooccurence.py"

rule count_journal_cooccurence_null:
    input: rules.dl_gbq_journal_refs.output
    output: NULL_JOURNAL_COOCCURENCE_COUNTS
    params: shuffle = True
    script: "../scripts/processing/count_journal_cooccurence.py"


rule count_journal_cooccurence_observed:
    input: rules.dl_gbq_journal_refs.output
    output: OBSERVED_JOURNAL_COOCCURENCE_COUNTS
    params: shuffle = False
    script: "../scripts/processing/count_journal_cooccurence.py"

rule count_journal_cooccurence_null:
    input: rules.dl_gbq_journal_refs.output
    output: NULL_JOURNAL_COOCCURENCE_COUNTS
    params: shuffle = True
    script: "../scripts/processing/count_journal_cooccurence.py"

rule calculate_novelty_zscores:
    input: 
        observed = rules.count_journal_cooccurence_observed.output,
        null = lambda wc: collect(
                    rules.count_journal_cooccurence_null.output,
                    iter = NULL_ITERS,
                    year = wc.year
                )
    output: NOVELTY_ZSCORES
    script: "../scripts/processing/calculate_novelty_zscores.py"

rule agg_novelty_zscores:
    input:
        ancient(
            expand(
                rules.calculate_novelty_zscores.output, 
                year = YEARS_OF_STUDY
            )
        )
    output: AGG_NOVELTY_ZSCORES
    script: "../scripts/processing/agg_novelty_zscores.py"