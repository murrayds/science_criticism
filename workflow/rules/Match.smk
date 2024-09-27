rule match_papers_for_metric_density_comparison:
    input: 
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: MATCHED_PAPERS_FOR_METRIC_DENSITY_COMPARISON
    conda: "../envs/r-conda.yaml"
    script: "../scripts/matching/match_papers_for_metric_density_comparison.R"

rule _match_papers_for_impact_comparison_split:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: temp(_MATCHED_PAPERS_FOR_IMPACT_COMPARISON_SPLIT)
    resources:
        mem_mb = 8000,
        cpus = 8
    conda: "../envs/r-conda.yaml"
    script: "../scripts/matching/match_papers_for_impact_comparison.R"

rule match_papers_for_impact_comparison:
    input: lambda wc:
        expand(
            rules._match_papers_for_impact_comparison_split.output,
            venue = get_venues(config),
            delay = wc.delay,
            cite_tolerance = wc.cite_tolerance,
            year_tolerance = wc.year_tolerance
    )
    output: MATCHED_PAPERS_FOR_IMPACT_COMPARISON
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/gather.py"

# Process each venue individually for better paralellization
rule _match_authors_split:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.dl_gbq_career_histories.output
    output:
        temp(_MATCHED_AUTHORS_SPLIT)
    resources:
        mem_mb = 8000,
        cpus = 8
    conda: "../envs/r-conda.yaml"
    script: "../scripts/matching/match_authors.R"

rule match_authors:
    input: lambda wc:
        expand(
            rules._match_authors_split.output,
            venue = get_venues(config),
            authorship = wc.authorship,
            cite_tolerance = wc.cite_tolerance,
            prod_tolerance = wc.prod_tolerance
    )
    output: MATCHED_AUTHORS
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/gather.py"

rule agg_matched_paper_impact_diagnostics:
    input: 
        rules.agg_letters.output,
        collect(
            rules.match_papers_for_impact_comparison.output,
            delay = config["matching"]["impact_delay"],
            cite_tolerance = config["matching"]["cite_tolerance"],
            year_tolerance = config["matching"]["year_tolerance"]
        )
    output: AGG_PAPER_IMPACT_MATCHED_DIAGNOSTICS
    conda: "../envs/r-conda.yaml"
    script: "../scripts/matching/agg_paper_impact_match_diagnostics.R"

rule agg_matched_author_diagnostics:
    input: 
        rules.agg_letters.output,
        collect(
            rules.match_authors.output,
            authorship = config["matching"]["authorship"],
            cite_tolerance = config["matching"]["cite_tolerance"],
            prod_tolerance = config["matching"]["productivity_tolerance"],
        )
    output: AGG_AUTHOR_MATCHED_DIAGNOSTICS
    conda: "../envs/r-conda.yaml"
    script: "../scripts/matching/agg_author_match_diagnostics.R"

rule table_paper_impact_match_diagnostics:
    input: rules.agg_matched_paper_impact_diagnostics.output
    output: PAPER_IMPACT_MATCH_DIAGNOSTICS_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/tables/table_paper_match_diagnostics.R"

rule table_author_match_diagonstics:
    input: rules.agg_matched_author_diagnostics.output 
    output: AUTHOR_MATCH_DIAGNOSTICS_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/tables/table_author_match_diagnostics.R"