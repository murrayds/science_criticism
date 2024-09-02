rule match_papers_for_metric_density_comparison:
    input: 
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: MATCHED_PAPERS_FOR_METRIC_DENSITY_COMPARISON
    script: "../scripts/processing/match_papers_for_metric_density_comparison.R"

rule match_papers_for_impact_comparison:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: MATCHED_PAPERS_FOR_IMPACT_COMPARISON
    script: "../scripts/processing/match_papers_for_impact_comparison.R"

rule match_authors:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.dl_gbq_career_histories.output
    output: MATCHED_AUTHORS
    script: "../scripts/processing/match_authors.R"