rule author_marginal_probabilities:
    input:
        authors = rules.dl_gbq_author_profiles.output,
        histories = rules.dl_gbq_career_histories.output,
        matched = expand(
            rules.match_papers_for_metric_density_comparison.output,
            cite_tolerance = 0.05,
            year_tolerance = 1
        )
    params: ELITE_UNIVERSITIES
    output: AUTHOR_MARGINAL_PROBABILITIES
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/gen_author_marginal_probabilities.R"

rule plot_impact_rank_histogram:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: IMPACT_RANK_HISTOGRAM
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_impact_rank_histogram.R"

rule plot_impact_likelihood_scatter:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output
    output: IMPACT_LIKELIHOOD_SCATTER
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_impact_likelihood_scatter.R"

rule plot_paper_metrics_density:
    input: 
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output,
        rules.dl_gbq_paper_field_refs.output,
        rules.dl_gbq_paper_field_cites.output,
        rules.calculate_paper_novelty.output, # TODO: Replace with novelty once complete...
        expand(
            rules.match_papers_for_metric_density_comparison.output,
            cite_tolerance = 0.05,
            year_tolerance = 3
        )
    output: 
        PAPER_METRICS_DENSITY_PLOT,
        PAPER_METRICS_DENSITY_PLOT_TEST_1SAMPLE_TABLE,
        PAPER_METRICS_DENSITY_PLOT_TEST_2SAMPLE_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_metrics_density.R"


rule table_fit_by_venue:
    input: 
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: FIT_BY_VENUE_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/tables/table_fit_by_venue.R"


rule plot_impact_comparison:
    input: 
        expand(
            rules.match_papers_for_impact_comparison.output,
            delay = 3,
            cite_tolerance = 0.10,
            year_tolerance = 1
        )
    output:
        POOLED_IMPACT_COMPARISON_PLOT,
        PAIRWISE_IMPACT_COMPARISON_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_impact_comparison.R"

rule plot_match_diagnostic_impact:
    input:
        expand(
            rules.match_papers_for_impact_comparison.output,
            delay = 3,
            cite_tolerance = 0.10,
            year_tolerance = 1
        )
    output:
        MATCHING_DIAGNOSTIC_IMPACT_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_match_diagnostic_impact.R"

rule plot_field_representation:
    input: 
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: FIELD_REPRESENTATION_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_field_representation.R"


rule plot_author_comparison:
    input: lambda wc: 
        expand(
            rules.match_authors.output,
            cite_tolerance = 0.15,
            prod_tolerance = 1.0,
            authorship = wc.authorship,
            metric = wc.metric
        )
    output: PAIRWISE_AUTHOR_COMPARISON_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_matched_author_comparison.R"

rule plot_cite_ratio:
    input:
        letters = rules.agg_letters.output,
        traj = rules.agg_dual_cite_trajectories.output 
    output: CITE_RATIO_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_cite_ratio.R"

rule paper_author_demographic_ame:
    input: rules.author_marginal_probabilities.output
    output: PAPER_AUTHOR_DEMOGRAPHIC_AME
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_author_demographics_ame.R"
