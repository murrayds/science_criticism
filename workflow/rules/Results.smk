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

rule table_counts:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output
    output: COUNTS_TABLE 
    conda: "../envs/r-conda.yaml"
    script: "../scripts/tables/table_counts.R"

rule table_fit_by_venue:
    input: 
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.agg_fields.output
    output: FIT_BY_VENUE_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/tables/table_fit_by_venue.R"

rule plot_match_diagnostic_impact:
    input:
        expand(
            rules.match_papers_for_impact_comparison.output,
            delay = 3,
            cite_tolerance = 0.05,
            year_tolerance = 2
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

rule plotdata_paper_comparison:
    input: 
        expand(
            rules.match_papers_for_impact_comparison.output,
            delay = 3,
            cite_tolerance = 0.05,
            year_tolerance = 1
        )
    output:
        PAPER_IMPACT_COMPARISON_PLOTDATA
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plotdata_paper_impact_comparison.R"

rule plotdata_author_comparison:
    input: lambda wc: 
        expand(
            rules.match_authors.output,
            cite_tolerance = 0.10,
            prod_tolerance = 0.10,
            authorship = wc.authorship,
            metric = wc.metric
        )
    output: AUTHOR_COMPARISON_PLOTDATA
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plotdata_matched_author_comparison.R"

rule plot_matched_comparison:
    input:
        rules.plotdata_paper_comparison.output,
        expand(
            rules.plotdata_author_comparison.output,
            authorship = config["matching"]["authorship"],
            metric = ["citations", "prod"]
        )
    output: MATCHED_COMPARISON_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_matched_comparison.R"

rule plot_paper_lagtype_comparison:
    input:
        rules.plotdata_paper_comparison.output,
    output: PAPER_COMPARISON_LAGTYPE_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_comparison_lag.R"

rule impact_temporal_comparison:
    input:
        letters = rules.agg_letters.output,
        traj = rules.agg_dual_cite_trajectories.output,
        titles = rules.agg_paper_titles.output
    output: IMPACT_TEMPORAL_COMPARISON
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_impact_temporal_comparison.R"

rule paper_author_demographic_ame:
    input: rules.author_marginal_probabilities.output
    output: PAPER_AUTHOR_DEMOGRAPHIC_AME
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_author_demographics_ame.R"

rule paper_tweet_distribution:
    input:
        rules.agg_letters.output,
        rules.agg_nonletters.output,
        rules.dl_gbq_sciscinet_features.output
    output: PAPER_TWEET_DISTRIBUTION
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_tweet_distribution.R"

rule paper_news_tweet_distribution:
    input:
        expand(
            rules.match_papers_for_metric_density_comparison.output,
            cite_tolerance = 0.05,
            year_tolerance = 3
        ),
        rules.dl_gbq_sciscinet_features.output
    output: PAPER_NEWS_TWEET_PROPORTION
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_paper_news_tweet_proportion.R"
    
rule embedding_cocite_density_plot:
    input:
        letters = rules.agg_letters.output,
        titles = rules.agg_paper_titles.output,
        traj = rules.agg_dual_cite_trajectories.output,
        emb = rules.agg_paper_title_similarities.output
    output:
        plot = EMBEDDING_COCITE_DENSITY_PLOT,
        table = EMBEDDING_COCITE_DENSITY_PLOT_TESTS_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_embedding_cocite_density.R"

rule plot_letter_altmetric_comparison:
    input:
        rules.agg_letters.output,
        rules.dl_gbq_sciscinet_features.output
    output: LETTER_ALTMETRIC_COMPARISON_PLOT
    conda: "../envs/r-conda.yaml"
    script: "../scripts/plotting/plot_letter_altmetric_comparison.R"

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

rule table_author_match_quality:
    input: lambda wc: 
        expand(
            rules.match_authors.output,
            cite_tolerance = 0.10,
            prod_tolerance = 0.10,
            authorship = wc.authorship,
        )
    output: AUTHOR_MATCH_QUALITY_TABLE
    conda: "../envs/r-conda.yaml"
    script: "../scripts/tables/table_author_match_quality.R"