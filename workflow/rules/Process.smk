from workflow_helpers import get_venues

rule agg_letters:
    input:
        letters=expand(rules.dl_gbq_letter_metadata.output, venue = get_venues(config)),
        impacts=expand(rules.dl_gbq_paper_impact.output, venue = get_venues(config)),
    output: AGG_LETTERS
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/agg_letters.py"

rule agg_nonletters:
    input: 
        letters=expand(rules.dl_gbq_letter_metadata.output, venue = get_venues(config)),
        impacts=expand(rules.dl_gbq_paper_impact.output, venue = get_venues(config)),
    output: AGG_NONLETTERS
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/agg_nonletters.py"
        
rule agg_fields:
    input:
        fields = expand(
            rules.dl_gbq_paper_fields.output, 
            venue = get_venues(config)
        ),
    params: FIELD_HIERARCHY
    output: AGG_FIELDS
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/agg_fields.py"

rule agg_dual_cite_trajectories:
    input: 
        traj = expand(
            rules.dl_gbq_dual_citation_trajectories.output,
            venue = get_venues(config)
        )
    output: AGG_DUAL_CITE_TRAJECTORIES
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/gather.py"

rule agg_paper_titles:
    input:
        titles = expand(
            rules.dl_gbq_citing_paper_titles.output,
            venue = get_venues(config)
        )
    output: AGG_PAPER_TITLES
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/gather.py"


rule calculate_paper_novelty:
    input: 
        zscores = ancient(rules.agg_novelty_zscores.output),
        letters = rules.agg_letters.output,
        nonletters = rules.agg_nonletters.output,
        refs = rules.dl_gbq_target_journal_refs.output
    output: PAPER_NOVELTY_SCORES
    conda: "../envs/python-minimal.yaml"
    params: 
        year_min = config["novelty"]["start_year"],
        year_max = config["novelty"]["end_year"] - 1
    script: "../scripts/processing/calculate_paper_novelty.py"

rule title_embedding_similarities:
    input: rules.dl_gbq_citing_paper_titles.output
    output: TITLE_EMBEDDING_SIMILARITIES
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/embed_titles_and_calculate_similarity.py"

rule agg_paper_title_similarities:
    input: 
        emb = expand(
            rules.title_embedding_similarities.output,
            venue = get_venues(config)
        )
    output: AGG_TITLE_EMBEDDING_SIMILARITIES
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/gather.py"

rule agg_scite_metrics:
    input: 
        metrics = expand(
            rules.dl_gbq_scite_metrics.output,
            venue = get_venues(config)
        )
    output: AGG_SCITE_METRICS
    conda: "../envs/python-minimal.yaml"
    script: "../scripts/processing/gather.py"