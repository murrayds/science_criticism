rule dl_gbq_letter_metadata:
    input: LETTER_IDS
    output: LETTER_METADATA,
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_letter_metadata.py"

rule dl_gbq_paper_impact:
    output: PAPER_IMPACT,
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_paper_impact.py"

rule dl_gbq_paper_fields:
    output: PAPER_FIELDS
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_paper_fields.py"

rule dl_gbq_journal_refs:
    output: JOURNAL_REFS
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_journal_refs.py"

rule dl_gbq_target_journal_refs:
    output: TARGET_JOURNAL_REFS
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_target_journal_refs.py"

rule dl_gbq_paper_field_refs:
    output: PAPER_FIELD_REFS
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_paper_field_refs.py"

rule dl_gbq_paper_field_cites:
    output: PAPER_FIELD_CITES
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_paper_field_cites.py"

rule dl_gbq_career_histories:
    output: CAREER_HISTORIES
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_career_histories.py"

rule dl_gbq_dual_citation_trajectories:
    input: LETTER_IDS
    output: DUAL_CITE_TRAJECTORIES,
    threads: workflow.cores * 0.5 if workflow.cores >= 2 else 1
    conda: "../envs/python-gbq.yaml"
    script: "../scripts/download/dl_gbq_dual_citation_trajectories.py"