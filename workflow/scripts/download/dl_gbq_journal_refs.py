from google.cloud import bigquery
from google.api_core.exceptions import NotFound

from dl_helpers import extract_data_to_local_file, gen_random_sequence

# Retreive the year wildcard passed into the file
YEAR = int(snakemake.wildcards.year)

mag = snakemake.config["bigquery"]["mag_path"]
aps = snakemake.config["bigquery"]["aps_path"]
temp = snakemake.config["bigquery"]["temp_path"]

# Check if the table exists
BUILD_TABLE_QUERY = f"""
WITH target_journals AS (
  SELECT 
    JournalID
  FROM `{mag}.Journals`
  WHERE PaperCount > 1000
),
papers AS (
  SELECT
    p.PaperId as CitingPaperId,
    p.Year as CitingPaperYear
  FROM `{mag}.Papers` p
  INNER JOIN target_journals j on j.JournalId = p.JournalId
  WHERE DocType = "Journal"
  AND Doi IS NOT NULL
  AND ReferenceCount >= 5
  AND DocSubTypes = ""
),
aps_mag_citations AS (
  -- MAG has poor coverage of APS citations...lets supplement with info from APS
  SELECT 
    citing.PaperId as PaperId,
    cited.PaperId as PaperReferenceId
  FROM `{aps}.citations` aps
  LEFT JOIN `{mag}.Papers` citing on citing.Doi = UPPER(aps.citing_doi)
  LEFT JOIN `{mag}.Papers` cited on cited.Doi = UPPER(aps.cited_doi)
),
all_refs AS (
  SELECT DISTINCT *
  FROM (
    SELECT
      *
    FROM `{mag}.PaperReferences`
    UNION ALL 
    SELECT 
      *
    FROM aps_mag_citations
  )
),
refs AS (
  select
    p.*,
    r.PaperReferenceId,
    p2.JournalId as ReferenceJournalId,
    p2.Year as ReferenceYear
  FROM all_refs r
  INNER JOIN papers as p ON p.CitingPaperId = r.PaperId
  LEFT JOIN `{mag}.Papers` p2 on p2.PaperId = r.PaperReferenceId
  INNER JOIN target_journals j on j.JournalId = p2.JournalId
  WHERE p2.Doi IS NOT NULL
  AND p2.DocType = "Journal"
)
SELECT 
  CitingPaperId,
  CitingPaperYear,
  PaperReferenceId,
  ReferenceJournalId,
  ReferenceYear
FROM refs
"""

client = bigquery.Client()
# First, check to see if the table exists...if not create it
# Additionally, check if the table is more than 1 hour old, and if it is, re-create the table
MAIN_REF_TABLE = f"{temp}.temp_refs_table"

try:
    table = client.get_table(MAIN_REF_TABLE)
    print(f"Table {MAIN_REF_TABLE} already exists.")
    # Check if the table is more than 1 hour old
    from datetime import datetime, timedelta
    if datetime.now() - table.created > timedelta(hours=1):
        print(f"Table {MAIN_REF_TABLE} is more than 1 hour old. Re-creating it now.")
        # Set up the query job configuration
        job_config = bigquery.QueryJobConfig(
            destination=MAIN_REF_TABLE,
            write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
            use_query_cache=False
        )

        # Run the query and save results to the new table
        query_job = client.query(BUILD_TABLE_QUERY, job_config=job_config)
        query_job.result()  # Wait for the query to complete

        print(f"Query results re-exported to: {MAIN_REF_TABLE}")
    else:
        print(f"Table {MAIN_REF_TABLE} is less than 1 hour old. No re-creation needed.")
except NotFound:
    print(f"Table {MAIN_REF_TABLE} is not found. Creating it now.")
    
    # Set up the query job configuration
    job_config = bigquery.QueryJobConfig(
        destination=MAIN_REF_TABLE,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
        use_query_cache=False
    )

    # Run the query and save results to the new table
    query_job = client.query(BUILD_TABLE_QUERY, job_config=job_config)
    query_job.result()  # Wait for the query to complete

    print(f"Query results exported to: {MAIN_REF_TABLE}")


# To make downloading and downstream processing easier, we will query only 
# one year at a time...
SELECT_YEAR_QUERY = f"""
SELECT 
    CitingPaperId,
    ReferenceJournalId,
    ReferenceYear
FROM {MAIN_REF_TABLE}
WHERE CitingPaperYear = {YEAR}
"""

# Execute the query
query_job = client.query(SELECT_YEAR_QUERY)

random_seq = gen_random_sequence()
SELECTED_TABLE_REF = f"{temp}.temp_{random_seq}"

# Set up the query job configuration
job_config = bigquery.QueryJobConfig(
    destination=SELECTED_TABLE_REF,
    write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
    use_query_cache=False
)

# Run the query and save results to the new table
query_job = client.query(SELECT_YEAR_QUERY, job_config=job_config)
query_job.result()  # Wait for the query to complete


result = extract_data_to_local_file(
    table = SELECTED_TABLE_REF,
    local_filename = snakemake.output[0],
    client = client
) 


# Delete the temporary GBQ tables...
client.delete_table(SELECTED_TABLE_REF, not_found_ok=True)
print(f"Table '{SELECTED_TABLE_REF}' deleted.")


