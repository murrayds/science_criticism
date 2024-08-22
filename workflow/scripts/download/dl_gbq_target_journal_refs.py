from google.cloud import bigquery
from google.api_core.exceptions import NotFound

from dl_helpers import extract_data_to_local_file 

QUERY = """
WITH target_journals AS (
  SELECT 
    JournalID
  FROM `ccnr-success.mag.Journals`
  WHERE PaperCount > 1000
),
papers AS (
  SELECT
    p.PaperId as CitingPaperId,
    p.Year as CitingPaperYear
  FROM `ccnr-success.mag.Papers` p
  WHERE JournalId IN (
    24807848,
    3880285,
    137773608,
    125754415
  ) 
  AND DocType = "Journal"
  AND year >= 2000
),
refs AS (
  select
    p.*,
    r.PaperReferenceId,
    p2.JournalId as ReferenceJournalId,
    p2.Year as ReferenceYear
  FROM `ccnr-success.mag.PaperReferences` as r
  INNER JOIN papers as p ON p.CitingPaperId = r.PaperId
  LEFT JOIN `ccnr-success.mag.Papers` p2 on p2.PaperId = r.PaperReferenceId
  INNER JOIN target_journals j on j.JournalId = p2.JournalId
  WHERE p2.Doi IS NOT NULL
  AND p2.DocType = "Journal"
)
SELECT 
  CitingPaperId,
  ReferenceJournalId,
FROM refs
"""

client = bigquery.Client()

# Execute the query
TEMP_TABLE_REF = "ccnr-success.dmurray.temp_journal_paper_refs"

# Set up the query job configuration
job_config = bigquery.QueryJobConfig(
    destination=TEMP_TABLE_REF,
    write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
    use_query_cache=False
)

# Run the query and save results to the new table
query_job = client.query(QUERY, job_config=job_config)
query_job.result()  # Wait for the query to complete

result = extract_data_to_local_file(
    table = TEMP_TABLE_REF,
    local_filename = snakemake.output[0],
    client = client
) 

# Delete the temporary GBQ tables...
client = bigquery.Client()
client.delete_table(TEMP_TABLE_REF, not_found_ok=True)
print(f"Table '{TEMP_TABLE_REF}' deleted.")



