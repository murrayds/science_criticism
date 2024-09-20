from google.cloud import bigquery
from google.api_core.exceptions import NotFound

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
temp = snakemake.config["bigquery"]["temp_path"]

QUERY = f"""
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
  WHERE JournalId IN (
    {",".join(map(str, snakemake.config["venues"].values()))}
  ) 
  AND DocType = "Journal"
  AND DocSubTypes = ""
  AND year >= 2000
),
refs AS (
  select
    p.*,
    r.PaperReferenceId,
    p2.JournalId as ReferenceJournalId,
    p2.Year as ReferenceYear
  FROM `{mag}.PaperReferences` as r
  INNER JOIN papers as p ON p.CitingPaperId = r.PaperId
  LEFT JOIN `{mag}.Papers` p2 on p2.PaperId = r.PaperReferenceId
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
random_seq = gen_random_sequence()
TEMP_TABLE_REF = f"{temp}.temp_{random_seq}"

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
client.delete_table(TEMP_TABLE_REF, not_found_ok=True)
print(f"Table '{TEMP_TABLE_REF}' deleted.")



