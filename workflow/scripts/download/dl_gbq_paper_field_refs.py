# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file 

QUERY = """
WITH papers AS (
  SELECT 
    PaperId
  FROM ccnr-success.mag.Papers
  WHERE JournalId IN (
    24807848,
    3880285,
    137773608,
    125754415
  )
  AND year >= 2000
  AND DocType = "Journal"
),
refs AS (
  SELECT
    r.*
  FROM ccnr-success.mag.PaperReferences r
  INNER JOIN papers as p on p.PaperId = r.PaperId
),
ref_fields AS (
  SELECT 
    r.*,
    f.FieldOfStudyId,
    fos.FieldOfStudyId as ParentField
  FROM ccnr-success.mag.PaperFieldsOfStudy f
  INNER JOIN refs r on r.PaperReferenceId = f.PaperId
  LEFT JOIN ccnr-success.mag.FieldOfStudyChildren fos on fos.ChildFieldOfStudyId = f.FieldOfStudyId
  WHERE f.Score > 0
)
SELECT 
  base.PaperId,
  base.field
FROM (
  SELECT 
    r.PaperId,
    r.PaperReferenceId,
    COALESCE(ParentField, FieldOfStudyId) AS field
  FROM ref_fields r
) base
LEFT JOIN ccnr-success.mag.FieldsOfStudy fos on fos.FieldOfStudyId = field
WHERE Level = 0
ORDER BY PaperId, PaperReferenceId
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
