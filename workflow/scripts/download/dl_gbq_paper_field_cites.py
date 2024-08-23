# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file 

QUERY = """
WITH papers AS (
  SELECT 
    PaperId,
    Year
  FROM ccnr-success.mag.Papers
  WHERE JournalId IN (
    24807848,
    3880285,
    137773608,
    125754415
  )
  AND Year >= 2000
  AND DocType = "Journal"
),
cites AS (
  SELECT
    r.PaperReferenceId as CitedPaperId,
    r.PaperId as CitingPaperId
  FROM ccnr-success.mag.PaperReferences r
  INNER JOIN papers as p1 on p1.PaperId = r.PaperReferenceId
  LEFT JOIN ccnr-success.mag.Papers as p2 on p2.PaperId = r.PaperId
  WHERE (p2.Year - p1.Year) <= 5
),
ref_fields AS (
  SELECT 
    r.*,
    f.FieldOfStudyId,
    fos.FieldOfStudyId as ParentField
  FROM ccnr-success.mag.PaperFieldsOfStudy f
  INNER JOIN cites r on r.CitingPaperId = f.PaperId
  LEFT JOIN ccnr-success.mag.FieldOfStudyChildren fos on fos.ChildFieldOfStudyId = f.FieldOfStudyId
  WHERE f.Score > 0
)
SELECT DISTINCT
  base.CitedPaperId,
  base.CitingPaperId,
  base.field
FROM (
  SELECT 
    r.CitedPaperId,
    r.CitingPaperId,
    COALESCE(ParentField, FieldOfStudyId) AS field
  FROM ref_fields r
) base
LEFT JOIN ccnr-success.mag.FieldsOfStudy fos on fos.FieldOfStudyId = field
WHERE Level = 0
ORDER BY CitedPaperId, CitingPaperId
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
