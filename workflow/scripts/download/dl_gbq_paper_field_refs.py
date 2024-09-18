# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
aps = snakemake.config["bigquery"]["aps_path"]
temp = snakemake.config["bigquery"]["temp_path"]

QUERY = f"""
WITH papers AS (
  SELECT 
    PaperId
  FROM {mag}.Papers
  WHERE JournalId IN (
    {",".join(map(str, snakemake.config["venues"].values()))}
  )
  AND year >= 2000
  AND DocType = "Journal"
  AND DocSubTypes = ""
),
aps_mag_citations AS (
  -- MAG has poor coverage of APS citations...lets supplement with info from APS
  SELECT 
    citing.PaperId as PaperId,
    cited.PaperId as PaperReferenceId
  FROM `{aps}.citations` aps
  LEFT JOIN `{mag}.Papers` citing on citing.Doi = aps.citing_doi
  LEFT JOIN `{mag}.Papers` cited on cited.Doi = aps.cited_doi
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
  SELECT
    r.*
  FROM all_refs r
  INNER JOIN papers as p on p.PaperId = r.PaperId
),
ref_fields AS (
  SELECT 
    r.*,
    f.FieldOfStudyId,
    fos.FieldOfStudyId as ParentField
  FROM {mag}.PaperFieldsOfStudy f
  INNER JOIN refs r on r.PaperReferenceId = f.PaperId
  LEFT JOIN {mag}.FieldOfStudyChildren fos on fos.ChildFieldOfStudyId = f.FieldOfStudyId
  WHERE f.Score > 0
)
SELECT 
  base.PaperId,
  base.PaperReferenceId,
  base.field
FROM (
  SELECT 
    r.PaperId,
    r.PaperReferenceId,
    COALESCE(ParentField, FieldOfStudyId) AS field
  FROM ref_fields r
) base
LEFT JOIN {mag}.FieldsOfStudy fos on fos.FieldOfStudyId = field
WHERE Level = 0
ORDER BY PaperId, PaperReferenceId
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
