# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
aps = snakemake.config["bigquery"]["aps_path"]
temp = snakemake.config["bigquery"]["temp_path"]

QUERY = f"""
WITH papers AS (
  SELECT 
    PaperId,
    Year
  FROM {mag}.Papers
  WHERE JournalId IN (
    {",".join(map(str, snakemake.config["venues"].values()))}
  )
  AND Year >= 2000
  AND DocType = "Journal"
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
cites AS (
  SELECT
    r.PaperReferenceId as CitedPaperId,
    r.PaperId as CitingPaperId
  FROM all_refs r
  INNER JOIN papers as p1 on p1.PaperId = r.PaperReferenceId
  LEFT JOIN {mag}.Papers as p2 on p2.PaperId = r.PaperId
  WHERE (p2.Year - p1.Year) <= 5
),
ref_fields AS (
  SELECT 
    r.*,
    f.FieldOfStudyId,
    fos.FieldOfStudyId as ParentField
  FROM {mag}.PaperFieldsOfStudy f
  INNER JOIN cites r on r.CitingPaperId = f.PaperId
  LEFT JOIN {mag}.FieldOfStudyChildren fos on fos.ChildFieldOfStudyId = f.FieldOfStudyId
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
LEFT JOIN {mag}.FieldsOfStudy fos on fos.FieldOfStudyId = field
WHERE Level = 0
ORDER BY CitedPaperId, CitingPaperId
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
