"""
This script downloads data from Google BigQuery related to the impact of papers in a specific venue.
It executes a query that calculates the impact of papers over different time periods and saves the results to a local file.
"""
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
aps = snakemake.config["bigquery"]["aps_path"]
temp = snakemake.config["bigquery"]["temp_path"]

venue = snakemake.config["venues"][snakemake.wildcards.venue]

QUERY = f"""
WITH aps_mag_citations AS (
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
reference_count AS (
  SELECT
    PaperId,
    COALESCE(ReferenceCount, 0) as ReferenceCount
  FROM (
    SELECT 
      r.PaperId,
      COUNT(r.PaperReferenceId) as ReferenceCount
    FROM all_refs r
    LEFT JOIN {mag}.Papers p on r.PaperId = p.PaperId
    GROUP BY r.PaperId
  )
)
SELECT
  cited.PaperId as id,
  ANY_VALUE(cited.Year) as year,
  ANY_VALUE(EXTRACT(MONTH FROM cited.Date)) AS month,
  ANY_VALUE(rc.ReferenceCount) as ReferenceCount,
  "{snakemake.wildcards.venue}" as venue,
  SUM(IF((citing.Year - cited.Year) <= 1, 1, 0)) AS impact_1year,
  SUM(IF((citing.Year - cited.Year) <= 2, 1, 0)) AS impact_2year,
  SUM(IF((citing.Year - cited.Year) <= 3, 1, 0)) AS impact_3year,
  SUM(IF((citing.Year - cited.Year) <= 4, 1, 0)) AS impact_4year,
  SUM(IF((citing.Year - cited.Year) <= 5, 1, 0)) AS impact_5year,
  SUM(IF((citing.Year - cited.Year) <= 6, 1, 0)) AS impact_6year,
  SUM(IF((citing.Year - cited.Year) <= 7, 1, 0)) AS impact_7year,
  SUM(IF((citing.Year - cited.Year) <= 8, 1, 0)) AS impact_8year,
  SUM(IF((citing.Year - cited.Year) <= 9, 1, 0)) AS impact_9year,
  SUM(IF((citing.Year - cited.Year) <= 10, 1, 0)) AS impact_10year,
  SUM(IF((citing.Year - cited.Year) <= 11, 1, 0)) AS impact_11year,
  SUM(IF((citing.Year - cited.Year) <= 12, 1, 0)) AS impact_12year,
  SUM(IF((citing.Year - cited.Year) <= 13, 1, 0)) AS impact_13year,
  SUM(IF((citing.Year - cited.Year) <= 14, 1, 0)) AS impact_14year,
  SUM(IF((citing.Year - cited.Year) <= 15, 1, 0)) AS impact_15year,
FROM all_refs r 
LEFT JOIN `{mag}.Papers` citing on citing.PaperId = r.PaperId
LEFT JOIN `{mag}.Papers` cited on cited.PaperId = r.PaperReferenceId
LEFT JOIN reference_count rc on rc.PaperId = cited.PaperId
WHERE 
  cited.JournalId = {venue}
  AND citing.DocType = "Journal"
  AND cited.DocType = "Journal"
  AND cited.Year >= 2000
  AND cited.DocSubTypes = ""
GROUP BY cited.PaperId
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
    client = client,
    random_seq = random_seq
) 

# Delete the temporary GBQ tables...
client.delete_table(TEMP_TABLE_REF, not_found_ok=True)
print(f"Table '{TEMP_TABLE_REF}' deleted.")
