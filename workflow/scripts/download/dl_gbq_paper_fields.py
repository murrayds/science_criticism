"""
This script downloads data from Google BigQuery related to the fields assigned to papers in a specific venue.
It executes a query that calculates the impact of papers over different time periods and saves the results to a local file.
"""
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
temp = snakemake.config["bigquery"]["temp_path"]

venue = snakemake.config["venues"][snakemake.wildcards.venue]

QUERY = f"""
WITH target_papers AS (
  SELECT 
    PaperId 
  FROM {mag}.Papers 
  WHERE 
    JournalId = {venue}
    AND DocType = "Journal"
    AND Year >= 2000
    AND DocSubTypes = ""
)

SELECT
  p.PaperId AS id,
  f.FieldOfStudyId AS field,
  fos.Level as level
FROM {mag}.PaperFieldsOfStudy f
INNER JOIN target_papers p on p.PaperId = f.PaperId
LEFT JOIN {mag}.FieldsOfStudy fos on fos.FieldOfStudyId = f.FieldOfStudyId
WHERE 
  (fos.Level = 0 OR fos.Level = 1)
  AND f.Score > 0
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
