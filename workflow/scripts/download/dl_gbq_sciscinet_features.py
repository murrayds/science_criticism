# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

sciscinet = snakemake.config["bigquery"]["sciscinet_path"]
temp = snakemake.config["bigquery"]["temp_path"]

QUERY = f"""
SELECT 
  PaperId as id,
  doi,
  Atyp_10pct_Z,
  Tweet_Count, 
  Newsfeed_Count,
  CAST(Team_Size AS INT64) as Team_Size
FROM `{sciscinet}.SciSciNetPapers`
WHERE JournalId IN (
  {",".join(map(str, snakemake.config["venues"].values()))}
)
AND year >= 2000
AND DocType = "Journal"
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



