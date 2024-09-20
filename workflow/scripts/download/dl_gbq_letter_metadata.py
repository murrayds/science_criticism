"""
This script downloads data from Google BigQuery related to the impact of papers in a specific venue.
It executes a query that calculates the impact of papers over different time periods and saves the results to a local file.
"""
import pandas as pd

from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
temp = snakemake.config["bigquery"]["temp_path"]

venue = snakemake.config["venues"][snakemake.wildcards.venue]

letter_ids = pd.read_csv(snakemake.input[0])

client = bigquery.Client()

SOURCE_DATA_TABLE_REF = f"{temp}.temp_{venue}_letter_ids"
# Upload data from the file at the path given by the variable "letter_ids" to a temporary table on Google Big Query
job_config = bigquery.LoadJobConfig(
    source_format=bigquery.SourceFormat.CSV,
    skip_leading_rows=1,
    autodetect=True,
    schema=[
        bigquery.SchemaField("original_doi", "STRING"),
        bigquery.SchemaField("letter_doi", "STRING")
    ],
    write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE
)

job = client.load_table_from_dataframe(
    letter_ids,
    SOURCE_DATA_TABLE_REF,
    job_config=job_config
)

job.result()  # Wait for the job to complete

print(f"Data uploaded to temporary table: {SOURCE_DATA_TABLE_REF}")



# Now, we will proceed to gather the 
QUERY = f"""
SELECT
  orig.PaperId as original_id,
  lett.PaperId as letter_id,
  orig.Year as original_year,
  lett.Year as letter_year,
  IF(orig.DocSubTypes = "", 0, 1) AS retracted,
  "{snakemake.wildcards.venue}" as venue
FROM `{SOURCE_DATA_TABLE_REF}` s 
LEFT JOIN `{mag}.Papers` orig on orig.Doi = s.original_doi
LEFT JOIN `{mag}.Papers` lett on lett.Doi = s.letter_doi
"""

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
client.delete_table(SOURCE_DATA_TABLE_REF, not_found_ok=True)
print(f"Table '{SOURCE_DATA_TABLE_REF}' deleted.")
client.delete_table(TEMP_TABLE_REF, not_found_ok=True)
print(f"Table '{TEMP_TABLE_REF}' deleted.")
