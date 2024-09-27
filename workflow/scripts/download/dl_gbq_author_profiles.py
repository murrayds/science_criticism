# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

mag = snakemake.config["bigquery"]["mag_path"]
aps = snakemake.config["bigquery"]["aps_path"]
temp = snakemake.config["bigquery"]["temp_path"]

QUERY = f"""
WITH target_authors AS (
  SELECT DISTINCT 
    AuthorId, AffiliationId, PaperId
  FROM (
    SELECT 
      paa.AuthorId,
      paa.AffiliationId,
      paa.PaperId
    FROM `{mag}.PaperAuthorAffiliations` paa
    LEFT JOIN `{mag}.Papers` p on p.PaperId = paa.PaperId
    WHERE p.JournalId in (
      {",".join(map(str, snakemake.config["venues"].values()))}
    )
    AND p.Year Between 1998 and 2020
    AND p.DocType = "Journal"
    AND DocSubTypes = ""
  )
),
author_order AS (
  SELECT DISTINCT
    paa.PaperId,
    paa.AuthorId,
    CASE paa.AuthorSequenceNumber
      WHEN 1 THEN 'f'
      WHEN seq.num_authors then 'l'
      ELSE 'm'
      END AS author_position,
    seq.num_authors
  FROM `{mag}.PaperAuthorAffiliations` paa 
  LEFT JOIN (
    SELECT 
      paa.PaperId,
      MAX(paa.AuthorSequenceNumber) as num_authors
    FROM `{mag}.PaperAuthorAffiliations` paa
    GROUP BY paa.PaperId
  ) seq on seq.PaperId = paa.PaperId
)

SELECT DISTINCT
  t.*, 
  ord.author_position,
  ord.num_authors,
  aff.GridId,
  aff.NormalizedName as affiliation_name,
  aff.Iso3166Code AS affiliation_country,
  gender.Gender
FROM target_authors t
LEFT JOIN {mag}.Affiliations aff ON aff.AffiliationId = t.AffiliationId
LEFT JOIN author_order ord on ord.PaperId = t.PaperId and ord.AuthorId = t.AuthorId
LEFT JOIN {mag}.Authors auth ON auth.AuthorId = t.AuthorId
LEFT JOIN {mag}.GenderNames gender on gender.NormalizedName = auth.NormalizedName
WHERE ord.author_position IN ("f", "l")
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
