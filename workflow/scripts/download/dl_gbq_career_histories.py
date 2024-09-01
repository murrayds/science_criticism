# Check if the table exists
from google.cloud import bigquery

from dl_helpers import extract_data_to_local_file, gen_random_sequence

QUERY = f"""
WITH citations AS (
  SELECT 
    cited.PaperId,
    cited.Year,
    COUNT(*) as impact_3year
  FROM `ccnr-success.mag.PaperReferences` r 
  LEFT JOIN `ccnr-success.mag.Papers` citing on citing.PaperId = r.PaperId -- merge paper metadata
  LEFT JOIN `ccnr-success.mag.Papers` cited on cited.PaperId = r.PaperReferenceId -- merge reference metadata
  WHERE citing.Year BETWEEN cited.Year AND cited.Year + 3
    AND citing.DocType = "Journal" AND cited.DocType = "Journal"
    AND cited.Year >= 1980
  GROUP BY cited.PaperId, cited.Year
),
field_averages AS (
  SELECT
    pfos.FieldOfStudyId,
    cited.Year,
    AVG(impact_3year) as AvgFieldCitations
  FROM `ccnr-success.mag.PaperFieldsOfStudy` pfos 
  LEFT JOIN `ccnr-success.mag.FieldsOfStudy` finfo on finfo.FieldOfStudyId = pfos.FieldOfStudyId
  LEFT JOIN citations cited on cited.PaperId = pfos.PaperId
  WHERE pfos.Score > 0 AND finfo.Level = 1 -- only calculate at level 1 field...
  GROUP BY pfos.FieldOfStudyId, cited.Year
),
normalized_impact AS (
  SELECT 
    PaperId,
    AVG(impact_3year_norm) as impact_3year_norm 
    FROM (
      SELECT
        cited.PaperId,
        favg.FieldOfStudyId,
        cited.impact_3year / favg.AvgFieldCitations as impact_3year_norm
      FROM `ccnr-success.mag.PaperFieldsOfStudy` pfos 
      LEFT JOIN citations cited on cited.PaperId = pfos.PaperId
      INNER JOIN field_averages favg on favg.FieldOfStudyId = pfos.FieldOfStudyId and favg.Year = cited.Year
    )
    GROUP BY PaperId
),
level0_field AS (
  SELECT
    f.PaperId,
    f.field,
    fos.NormalizedName
  FROM (
    SELECT 
      pfos.PaperId,
      COALESCE(fosc.FieldOfStudyId, pfos.FieldOfStudyId) as field
    FROM `ccnr-success.mag.PaperFieldsOfStudy` pfos
    LEFT JOIN `ccnr-success.mag.FieldOfStudyChildren` fosc on fosc.ChildFieldOfStudyId = pfos.FieldOfStudyId
    WHERE pfos.Score > 0 
  ) f
  LEFT JOIN ccnr-success.mag.FieldsOfStudy fos on fos.FieldOfStudyId = f.field
  WHERE fos.Level = 0
),
level0_field_top AS (
  SELECT DISTINCT
    PaperId,
    FIRST_VALUE(field) OVER (PARTITION BY PaperId ORDER BY occurences DESC) field,
    FIRST_VALUE(NormalizedName) OVER (PARTITION BY PaperId ORDER BY occurences DESC) name
    FROM (
      SELECT 
        PaperId, 
        field, 
        NormalizedName,
        COUNT (*) AS occurences
      FROM level0_field
      GROUP BY PaperId, field, NormalizedName
    )
),
-- Now we move on to constructing the career histories...
-- We need to get all authors who published at least once in our target papers...
target_authors AS (
  SELECT DISTINCT 
    AuthorId
  FROM (
    SELECT 
      paa.AuthorId, 
    FROM `ccnr-success.mag.PaperAuthorAffiliations` paa
    LEFT JOIN `ccnr-success.mag.Papers` p on p.PaperId = paa.PaperId
    WHERE p.JournalId in (
      {",".join(map(str, snakemake.config["venues"].values()))}
    )
    AND p.Year Between 1998 and 2020
    AND p.DocType = "Journal"
    AND DocSubTypes = ""
  )
),
-- Now we can select all papers by the relevant authors...
-- One thing we should get, though, is career history...
career_histories AS (
  SELECT DISTINCT
    ta.AuthorId,
    p.PaperId,
    p.Year
  FROM `ccnr-success.mag.PaperAuthorAffiliations` paa
  INNER JOIN target_authors ta on ta.AuthorId = paa.AuthorId
  LEFT JOIN `ccnr-success.mag.Papers` p on p.PaperId = paa.PaperId
  WHERE p.DocType = "Journal" 
  AND p.Year Between 1980 and 2020
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
  FROM `ccnr-success.mag.PaperAuthorAffiliations` paa 
  LEFT JOIN (
    SELECT 
      paa.PaperId,
      MAX(paa.AuthorSequenceNumber) as num_authors
    FROM `ccnr-success.mag.PaperAuthorAffiliations` paa
    GROUP BY paa.PaperId
  ) seq on seq.PaperId = paa.PaperId
)

SELECT 
  hist.*, 
  ord.author_position,
  ord.num_authors,
  impact.impact_3year_norm,
  f.field
from career_histories hist
LEFT JOIN author_order ord on ord.PaperId = hist.PaperId and ord.AuthorId = hist.AuthorId
LEFT JOIN normalized_impact impact ON impact.PaperId = hist.PaperId 
LEFT JOIN level0_field_top f on hist.PaperId = f.PaperId
ORDER BY AuthorId, PaperId, Year
"""

client = bigquery.Client()

# Execute the query
random_seq = gen_random_sequence()
TEMP_TABLE_REF = f"ccnr-success.dmurray.temp_{random_seq}"

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
