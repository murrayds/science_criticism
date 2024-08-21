from google.cloud import bigquery
from google.api_core.exceptions import NotFound
# Initialize Google Cloud Storage client
from google.cloud import storage
from concurrent.futures import ThreadPoolExecutor

import os
import shutil
import glob
import gzip

# Retreive the year wildcard passed into the file
YEAR = int(snakemake.wildcards.year)


# Check if the table exists
BUILD_TABLE_QUERY = """
WITH target_journals AS (
  SELECT 
    JournalID
  FROM `ccnr-success.mag.Journals`
  WHERE PaperCount > 1000
),
papers AS (
  SELECT
    p.PaperId as CitingPaperId,
    p.Year as CitingPaperYear
  FROM `ccnr-success.mag.Papers` p
  INNER JOIN target_journals j on j.JournalId = p.JournalId
  WHERE DocType = "Journal"
  AND Doi IS NOT NULL
  AND ReferenceCount >= 5
),
refs AS (
  select
    p.*,
    r.PaperReferenceId,
    p2.JournalId as ReferenceJournalId,
    p2.Year as ReferenceYear
  FROM `ccnr-success.mag.PaperReferences` as r
  INNER JOIN papers as p ON p.CitingPaperId = r.PaperId
  LEFT JOIN `ccnr-success.mag.Papers` p2 on p2.PaperId = r.PaperReferenceId
  INNER JOIN target_journals j on j.JournalId = p2.JournalId
  WHERE p2.Doi IS NOT NULL
  AND p2.DocType = "Journal"
)
SELECT 
  CitingPaperId,
  CitingPaperYear,
  PaperReferenceId,
  ReferenceJournalId,
  ReferenceYear
FROM refs
"""

client = bigquery.Client()

print(f"Beggining select table query.")

# First, check to see if the table exists...if not create it
MAIN_REF_TABLE = "ccnr-success.mag.dmurray_selected_refs"
try:
    client.get_table(MAIN_REF_TABLE)
    print(f"Table {MAIN_REF_TABLE} already exists.")
except NotFound:
    print(f"Table {MAIN_REF_TABLE} is not found. Creating it now.")
    
    # Set up the query job configuration
    job_config = bigquery.QueryJobConfig(
        destination=MAIN_REF_TABLE,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
        use_query_cache=False
    )

    # Run the query and save results to the new table
    query_job = client.query(BUILD_TABLE_QUERY, job_config=job_config)
    query_job.result()  # Wait for the query to complete

    print(f"Query results exported to: {MAIN_REF_TABLE}")




# To make downloading and downstream processing easier, we will query only 
# one year at a time...
SELECT_YEAR_QUERY = f"""
SELECT 
    CitingPaperId,
    ReferenceJournalId,
    ReferenceYear
FROM {MAIN_REF_TABLE}
WHERE CitingPaperYear = {YEAR}
"""


# Execute the query
query_job = client.query(SELECT_YEAR_QUERY)
SELECTED_TABLE_REF = f"ccnr-success.dmurray.selected_refs_{YEAR}"

# Set up the query job configuration
job_config = bigquery.QueryJobConfig(
    destination=SELECTED_TABLE_REF,
    write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
    use_query_cache=False
)

# Run the query and save results to the new table
query_job = client.query(SELECT_YEAR_QUERY, job_config=job_config)
query_job.result()  # Wait for the query to complete

print(f"Query results exported to table {SELECTED_TABLE_REF}")

# Now we will export the temporary table to a google cloud storagebucket...
destination_uri = "gs://{}/{}".format("dmurray_temp", f"selected_refs_{YEAR}_*.csv.gz")

# setup the job config
job_config = bigquery.job.ExtractJobConfig()
job_config.compression = bigquery.Compression.GZIP

extract_job = client.extract_table(
    SELECTED_TABLE_REF,
    destination_uri,
    # Location must match that of the source table.
    location="US",
    job_config = job_config
)  # API request
extract_job.result()  # Waits for job to complete.

print("Export job completed...")

#
# Now extract the files from the google cloud bucket...
#

# Initialize Google Cloud Storage client
storage_client = storage.Client()

# Get the bucket
# Assume that one already exists to hold the results...
bucket = storage_client.get_bucket("dmurray_temp")


# Create a local directory to store the downloaded files
local_dir = f"temp_{YEAR}"
if os.path.exists(local_dir):
    shutil.rmtree(local_dir)
# Create the directory
os.makedirs(local_dir)

# List all blobs with the specified prefix
blobs = list(bucket.list_blobs(prefix=f"selected_refs_{YEAR}_"))

# Download all matching files in parallel
def download_blob(blob):
    if blob.name.endswith(".csv.gz"):
        local_path = os.path.join(local_dir, os.path.basename(blob.name))
        blob.download_to_filename(local_path)
        print(f"Downloaded {blob.name}")
        blob.delete()
        print(f"Deleted {blob.name} from bucket")

with ThreadPoolExecutor(max_workers=10) as executor:
    executor.map(download_blob, blobs)

print("All files have been downloaded and removed from the bucket.")

# Now unzip each of these compressed files
print("Unzipping downloaded files...")
for filename in os.listdir(local_dir):
    if filename.endswith(".csv.gz"):
        gz_path = os.path.join(local_dir, filename)
        csv_path = os.path.join(local_dir, filename[:-3])  # Remove .gz extension
        
        with gzip.open(gz_path, 'rb') as gz_file:
            with open(csv_path, 'wb') as csv_file:
                csv_file.write(gz_file.read())
        
        # Remove the original .gz file
        os.remove(gz_path)

print("All files have been unzipped and .gz files removed.")

#
# Combine all .csv files into a single file...
#
allFiles = glob.glob(f"{local_dir}/*.csv")
allFiles.sort()  # glob lacks reliable ordering, so impose your own if output order matters
with open(snakemake.output[0], 'wb') as outfile:
    for i, fname in enumerate(allFiles):
        with open(fname, 'rb') as infile:
            if i != 0:
                infile.readline()  # Throw away header on all but first file
            # Block copy rest of file from input to output without parsing
            shutil.copyfileobj(infile, outfile)
            print(fname + " has been imported.")

#
# CLEANUP
#

# Delete the local directory that contained the individual .csv files
shutil.rmtree(local_dir)
print(f"Temporary directory {local_dir} has been deleted.")

# Delete the temporary GBQ tables...
client = bigquery.Client()
client.delete_table(SELECTED_TABLE_REF, not_found_ok=True)
print(f"Table '{SELECTED_TABLE_REF}' deleted.")


