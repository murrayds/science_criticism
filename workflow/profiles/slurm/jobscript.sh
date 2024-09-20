#!/bin/bash
#SBATCH --job-name={rule}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time={resources.runtime}:00
#SBATCH --mem={resources.mem_mb}M
#SBATCH --cpus-per-task={resources.cpus}
#SBATCH --partition=netsi_standard
#SBATCH --account=d.murray

# Load the anaconda module to make conda available
module load anaconda3/2022.05

echo "Jobscript called for rule {rule}, job ID: ${SLURM_JOB_ID}"
echo "SLURM partition: netsi_standard"
echo "SLURM account: d.murray"

# Run Snakemake command
{exec_job}
