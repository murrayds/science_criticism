#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=512GB
#SBATCH --ntasks 50
#SBATCH --time=4:00:00
#SBATCH --job-name=novelty_workflow
#SBATCH --partition=netsi_standard
#SBATCH -e error_%j.txt # Standard error file

make novelty-hpc
