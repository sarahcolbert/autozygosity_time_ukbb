#!/bin/bash
#SBATCH --job-name=qc-%a
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=../logs/qc-%a.out
#SBATCH --error=../logs/qc-%a.err
#SBATCH --array=0-4

date
hostname

## set job id variable
job_id=$SLURM_ARRAY_TASK_ID

## select which population to process, based on array id
pops=(afr amr csa eas mid eur)
pop=${pops["$job_id"]}

## set variables
data_dir=#redacted
plink_dir="$data_dir"/raw/ancestry


/redacted/plink1.9/plink --bfile "$data_dir"/processed/"$pop"/"$pop"_qc1_unrel_merged \
--mind 0.03 \
--make-bed \
--out "$data_dir"/processed/"$pop"/"$pop"_qc_unrel_merged

date
