#!/bin/bash
#SBATCH --job-name=merge-%a
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --output=../logs/merge-%a.out
#SBATCH --error=../logs/merge-%a.err
#SBATCH --array=0-5

date
hostname

## set job id variable
job_id=$SLURM_ARRAY_TASK_ID

## select which population to process, based on array id
pops=(afr amr csa eas eur mid)
pop=${pops["$job_id"]}


## set variables
analysis_dir=#redacted
plink_dir="$analysis_dir"/data/raw/plink

cd "$analysis_dir"/data/raw/ancestry/"$pop"

rm mergelist.txt

## merge with plink
for i in {1..22}
do
echo "$pop"_unrel_chr$i >> mergelist.txt
done

/redacted/plink1.9/plink \
--merge-list mergelist.txt \
--make-bed \
--out "$pop"_unrel_merged

date
