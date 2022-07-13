#!/bin/bash
#SBATCH --job-name=roh-%a
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=../logs/roh-%a.out
#SBATCH --error=../logs/roh-%a.err
#SBATCH --array=0-5

date
hostname

## set job id variable
job_id=$SLURM_ARRAY_TASK_ID

## select which population to process, based on array id
pops=(afr amr csa eas mid eur)
pop=${pops["$job_id"]}

## set variables
qc_data=#redacted/data/processed/"$pop"/"$pop"_qc_unrel_merged
out_dir=#redacted/data/processed/roh_estimates/

## call ROHs in plink
/redacted/plink1.9/plink --bfile ${qc_data} \
--homozyg-window-snp 50 \
--homozyg-snp 50  \
--homozyg-kb 1500  \
--homozyg-gap 1000  \
--homozyg-density 50 \
--homozyg-window-missing 5 \
--homozyg-window-het 1 \
--out ${out_dir}${pop}_roh

## estimate F_uni in plink
/ref/aalab/software/plink1.9/plink --bfile ${qc_data} \
--autosome \
--ibc \
--out ${out_dir}${pop}_funi

date
