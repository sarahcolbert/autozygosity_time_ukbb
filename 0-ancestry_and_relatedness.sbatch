#!/bin/bash

#SBATCH -J qc_and_prep
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --output=../logs/qc_and_prep-%a.out
#SBATCH --error=../logs/qc_and_prep-%a.err
#SBATCH --array=1-22,23-44,45-66,67-88,89-110,111-132

date
hostname

## set populations
pops=(afr amr csa eas mid eur)

## set job_id variable
job_id="${SLURM_ARRAY_TASK_ID}"
for i in {1..22}; do
    if [[ "$i" = "$job_id" ]]; then
       pop=${pops[0]}
        cc=$job_id
        break
    fi
done
for i in {23..44}; do
    if [[ "$i" = "$job_id" ]]; then
        pop=${pops[1]}
        cc=$(( $job_id - 22 ))
        break
    fi
done
for i in {45..66}; do
    if [[ "$i" = "$job_id" ]]; then
        pop=${pops[2]}
        cc=$(( $job_id - 44 ))
        break
    fi
done
for i in {67..88}; do
    if [[ "$i" = "$job_id" ]]; then
        pop=${pops[3]}
        cc=$(( $job_id - 66 ))
        break
    fi
done
for i in {89..110}; do
    if [[ "$i" = "$job_id" ]]; then
        pop=${pops[4]}
        cc=$(( $job_id - 88 ))
        break
    fi
done
for i in {111..132}; do
    if [[ "$i" = "$job_id" ]]; then
        pop=${pops[5]}
        cc=$(( $job_id - 110 ))
        break
    fi
done

## set variables
analysis_dir=#redacted
plink_dir="$analysis_dir"/data/raw/plink

## use plink to subset each ancestry and only keep unrelated samples
/redacted/plink1.9/plink --bed "$plink_dir"/ukb_cal_chr"$cc"_v2.bed \
--bim "$plink_dir"/ukb_snp_chr"$cc"_v2.bim \
--fam "$plink_dir"/ukb_cal_chr22_v2.fam \
--keep "$analysis_dir"/data/raw/panukb/panukb_"$pop"_ids_unrel.txt \
--make-bed \
--out "$analysis_dir"/data/raw/ancestry/"$pop"/"$pop"_unrel_chr"$cc"

date
