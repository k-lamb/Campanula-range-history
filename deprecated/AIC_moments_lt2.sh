#!/usr/bin/bash
#SBATCH -J AIC # A single job name for the array
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH -o /home/ksl2za/automation_dadi/All/AIC/out_%j.out
#SBATCH -e /home/ksl2za/automation_dadi/All/AIC/error_%j.error
#SBATCH -p standard
#SBATCH -a lfg_lab

### sbatch --array=1-4 AIC_moments_lt2.sh

echo "began at"  `date`

#Load conda module
module purge
module load anaconda/2020.11-py3.8

#Load the metadata object into memory
metadata=/home/ksl2za/automation_dadi/All/AIC_metadat.txt #Address to the metadata. 

#Mining the metadata file
model=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $1 }' )

cd ~/automation_dadi/All/
python AIC_moments_lt2.py $model

#Print the time
echo "ended at"  `date`