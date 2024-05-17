#!/usr/bin/bash
#SBATCH -J am_65 # A single job name for the array
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH --mem=8G
#SBATCH -o /scratch/ksl2za/demo_infer/65_linreg/ancient_migration/out/am_out_%j.txt
#SBATCH -e /scratch/ksl2za/demo_infer/65_linreg/ancient_migration/error/am_error_%j.txt
#SBATCH -p standard
#SBATCH -a lfg_lab

### sbatch --array=1-3 am_extended.sh

#Load conda module
module purge
module load anaconda

#Activate dadi kernel
source activate dadi_env

#Load the metadata object into memory
metadata=/scratch/ksl2za/demo_infer/65_linreg/All_metadat_65_linreg.txt #Address to the metadata. 

#Pair_names=> the names of the two pools being considered. Separated by a predictable delimiter like "|" (do not use "_")!!!
#fs=> the address of the file containing the 2D SFS
#iterations=> number of times for the script to loop

#===> Want to fix SLURM_ARRAY_TASK_ID? this is useful for debugging <====
#=# SLURM_ARRAY_TASK_ID=120

#Mining the metadata file
iter=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $1 }' )
Pair=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $2 }' )
pop1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $3 }' )
pop2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $4 }' )
proj1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $5 }' )
proj2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $6 }' )

echo "Now Processing" $Pair
echo $Pair "is running" $iter "iterations"

#Run DADI
# this script takes 6 arguments
#iterations = sys.argv[1] ==> number of runs ... a number
#pair_name = sys.argv[2] ==> name of the pair -> $Pair
#pop1 = sys.argv[3] ==> pop1 name -> $pop1
#pop2 = sys.argv[4] ==> pop2 name -> $pop2
#projection1 = int(sys.argv[5]) ==> projection for pop1 -> $proj1
#projection2 = int(sys.argv[6]) ==> projection for pop2 -> $proj2

cd /scratch/ksl2za/demo_infer/65_linreg/ancient_migration
python am_extended.py \
$iter \
$Pair \
$pop1 \
$pop2 \
$proj1 \
$proj2
#De-Activate moments kernel
conda deactivate

#Print the time
echo "ended at"  `date`