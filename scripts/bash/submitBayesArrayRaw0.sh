#!/bin/bash
#
#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=./outputTextRaw0/r_output_%J_%a.txt
#SBATCH --error=./errorTextRaw0/r_error_%J_%a.txt
#SBATCH --time=36:00:00
#SBATCH --job-name=no_change_point_bayes
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/pinesRep
#SBATCH --array=1-288
#SBATCH --array=1,2,3,9,10,11,17,18,19,25,26,27,73,74,75,81,82,83,89,90,91,97,98,99,145,146,147,153,154,155,161,162,163,169,170,171,217,218,219,225,226,227,233,234,235,241,242,243

module purge
setup_accre_software_stack
module load r/4.4.0

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
date
start=`date +%s`
Rscript ./scripts/rCode/sysCallRaw0.r ${SLURM_ARRAY_TASK_ID}
end=`date +%s`
runtime=$((end-start))
date
echo " Total runtime:" 
echo ${runtime}
echo "Done"