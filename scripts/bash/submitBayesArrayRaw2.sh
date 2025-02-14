#!/bin/bash
#
#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=./outputTextRaw2/r_output_%J_%a.txt
#SBATCH --error=./errorTextRaw2/r_error_%J_%a.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=change_point_bayes
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/pinesRep
#SBATCH --array=1-288
#SBATCH --array=1,2,3,9,10,11,17,18,19,25,26,27,73,74,75,81,82,83,89,90,91,97,98,99,145,146,147,153,154,155,161,162,163,169,170,171,217,218,219,225,226,227,233,234,235,241,242,243


module purge
# module load GCC/11.3.0
# module load OpenMPI/4.1.4
# module load R/4.2.1

module load r/4.4.0

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
date
start=`date +%s`
Rscript ./scripts/rCode/sysCallRaw2.r ${SLURM_ARRAY_TASK_ID}
echo "Done"