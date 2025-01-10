#!/bin/bash
#
#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=./outputTextRaw/r_output_%J_%a.txt
#SBATCH --error=./errorTextRaw/r_error_%J_%a.txt
#SBATCH --time=36:00:00
#SBATCH --job-name=change_point_bayes
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/pinesRep
#SBATCH --array=1-288

module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
date
start=`date +%s`
Rscript ./scripts/rCode/sysCallRaw.r ${SLURM_ARRAY_TASK_ID}
echo "Done"