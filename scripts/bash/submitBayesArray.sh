#!/bin/bash
#
#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=./outputText/r_output_%J_%a.txt
#SBATCH --error=./errorText/r_error_%J_%a.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=change_point_bayes
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/pinesRep
#SBATCH --array=1-24

module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
date
start=`date +%s`
Rscript ./scripts/rCode/sysCall.r ${SLURM_ARRAY_TASK_ID}
echo "Done"
end=`date +%s`
runtime=$((end-start))
date
echo " Total runtime:" 
echo ${runtime}