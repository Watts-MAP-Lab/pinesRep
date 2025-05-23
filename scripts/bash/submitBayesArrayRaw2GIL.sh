#!/bin/bash
#
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --output=./outputTextRaw2G/r_output_%J_%a.txt
#SBATCH --error=./errorTextRaw2G/r_error_%J_%a.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=twoGaus_change_point
#SBATCH --mail-user=adon.rosen@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/rosena/pinesRep
#SBATCH --array=1-48

module purge
setup_accre_software_stack
module load r/4.4.0

echo ${SLURM_ARRAY_TASK_ID}
echo "Submitting job"
date
start=`date +%s`
Rscript ./scripts/rCode/sysCallRaw2GausInvLogit.r ${SLURM_ARRAY_TASK_ID}
end=`date +%s`
runtime=$((end-start))
date
echo " Total runtime:" 
echo ${runtime}
echo "Done"