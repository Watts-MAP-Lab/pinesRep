#!/bin/bash
#
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --output=./outputTextRaw0/r_output_%J_%a.txt
#SBATCH --error=./errorTextRaw0/r_error_%J_%a.txt
#SBATCH --time=60:00:00
#SBATCH --job-name=no_change_point_bayes
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
Rscript ./scripts/rCode/sysCallRaw0.r ${SLURM_ARRAY_TASK_ID}
end=`date +%s`
runtime=$((end-start))
date
echo " Total runtime:" 
echo ${runtime}
echo "Done"