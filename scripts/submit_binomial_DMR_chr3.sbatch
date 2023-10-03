#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=binom_DMR_3
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/binom_DMR_3.out
#SBATCH -t 16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=40G



date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load sas/9.4

PROJ=/blue/concannon/share/jnewman/mingqi_arab

sas -work $TMPDIR -memsize 40G -sysin $PROJ/07_run_binomial_DMR_model_arab.sas -set chrom 3 -log $PROJ/binomial_DMR_3.log
