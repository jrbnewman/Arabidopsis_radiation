#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=brass_DAR_A07
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/brass_DAR_longday_replication_A07.out
#SBATCH -t 16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=8G


date;hostname

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load sas/9.4

PROJ=/blue/concannon/share/jnewman/mingqi_arab

sas -work $TMPDIR -memsize 8G -sysin $PROJ/run_binomial_DAR_model_brass_LongDay_replication.sas -set chrom A07 -log $PROJ/binomial_DAR_longday_replication_A07.log



