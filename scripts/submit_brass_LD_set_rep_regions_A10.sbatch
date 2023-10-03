#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=brass_DAR_A10
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/brass_DAR_longday_replication_A10.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G


date;hostname

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load sas/9.4

PROJ=/blue/concannon/share/jnewman/mingqi_arab


sas -work $TMPDIR -memsize 20G -sysin $PROJ/run_set_replication_regions.sas -set chrom A10 -log $PROJ/run_set_replication_regions_A10.log

