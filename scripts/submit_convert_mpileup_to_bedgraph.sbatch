#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=make_bedGraphs
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/make_bedgraphs_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=40G


date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load sas/9.4


PROJ=/blue/concannon/share/jnewman/mingqi_arab

DATAIN=$PROJ/bwa_mem_aln_pe/mpileups_genome_tair10

cd $DATAIN

for MPILEUP in *.mpileup
do
   TSV=${MPILEUP/mpileup/tsv}
   cut -f1,2,4 $MPILEUP > $TSV
done

sas -work $TMPDIR -memsize 40G -sysin $PROJ/scripts/mpileup_to_bedgraph.sas


