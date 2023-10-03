#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=make_bigWig
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/bigWig%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=8G


date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load ucsc/20210803


PROJ=/blue/concannon/share/jnewman/mingqi_arab
CHROMSIZE=/blue/concannon/share/jnewman/references/tair10/chrom.sizes

DATAIN=$PROJ/bwa_mem_aln_pe/mpileups_genome_tair10

cd $DATAIN

for BEDGRAPH in *.bedGraph
do
   BW=${BEDGRAPH/bedGraph/bw}
   grep -v VAR2 $BEDGRAPH > ./roz_bg
   bedGraphToBigWig ./roz_bg $CHROMSIZE $BW
done


