#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=scanGenome
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.scanMotif.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48gb


module load homer

PROJ=/blue/concannon/share/jnewman/mingqi_arab
MOTIF=$PROJ/top_motifs.motif
GENOME=tair10
OUT=$PROJ/top_motifs_genomewide.txt

scanMotifGenomeWide.pl $MOTIF $GENOME > $OUT

