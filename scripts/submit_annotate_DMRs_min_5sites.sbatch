#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=annot_DMR
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.annot_DMR.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48gb


module load homer

PROJ=/blue/concannon/share/jnewman/mingqi_arab
GTF=/blue/concannon/share/jnewman/references/tair10/Arabidopsis_thaliana.TAIR10.37.gtf
# Base annotations

annotatePeaks.pl ${PROJ}/DMRs_min_5_sites_for_HOMER.bed tair10 -gtf $GTF  > ${PROJ}/DMRs_min_5_sites_for_HOMER_annotation.txt

