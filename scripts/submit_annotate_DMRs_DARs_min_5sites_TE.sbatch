#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=annot_bed
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.annot_bed.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48gb


module load homer

PROJ=/blue/concannon/share/jnewman/mingqi_arab
GTF=$PROJ/TAIR10_transposable_elements.noChr.gtf
# Base annotations

annotatePeaks.pl ${PROJ}/DMRs_min_5_sites_for_HOMER.bed tair10 -gtf $GTF  > ${PROJ}/at_DMR_to_TE_min_5_sites.txt
annotatePeaks.pl ${PROJ}/DARs_min_5_sites_for_HOMER.bed tair10 -gtf $GTF  > ${PROJ}/at_DAR_to_TE_min_5_sites.txt


