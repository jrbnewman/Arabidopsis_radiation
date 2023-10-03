#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=annot_bed
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.annot_bed.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48gb


module load homer/4.10

PROJ=/blue/concannon/share/jnewman/mingqi_arab
GFF=/blue/concannon/share/jnewman/references/Brapa_15/R500_v1.6.gtf

GTF=$PROJ/Br_transposon_RepeatModeler.gtf
FASTA=/blue/concannon/share/jnewman/references/Brapa_15/Brapa_R500_V1.2.fasta
# Base annotations

#loadGenome.pl -name Brapa -org Brapa -fasta $FASTA -gff3 $GTF

#parseGTF.pl $GFF ann > $PROJ/Br_gff_annot.txt
#grep "exon" $PROJ/Br_gff_annot.txt > $PROJ/Br_gff_annot2.txt
#grep -v "exon" $PROJ/Br_gff_annot.txt >> $PROJ/Br_gff_annot2.txt

#parseGTF.pl $GTF ann > $PROJ/Br_TE_annot.txt
#grep "exon" $PROJ/Br_TE_annot.txt > $PROJ/Br_TE_annot2.txt
#grep -v "exon" $PROJ/Br_TE_annot.txt >> $PROJ/Br_TE_annot2.txt

#assignGenomeAnnotation $PROJ/Br_gff_annot2.txt $PROJ/Br_gff_annot2.txt -prioritize $PROJ/Br_gff_annot.final.txt  > $PROJ/gff_stats.txt
#assignGenomeAnnotation $PROJ/Br_TE_annot2.txt $PROJ/Br_TE_annot2.txt -prioritize $PROJ/Br_TE_annot.final.txt  > $PROJ/stats.txt

BED1=$PROJ/Br_merged_superregions_all_comparisons.bed
BED2=$PROJ/Br_merged_subregions_all_comparisons.bed
BED3=$PROJ/Br_merged_subregions_all_comparisons_200bp_ext.bed



annotatePeaks.pl ${BED1} $FASTA -gtf $GFF  > ${PROJ}/brassica_radiation_superregions_all_comparisons_HOMER_custom.txt
annotatePeaks.pl ${BED2} $FASTA -gtf $GFF  > ${PROJ}/brassica_radiation_merged_subregions_all_comparisons_HOMER_custom.txt
annotatePeaks.pl ${BED3} $FASTA -gtf $GFF  > ${PROJ}/brassica_radiation_merged_subregions_all_comparisons_200bp_ext_HOMER_custom.txt

#annotatePeaks.pl ${BED1} $FASTA -ann $PROJ/Br_gff_annot.final.txt > ${PROJ}/brassica_radiation_superregions_all_comparisons_HOMER.txt
#annotatePeaks.pl ${BED2} $FASTA -ann $PROJ/Br_gff_annot.final.txt > ${PROJ}/brassica_radiation_merged_subregions_all_comparisons_HOMER.txt
#annotatePeaks.pl ${BED3} $FASTA -ann $PROJ/Br_gff_annot.final.txt > ${PROJ}/brassica_radiation_merged_subregions_all_comparisons_200bp_ext_HOMER.txt


#annotatePeaks.pl ${BED1} $FASTA -ann $PROJ/Br_TE_annot.final.txt > ${PROJ}/brassica_radiation_superregions_all_comparisons_HOMER_TE.txt
#annotatePeaks.pl ${BED2} $FASTA -ann $PROJ/Br_TE_annot.final.txt > ${PROJ}/brassica_radiation_merged_subregions_all_comparisons_HOMER_TE.txt
#annotatePeaks.pl ${BED3} $FASTA -ann $PROJ/Br_TE_annot.final.txt > ${PROJ}/brassica_radiation_merged_subregions_all_comparisons_200bp_ext_HOMER_TE.txt





