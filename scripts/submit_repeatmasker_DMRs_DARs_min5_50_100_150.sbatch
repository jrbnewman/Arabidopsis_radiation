#!/bin/bash
#SBATCH --job-name=repeatMasker
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --cpus-per-task=4
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mem=12gb
#SBATCH -t 12:00:00
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/repeatMasker2_%j.out
#

module load repeatmasker samtools

# and here goes your command line

PROJ=/blue/concannon/share/jnewman/mingqi_arab

OUTPUT=$PROJ/repeatmasker

if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi


for REGION in min5_100bp_gc_01_dn min5_100bp_gc_01_up min5_100bp_gc_1_dn min5_100bp_gc_1_up min5_150bp_gc_01_dn min5_150bp_gc_01_up min5_150bp_gc_1_dn min5_150bp_gc_1_up min5_50bp_gc_01_dn min5_50bp_gc_01_up min5_50bp_gc_1_dn min5_50bp_gc_1_up
do

FASTA=$PROJ/dreme_sequences/${REGION}.fasta

if [ ! -e ${FASTA}.fai ]; then samtools faidx $FASTA; fi

SPECIES=arabidopsis

RepeatMasker $FASTA -s -species $SPECIES -html -gff -dir $OUTPUT -pa 4

done

