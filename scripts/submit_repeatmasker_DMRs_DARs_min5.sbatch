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
#SBATCH --array=1-40
#

module load repeatmasker samtools

# and here goes your command line

PROJ=/blue/concannon/share/jnewman/mingqi_arab

OUTPUT=$PROJ/repeatmasker

if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

DESIGN_FILE=$PROJ/design_files/min5_seq.txt
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

REGION=${ARRAY[0]}

FASTA=$PROJ/dreme_sequences/${REGION}.fasta

if [ ! -e ${FASTA}.fai ]; then samtools faidx $FASTA; fi

SPECIES=arabidopsis

RepeatMasker $FASTA -s -species $SPECIES -html -gff -dir $OUTPUT -pa 4

