#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=mark_dups
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.mark_dups.%j.%A.%a.out
#SBATCH -t 1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --array=1-12

## max array 8

#### Whole genome bisulfite sequencing analysis
#### Step 5) Mark PCR duplicates
####
#### In this script, the following will be performed:
####    (1) Using Picard to mark PCR duplicates in sorted BAM files (merged)
####    (2) Add read groups
####    (3) Use `samtools flagstat` to count SAM flags from output BAM file
####    (4) Count number of merge reads
####

### Load modules

module samtools
BED=/blue/concannon/share/jnewman/mingqi_arab/BSmap_output/merged_bams/nuclear_regions.bed




PROJ=/blue/concannon/share/jnewman/mingqi_arab


# Design file
    DESIGN_FILE=$PROJ/design_files/sample_design_file_nofq.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

        GROUP=${ARRAY[0]}
        SAMPLE=${ARRAY[1]}

INPUT=$PROJ/BSmap_output/merged_bams

ROZ=$PROJ/roz/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
if [ ! -e $ROZ ]; then mkdir -p $ROZ; fi


samtools view -o $INPUT/${SAMPLE}.nuclear.bam -L $BED $INPUT/${SAMPLE}.sorted.bam
samtools view $INPUT/${SAMPLE}.sorted.bam | cut -f2 | sort | uniq -c > $INPUT/${SAMPLE}.sorted.all.flagFreq.txt
samtools view $INPUT/${SAMPLE}.nuclear.bam | cut -f2 | sort | uniq -c > $INPUT/${SAMPLE}.sorted.nuclear.flagFreq.txt
rm $INPUT/${SAMPLE}.nuclear.bam



samtools view -o $INPUT/${SAMPLE}.nuclear.bam -L $BED $INPUT/${SAMPLE}.markdup.bam
samtools view $INPUT/${SAMPLE}.markdup.bam | cut -f2 | sort | uniq -c > $INPUT/${SAMPLE}.markdup.all.flagFreq.txt
samtools view $INPUT/${SAMPLE}.nuclear.bam | cut -f2 | sort | uniq -c > $INPUT/${SAMPLE}.markdup.nuclear.flagFreq.txt
rm $INPUT/${SAMPLE}.nuclear.bam



samtools view -o $INPUT/${SAMPLE}.nuclear.bam -L $BED $INPUT/${SAMPLE}.duprg.bam
samtools view $INPUT/${SAMPLE}.duprg.bam | cut -f2 | sort | uniq -c > $INPUT/${SAMPLE}.duprg.all.flagFreq.txt
samtools view $INPUT/${SAMPLE}.nuclear.bam | cut -f2 | sort | uniq -c > $INPUT/${SAMPLE}.duprg.nuclear.flagFreq.txt
rm $INPUT/${SAMPLE}.nuclear.bam

