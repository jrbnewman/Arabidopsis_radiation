#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=mark_dups
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.mark_dups.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
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
module load picard/2.18.3 samtools dibig_tools
export _JAVA_OPTIONS="-Xms10g -Xmx20g"

### Set directories

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

### (1) Mark duplicates with Picard

LOGS=$INPUT/picard_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi

picard MarkDuplicates \
       I=$INPUT/${SAMPLE}.sorted.bam \
       O=$INPUT/${SAMPLE}.markdup.bam \
       CREATE_INDEX=true \
       REMOVE_DUPLICATES=TRUE \
       ASSUME_SORTED=TRUE \
       METRICS_FILE=${LOGS}/${SAMPLE}-metrics.txt \
       USE_JDK_DEFLATER=true \
       USE_JDK_INFLATER=true \
       VALIDATION_STRINGENCY=SILENT \
       TMP_DIR=$ROZ

### (2) Add read groups




picard AddOrReplaceReadGroups \
       INPUT=$INPUT/${SAMPLE}.markdup.bam \
       OUTPUT=$INPUT/${SAMPLE}.duprg.bam \
       CREATE_INDEX=TRUE \
       RGPU=ICBRBI \
       RGID=ICBRBI \
       RGLB=${SAMPLE} \
       RGPL=illumina \
       RGSM=${SAMPLE} \
       USE_JDK_DEFLATER=true \
       USE_JDK_INFLATER=true \
       VALIDATION_STRINGENCY=SILENT \
       TMP_DIR=$ROZ


### (3) Count merged reads

COUNTSOUT=$INPUT/alignment_counts
if [ ! -e $COUNTSOUT ]; then mkdir -p $COUNTSOUT; fi

samtools idxstats $INPUT/${SAMPLE}.duprg.bam | grep -v ^\* | tcalc.py -q "do sum(C3)"> $COUNTSOUT/${SAMPLE}_BAM_counts.cnt

