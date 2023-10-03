#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=sort_bam
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.sort_index_bam.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --array=1-12

## max array 8

#### Whole genome bisulfite sequencing analysis
#### Step 4) Merge, sort, and index BAM files
####
#### In this script, the following will be performed:
####    (1) BAM files for samples from each sample are merged into one BAM file
####    (2) Merged BAM files are then sorted and indexed
####
#### BSmap version
####

### Load modules
module load samtools bamtools





### Set directories

PROJ=/blue/concannon/share/jnewman/mingqi_arab

# Design file
    DESIGN_FILE=$PROJ/design_files/sample_design_file_nofq.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

        GROUP=${ARRAY[0]}
        SAMPLE=${ARRAY[1]}

INPUT=$PROJ/BSmap_output

MERGEOUT=$INPUT/merged_bams
if [ ! -e $MERGEOUT ]; then mkdir -p $MERGEOUT; fi

### (1) Merge BAM files for each sample

bamtools merge -out $MERGEOUT/${SAMPLE}.bam \
               -in $INPUT/${SAMPLE}_r1.bam \
               -in $INPUT/${SAMPLE}_r2.bam \
               -in $INPUT/${SAMPLE}_r3.bam

### (2) Sort and index merged BAM files

samtools view -b -F 4 $MERGEOUT/${SAMPLE}.bam | samtools sort -@ 4 -O bam -T $MERGEOUT/tmp_${SAMPLE}_${SLURM_JOB_ID}.tmp - > $MERGEOUT/${SAMPLE}.sorted.bam

samtools index $MERGEOUT/${SAMPLE}.sorted.bam

