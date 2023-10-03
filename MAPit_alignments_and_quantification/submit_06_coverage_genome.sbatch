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



module load samtools


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


samtools coverage --depth 0 ${INPUT}/${SAMPLE}.duprg.bam >  ${INPUT}/${SAMPLE}.duprg.coverage.txt


samtools depth  -a ${INPUT}/${SAMPLE}.duprg.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${INPUT}/${SAMPLE}.duprg.coverage_average_whole_genome.txt
samtools depth -a ${INPUT}/${SAMPLE}.duprg.bam | grep -v ChrC | grep -v ChrM | awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${INPUT}/${SAMPLE}.duprg.coverage_average_nuclear_genome.txt



