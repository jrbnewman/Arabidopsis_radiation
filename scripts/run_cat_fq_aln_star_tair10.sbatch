#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=star_aln
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/star_aln_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --array=2-36


###36
#


date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load star

### 64 samples total

### FQsplitdups on raw reads
#Set directories

PROJ=/blue/concannon/share/jnewman/mingqi_arab
FASTQ=$PROJ/unzipped/fastqSplitDups

OUTPUT=$PROJ/star_aln

mkdir -p $OUTPUT

## Design file
DESIGN_FILE=$PROJ/design_files/arabidopsis_samples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}

cd ${FASTQ}

FQCAT=$PROJ/fq_cat
mkdir $FQCAT

READ1=$FQCAT/${SAMPLEID}_R1.fastq
READ2=$FQCAT/${SAMPLEID}_R2.fastq

for FILE in ${SAMPLEID}-*_R1_001_distinct.fq
        do
        cat $FILE >> $READ1
        cat ${FILE/_R1_/_R2_} >> $READ2
done



STAR --runThreadN 4 \
     --genomeDir $PROJ/references/TAIR_STAR \
     --sjdbGTFfile $PROJ/references/TAIR_STAR/Arabidopsis_thaliana.TAIR10.37.gtf \
     --sjdbOverhang 100 \
     --readFilesIn $READ1 $READ2 \
     --outSAMtype BAM Unsorted \
     --winAnchorMultimapNmax 200 \
     --outFilterMultimapNmax 100 \
     --outFileNamePrefix ${OUTPUT}/${SAMPLEID}_

