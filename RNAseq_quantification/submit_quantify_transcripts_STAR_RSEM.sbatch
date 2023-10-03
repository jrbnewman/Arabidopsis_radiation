#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --job-name=arab_quant_xs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOG_out.%j.%A.%a.out
#SBATCH -t 1-00:00:00
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --mem 48gb
#SBATCH --array=2-36%18

## 36 samples

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export ROZ=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# ICBR STAR-RSEM pipeline for Mingqi's arabidopsis data

## Load modules and set Java memory size

module load R/3.2.0 rsem/1.2.28 java/1.8.0_31 samtools/1.3.1 htseq/0.6.1p1 picard/2.5.0 trimmomatic/0.32 gcc/5.2.0 star/2.5.2b
export _JAVA_OPTIONS="-Xms8000m -Xmx16000m"

## Set directories

# Project folder
PROJ=/ufrc/concannon/share/jnewman/mingqi_arab

# Location of original data
ORIG=$PROJ/original_data

# Location to put concatenated FQ files
FASTQ=$PROJ/fastq_files
  if [ ! -e $FASTQ ]; then mkdir -p $FASTQ; fi

# Location to put trimmed FQs
FQTRIM=$PROJ/trimmed_fastq_files
   if [ ! -e $FQTRIM ]; then mkdir -p $FQTRIM; fi


## Design file
DESIGNFILE=$PROJ/design_files/arabidopsis_samples.csv


     DESIGN=$(cat ${DESIGNFILE} | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}
## Cat FQ files

cd $ORIG
FASTQCAT_R1=$FASTQ/${SAMPLE}_R1.fastq
FASTQCAT_R2=$FASTQ/${SAMPLE}_R2.fastq
FLAG=0
for READ1 in ${SAMPLE}*R1_001.fastq.gz;
    do
    if [ $FLAG == 0 ]
    then
    gunzip -c ${READ1} > ${FASTQCAT_R1}
    FLAG=1
    else
    gunzip -c ${READ1} >> ${FASTQCAT_R1}
    fi
done

FLAG=0
for READ2 in ${SAMPLE}*R2_001.fastq.gz;
    do 
    if [ $FLAG == 0 ] 
    then
    gunzip -c ${READ2} > ${FASTQCAT_R2}
    FLAG=1
    else
    gunzip -c ${READ2} >> ${FASTQCAT_R2}
    fi
done

## (1) Run trimmomatic

TRIMDIR=$PROJ/software/Trimmomatic-0.35

#java -jar ${TRIMDIR}/trimmomatic-0.35.jar \
#          PE \
#          -threads 8 \
#          -phred33 \
#          ${FASTQCAT_R1} \
#          ${FASTQCAT_R2} \
#          ${FQTRIM}/${SAMPLE}_paired_1.fq.gz \
 #         ${FQTRIM}/${SAMPLE}_unpaired_1.fq.gz \
 #         ${FQTRIM}/${SAMPLE}_paired_2.fq.gz \
  #        ${FQTRIM}/${SAMPLE}_unpaired_2.fq.gz \
   #        LEADING:20 \
    #       TRAILING:20 \
     #      SLIDINGWINDOW:4:15 \
      #     MINLEN:80

## (2) Align to genome with STAR and output transcriptome SAM file

# Set STAR output folder
STAROUT=$PROJ/aln_star_tair10
   if [ ! -e $STAROUT ]; then mkdir -p $STAROUT; fi

# Set STAR reference
#STARREF=$PROJ/references/TAIR10_STAR_2.5.2b
STARREF=$PROJ/references
echo "Running STAR"
STAR --runThreadN 8 \
     --genomeDir ${STARREF} \
     --genomeLoad LoadAndRemove \
     --readFilesIn ${FQTRIM}/${SAMPLE}_paired_1.fq.gz \
                   ${FQTRIM}/${SAMPLE}_paired_2.fq.gz \
     --readFilesCommand 'zcat' \
     --outFileNamePrefix ${STAROUT}/${SAMPLE}. \
     --outTmpDir $ROZ/tmp_${SLURM_JOB_ID} \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 48000000000 \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 3 \
     --alignIntronMin 7 \
     --alignIntronMax 40000 \
     --alignMatesGapMax 40000 \
     --outSAMattributes Standard \
     --outFilterIntronMotifs RemoveNoncanonical \
     --outSAMstrandField intronMotif \
     --quantMode TranscriptomeSAM


## (3) Quantify transcripts with RSEM
echo "Running RSEM"
# Set RSEM output folder
RSEMOUT=$PROJ/rsem_output/${SAMPLE}_expr
   if [ ! -e $RSEMOUT ]; then mkdir -p $RSEMOUT; fi

# Set RSEM reference
RSEMREF=$PROJ/references/TAIR10.rsem.tg

rsem-calculate-expression -p 8 \
                          --alignments \
                          --paired-end \
                          ${STAROUT}/${SAMPLE}*.Aligned.toTranscriptome.out.bam \
                          ${RSEMREF} \
                          ${RSEMOUT}


