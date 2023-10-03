#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=trim_qc
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.trim_qc.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --array=1-36

### Max array: 24


#### Whole genome bisulfite sequencing analysis
#### Step 1) Adapter-trimming and QC
####
#### In this script, the following will be performed:
####    (1) Count reads in FQ file
####    (2) Trim adapters using Trimmomatic
####    (3) Count reads in trimmed FQ file
####    (4) FASTQC on samples (pre and post trimming)

### Load modules

module load fastqc java trimmomatic

### Set directories

PROJ=/blue/concannon/share/jnewman/mingqi_arab


# Design file
    DESIGN_FILE=$PROJ/design_files/sample_design_file.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

        GROUP=${ARRAY[0]}
        SAMPLE=${ARRAY[1]}
        REP=${ARRAY[2]}
        ORIG=${ARRAY[3]}
        READ1=${ARRAY[4]}
        READ2=${ARRAY[5]}

    NAME=${SAMPLE}_${REP}


    FQIN1=${ORIG}/${READ1}
    FQIN2=${ORIG}/${READ2}

OUTDIR=$PROJ/qc_trim_output
   if [ ! -e $OUTLOG ]; then mkdir -p $OUTLOG; fi
OUTLOG=$OUTDIR/pretrim_counts
   if [ ! -e $OUTLOG ]; then mkdir -p $OUTLOG; fi

### (1) Count reads in input FQ file
OUTCOUNT=$OUTLOG/${NAME}-pretrim_counts.txt

NUMLINES=$(zcat ${FQIN1} | wc -l)
if [ "$NUMLINES" == "" ];
   then NUMSEQ=0
   else NUMSEQ=$((${NUMLINES} / 4))
fi

echo -e "$FQIN1     $NUMSEQ" > $OUTCOUNT

NUMLINES=$(zcat ${FQIN2} | wc -l)
if [ "$NUMLINES" == "" ];
   then NUMSEQ=0
   else NUMSEQ=$((${NUMLINES} / 4))
fi


echo -e "$FQIN2     $NUMSEQ" >> $OUTCOUNT

### (2) Pre-trim FASTQC

FASTQCOUT=$OUTDIR/fastqc_pretrim
   if [ ! -e $FASTQCOUT ]; then mkdir -p $FASTQCOUT; fi

fastqc -o $FASTQCOUT $FQIN1
fastqc -o $FASTQCOUT $FQIN2

### (3) Trimmomatic

ADAPT=TruSeq3-PE-2.fa
MINLEN=50
TRIMOUT=$OUTDIR

TRIMLOG=$OUTDIR/trimmomatic_logs
   if [ ! -e $TRIMLOG ]; then mkdir -p $TRIMLOG; fi

TRIMLOGOUT=$TRIMLOG/${NAME}_${GROUP}_trimmomatic.log

TRIMPAIRR1=$OUTDIR/${NAME}_R1.trimmed.paired.fastq.gz
TRIMUNPAIRR1=$OUTDIR/${NAME}_R1.trimmed.unpaired.fastq.gz
TRIMPAIRR2=$OUTDIR/${NAME}_R2.trimmed.paired.fastq.gz
TRIMUNPAIRR2=$OUTDIR/${NAME}_R2.trimmed.unpaired.fastq.gz

trimmomatic PE -threads 2 -trimlog ${TRIMLOGOUT} \
  ${FQIN1} ${FQIN2} $TRIMPAIRR1 $TRIMUNPAIRR1 $TRIMPAIRR2 $TRIMUNPAIRR2 \
  ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/${ADAPT}:2:30:10:2:true \
  LEADING:3 \
  TRAILING:3 \
  MAXINFO:$MINLEN:0.3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:$MINLEN



### (4) Post-trimming FASTQC

FASTQCOUT2=$OUTDIR/fastqc_posttrim
   if [ ! -e $FASTQCOUT2 ]; then mkdir -p $FASTQCOUT2; fi

fastqc -o $FASTQCOUT2 $TRIMPAIRR1
fastqc -o $FASTQCOUT2 $TRIMPAIRR2


### (5) Post-trimming counting

OUTLOG2=$OUTDIR/posttrim_counts
   if [ ! -e $OUTLOG2 ]; then mkdir -p $OUTLOG2; fi

### (1) Count reads in input FQ file
OUTCOUNT2=$OUTLOG2/${NAME}-posttrim_counts.txt

NUMLINES=$(zcat ${TRIMPAIRR1} | wc -l)
if [ "$NUMLINES" == "" ];
   then NUMSEQ=0
   else NUMSEQ=$(($NUMLINES / 4))
fi

echo -e "$TRIMPAIRR1     $NUMSEQ" > $OUTCOUNT2

NUMLINES=$(zcat ${TRIMPAIRR2} | wc -l)
if [ "$NUMLINES" == "" ];
   then NUMSEQ=0
   else NUMSEQ=$(($NUMLINES / 4))
fi

echo -e "$TRIMPAIRR2     $NUMSEQ" >> $OUTCOUNT2


