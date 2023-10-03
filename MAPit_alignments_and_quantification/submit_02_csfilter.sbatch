#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=csfilter
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.csfilter.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --array=1-36

### Max array: 24


#### Whole genome bisulfite sequencing analysis
#### Step 2) Incomplete bisulfite conversion filtering
####
#### In this script, the following will be performed:
####    (1) Remove unconverted/partially-converted reads
####    (2) Count reads after filtering
####
#### Note: although I am going to test out other mappers and methylation callers
#### (some of which can detect partial conversion), I decided to do this all up-front
#### so that the inputs for everything are the same. Plus, it might save time downstream


### Load modules

module load dibig_tools samtools

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



TRIMMEDFQ=$PROJ/qc_trim_output

TRIMR1=${NAME}_R1.trimmed.paired.fastq.gz
TRIMR2=${NAME}_R2.trimmed.paired.fastq.gz

#### (1) Run CSCaller
cd ${TRIMMEDFQ}
cscall -filter -1 ${TRIMR1} \
               -2 ${TRIMR2} \
               -amax 10000 \
               -fmax 1.0 \
               -s CG  GC

#### (2) Count reads after filtering

OUTLOG=$TRIMMEDFQ/postfilter_counts
   if [ ! -e $OUTLOG ]; then mkdir -p $OUTLOG; fi
OUTCOUNT=$OUTLOG/${NAME}-postfilter_counts.txt

NUMLINES=$(zcat ${TRIMMEDFQ}/good-${TRIMR1} | wc -l)
if [ "$NUMLINES" == "" ];
   then NUMSEQ=0
   else NUMSEQ=$(($NUMLINES / 4))
fi

echo -e "good-$FQIN1     $NUMSEQ" > $OUTCOUNT

NUMLINES=$(zcat ${TRIMMEDFQ}/good-${TRIMR2} | wc -l)
if [ "$NUMLINES" == "" ];
   then NUMSEQ=0
   else NUMSEQ=$(($NUMLINES / 4))
fi

echo -e "good-$FQIN2     $NUMSEQ" >> $OUTCOUNT


