#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=ercc
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/ercc_%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --array=2-721%20


###361
#


date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load bowtie/0.12.9
module load python/2.7.6

PROJ=/blue/concannon/share/jnewman/mingqi_arab
ORIG=$PROJ/original_data
FASTQ=$PROJ/unzipped
mkdir -p $FASTQ

REF=/blue/mcintyre/share/references/ERCC_Ambion/ERCC92

## Design file
DESIGN_FILE=$PROJ/design_files/qc_design_file.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

READ=${ARRAY[2]}
NAME=${ARRAY[1]}

gunzip -c $ORIG/${READ}.gz > $FASTQ/${READ}

NUMPROCS=2
OUTPUT=$PROJ/ercc_aln_se
mkdir -p $OUTPUT


LOGS=$OUTPUT/job_logs
ALN_ERROR_LOG=$LOGS/size_errors.txt
if [ ! -d $LOGS ]; then mkdir -p $LOGS; fi

MYLOG=$LOGS/${NAME}.log
printf "`date` $NAME PBS_ARRAYID:$PBS_ARRAYID HOSTNAME:$HOSTNAME\n" > $MYLOG


ALNLOGS=$OUTPUT/aln_logs
if [ ! -d $ALNLOGS ]; then mkdir -p $ALNLOGS; fi

UNALNDIR=$OUTPUT/unaln_reads
if [ ! -d $UNALNDIR ]; then mkdir -p $UNALNDIR; fi



source /blue/concannon/share/jnewman/qc_pipeline/scripts/alignment_functions.sh

READS=$FASTQ/${READ}

qual=`python /blue/mcintyre/share/python.git/identify_quality.py -i $READS`

if [ $qual == "phred64" ]; 
then
btqual="--phred64-quals"
lastqual="3"
else
# change to sanger format which is what all new illumina data is
btqual="--phred33-quals"
lastqual="1"
fi
bowtie_se_uniq

cat *.sam > $OUTPUT/${NAME}.sam 2>> $MYLOG

echo "START Combine Unaln FQ">>$MYLOG
cat *_unaln_bt.fq >$UNALNDIR/${NAME}_unaln.fq 2>> $MYLOG

rm $UNALNDIR/${NAME}_unaln.fq
rm $OUTPUT/${NAME}.sam
rm $FASTQ/${READ}

