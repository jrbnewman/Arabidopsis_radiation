#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=ercc
#SBATCH --account=concannon
#SBATCH --qos=brusko-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/ercc_%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --array=2-36

##36

date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load bowtie/0.12.9
module load python/2.7.6

PROJ=/blue/concannon/share/jnewman/mingqi_arab
ORIG=$PROJ/original_data
FQGEO=$PROJ/FQ_merged_for_GEO
mkdir -p $FQGEO


REF=/blue/mcintyre/share/references/ERCC_Ambion/ERCC92

DESIGN_FILE=$PROJ/design_files/arabidopsis_samples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"
SAMPLE=${ARRAY[0]}

if [ -f $FQGEO/${SAMPLE}_R1.fastq ]; then
rm $FQGEO/${SAMPLE}_R1.fastq
rm $FQGEO/${SAMPLE}_R2.fastq
fi


cd $ORIG
for FILE in ${SAMPLE}*_R1_001.fastq.gz
   do gunzip -c ${FILE} >> $FQGEO/${SAMPLE}_R1.fastq
done

for FILE in ${SAMPLE}*_R2_001.fastq.gz
   do gunzip -c ${FILE} >> $FQGEO/${SAMPLE}_R2.fastq
done


NUMPROCS=2
OUTPUT=$PROJ/ercc_aln_se
mkdir -p $OUTPUT




LOGS=$OUTPUT/job_logs
ALN_ERROR_LOG=$LOGS/size_errors.txt
if [ ! -d $LOGS ]; then mkdir -p $LOGS; fi

MYLOG=$LOGS/${SAMPLE}.log
printf "`date` $SAMPLE PBS_ARRAYID:$PBS_ARRAYID HOSTNAME:$HOSTNAME\n" > $MYLOG



ALNLOGS=$OUTPUT/aln_logs
if [ ! -d $ALNLOGS ]; then mkdir -p $ALNLOGS; fi

UNALNDIR=$OUTPUT/unaln_reads
if [ ! -d $UNALNDIR ]; then mkdir -p $UNALNDIR; fi



source /blue/concannon/share/jnewman/qc_pipeline/scripts/alignment_functions.sh

READS=$FQGEO/${SAMPLE}_R1.fastq
NAME=${SAMPLE}_R1

echo $SAMPLE
echo $READS
echo $NAME 

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


READS=$FQGEO/${SAMPLE}_R2.fastq
NAME=${SAMPLE}_R2
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


