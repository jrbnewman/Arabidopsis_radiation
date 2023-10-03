#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=checkFQ
#SBATCH --account=concannon
#SBATCH --qos=brusko-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/checkFQ_%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --array=1-36


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


gzip $FQGEO/${SAMPLE}_R1.fastq
gzip $FQGEO/${SAMPLE}_R2.fastq
 
 
 
export PATH=/home/jrbnewman/.conda/envs/fqutils/bin:$PATH

echo "Checking new R1 - RNA"
fastq_info ${FQGEO}/${SAMPLE}_R1.fastq.gz

echo "Checking new R2 - RNA"
fastq_info ${FQGEO}/${SAMPLE}_R2.fastq.gz

echo "Checking new R1+R2 - RNA"
fastq_info ${FQGEO}/${SAMPLE}_R1.fastq.gz ${FQGEO}/${SAMPLE}_R2.fastq.gz

