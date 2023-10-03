#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=mpileup
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/make_mpileup_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --array=1-36


###36
#


date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load python/2.7.6 bwa/0.7.15 samtools/1.3.1

### 64 samples total

### FQsplitdups on raw reads
#Set directories

PROJ=/blue/concannon/share/jnewman/mingqi_arab
FASTQ=$PROJ/unzipped/fastqSplitDups


OUTPUT=$PROJ/bwa_mem_aln_pe
SPLITSAM=$OUTPUT/bwa_parse

mkdir -p $OUTPUT
mkdir -p $SPLITSAM

## Design file
DESIGN_FILE=$PROJ/design_files/arabidopsis_samples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}

## set reference
REF=/blue/concannon/share/jnewman/references/tair10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

cd ${FASTQ}
for FILE in ${SAMPLEID}-*_R1_001_distinct.fq
	do
	READ1=${FILE/.fq/}
	READ2=${READ1/_R1_/_R2_}
	NAME=${READ1/_R1_001_distinct/}
	
	echo $READ1
	echo $READ2
	
	#### parse bwa alignment file
	
	#### run bwa-mem
	bwa mem -t 2 -M $REF $FASTQ/${READ1}.fq $FASTQ/${READ2}.fq > $OUTPUT/${NAME}.sam
	
	echo "************* PARSE BWA OUTPUT ************"
	
	python /blue/concannon/share/IKZF1_clones/scripts/BWASplitSAM_07mai.py -fq1 $FASTQ/${READ1}.fq -fq2 $FASTQ/${READ2}.fq -s $OUTPUT/${NAME}.sam --outdir $SPLITSAM
	
	## clean up
	
	#rm $OUTPUT/${NAME}.sam
        rm $SPLITSAM/*.fq	
done



