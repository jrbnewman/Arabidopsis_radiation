#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=fqDups
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/dataPrep_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-36


###36
#


date;hostname


module load python/2.7.14


### 64 samples total

### FQsplitdups on raw reads
#Set directories

PROJ=/blue/concannon/share/jnewman/mingqi_arab
ORIG=$PROJ/original_data

FQ=$PROJ/unzipped
    mkdir -p $FQ

LOGS=$PROJ/dataPrepLogs_PE
    mkdir -p $LOGS
FQSPLIT=$PROJ/unzipped/fastqSplitDups
    mkdir -p $FQSPLIT




## Design file
DESIGN_FILE=$PROJ/design_files/arabidopsis_samples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}


cd ${ORIG}
for FILE in ${SAMPLEID}-*_R1_001.fastq.gz
	do
	READ1=${FILE/.fastq.gz/}
	READ2=${READ1/_R1_/_R2_}
	NAME=${READ1/_R1_001/}

	echo $READ1
	echo $READ2

	gunzip -c ${READ1}.fastq.gz > $FQ/${READ1}.fastq
	gunzip -c ${READ2}.fastq.gz > $FQ/${READ2}.fastq

	python /blue/concannon/share/IKZF1_clones/scripts/fastqSplitDups_2MAI.py \
		-r1 $FQ/${READ1}.fastq \
		-r2 $FQ/${READ2}.fastq \
		--outdir $FQSPLIT \
		-o $LOGS/${NAME}_fqSplitDup_raw_summary.csv \
		-t $FQSPLIT/${NAME}_fqSplitDup_raw_table.tsv \
		-g $LOGS/${NAME}_fqSplitDup_raw.log

done
