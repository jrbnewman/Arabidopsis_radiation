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
	
	rm $OUTPUT/${NAME}.sam
	rm $SPLITSAM/${NAME}_gray.sam
	rm $SPLITSAM/${NAME}_ambiguous.sam
	rm $SPLITSAM/${NAME}_unrecognized.sam
        rm $SPLITSAM/*.fq	
done

## Merge all SAM files, convert to BAM then make a wiggle

OUTBAM=$PROJ/bam_out
OUTWIG=$OUTPUT/mpileups_genome_tair10

mkdir -p $OUTBAM
mkdir -p $OUTWIG

# Create LOG directory and start log
LOGS=$OUTBAM/logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
MYLOG=$OUTBAM/${SAMPLEID}_mpileip.log
printf "`date` $SAMPLEID  HOSTNAME:$HOSTNAME \n\n" > $MYLOG


#### Merge sam files 

ls $SPLITSAM/${SAMPLEID}*.sam
cat ${SPLITSAM}/${SAMPLEID}-*.sam > $OUTBAM/${SAMPLEID}.sam
rm ${OUTPUT}/${SAMPLEID}-*.sam

BAM=$OUTBAM/${SAMPLEID}
PILEUP=$OUTWIG/${SAMPLEID}.mpileup

printf "<-------------------- Convert SAM to BAM -------------------->" >> $MYLOG
echo `date`": Starting SAM to BAM conversion" >> $MYLOG
samtools view -ut $REF.fai -o $BAM.bam $OUTBAM/${SAMPLEID}.sam 2>> $MYLOG
samtools sort -T ${BAM}.tmp.sorted -o $BAM.sorted.bam $BAM.bam 2>> "${MYLOG}"
samtools index $BAM.sorted.bam >> $MYLOG
echo `date`": Finished SAM to BAM conversion" >> $MYLOG


#### Make mpielup

printf "<-------------------- Convert BAM to MPILEUP -------------------->" >> $MYLOG
echo `date`": Generating pileup" >> $MYLOG
samtools mpileup -d 1000000000 -f $REF $BAM.sorted.bam  > $PILEUP 2>> $MYLOG

rm ${BAM}.*


