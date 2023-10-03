#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=mpileup
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/make_mpileup_%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
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

OUTBAM=$PROJ/bam_out
OUTWIG=$OUTPUT/mpileups_genome_tair10_uniq

mkdir -p $OUTBAM
mkdir -p $OUTWIG

# Create LOG directory and start log
LOGS=$OUTBAM/logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
MYLOG=$OUTBAM/${SAMPLEID}_mpileip.log
#printf "`date` $SAMPLEID  HOSTNAME:$HOSTNAME \n\n" > $MYLOG


#### Merge sam files

#ls $SPLITSAM/${SAMPLEID}*.sam
#cat ${SPLITSAM}/${SAMPLEID}-*.sam > $OUTBAM/${SAMPLEID}.sam
#rm ${OUTPUT}/${SAMPLEID}-*.sam

BAM=$OUTBAM/${SAMPLEID}
PILEUP=$OUTWIG/${SAMPLEID}.mpileup

#printf "<-------------------- Convert SAM to BAM -------------------->" >> $MYLOG
#echo `date`": Starting SAM to BAM conversion" >> $MYLOG
#samtools view -ut $REF.fai -o $BAM.bam $OUTBAM/${SAMPLEID}.sam 2>> $MYLOG
#samtools sort -T ${BAM}.tmp.sorted -o $BAM.sorted.bam $BAM.bam 2>> "${MYLOG}"
#samtools index $BAM.sorted.bam >> $MYLOG
#echo `date`": Finished SAM to BAM conversion" >> $MYLOG

#### Make mpielup

#printf "<-------------------- Convert BAM to MPILEUP -------------------->" >> $MYLOG
#echo `date`": Generating pileup" >> $MYLOG
#samtools mpileup -d 1000000000 -f $REF $BAM.sorted.bam  > $PILEUP 2>> $MYLOG

#rm ${BAM}.bam

OUTCC=$PROJ/TE_elements_uniq
mkdir -p $OUTCC

#### COVERAGE COUNTS
    BED=/blue/concannon/share/jnewman/mingqi_arab/references/TE_fragments.bed
echo "************** FUSION COUNTS ****************"

    echo "Starting Coverage Counts for $NAME (fusions) `date`" 
    python $PROJ/scripts/rpkm_calculate.py \
        -b $BED \
        -m ${PILEUP} \
        -s $OUTBAM/${SAMPLEID}.sam \
        -n ${SAMPLEID} \
        --cv \
        -o $OUTCC/cvrg_cnts_${SAMPLEID}_uniq.csv
    echo "Finished Coverage Counts for $NAME (fusions) `date`"






### all alignments
OUTBAM=$PROJ/bam_out2
OUTWIG=$OUTPUT/mpileups_genome_tair10_all

mkdir -p $OUTBAM
mkdir -p $OUTWIG

# Create LOG directory and start log
LOGS=$OUTBAM/logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
#MYLOG=$OUTBAM/${SAMPLEID}_mpileip.log
printf "`date` $SAMPLEID  HOSTNAME:$HOSTNAME \n\n" > $MYLOG


#### Merge sam files

#ls $SPLITSAM/${SAMPLEID}*.sam
#cat ${OUTPUT}/${SAMPLEID}-*.sam > $OUTBAM/${SAMPLEID}.sam
#rm ${OUTPUT}/${SAMPLEID}-*.sam

BAM=$OUTBAM/${SAMPLEID}
PILEUP=$OUTWIG/${SAMPLEID}.mpileup

#printf "<-------------------- Convert SAM to BAM -------------------->" >> $MYLOG
#echo `date`": Starting SAM to BAM conversion" >> $MYLOG
#samtools view -ut $REF.fai -o $BAM.bam $OUTBAM/${SAMPLEID}.sam 2>> $MYLOG
#samtools sort -T ${BAM}.tmp.sorted -o $BAM.sorted.bam $BAM.bam 2>> "${MYLOG}"
#samtools index $BAM.sorted.bam >> $MYLOG
#echo `date`": Finished SAM to BAM conversion" >> $MYLOG

#### Make mpielup

#printf "<-------------------- Convert BAM to MPILEUP -------------------->" >> $MYLOG
#echo `date`": Generating pileup" >> $MYLOG
#samtools mpileup -d 1000000000 -f $REF $BAM.sorted.bam  > $PILEUP 2>> $MYLOG

#rm ${BAM}.*


OUTCC=$PROJ/TE_elements_all
mkdir -p $OUTCC
#### COVERAGE COUNTS
    BED=/blue/concannon/share/jnewman/mingqi_arab/references/TE_fragments.bed
echo "************** FUSION COUNTS ****************"

    echo "Starting Coverage Counts for $NAME (fusions) `date`" 
    python $PROJ/scripts/rpkm_calculate.py \
        -b $BED \
        -m ${PILEUP} \
        -s $OUTBAM/${SAMPLEID}.sam \
        -n ${SAMPLEID} \
        --cv \
        -o $OUTCC/cvrg_cnts_${SAMPLEID}_all.csv
    echo "Finished Coverage Counts for $NAME (fusions) `date`"


