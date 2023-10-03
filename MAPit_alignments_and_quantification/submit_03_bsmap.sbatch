#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=BSMap
#SBATCH --account=concannon
##SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.BSMap.%j.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --array=1-36

## max array 24

#### Whole genome bisulfite sequencing analysis
#### Step 3) Mapping
####
#### In this script, the following will be performed:
####    (1) Map reads using one of three mappers (BSMap, BS-Seeker2, BWA-Meth)
####    (2) Convert SAM to BAM if necessary
####    (3) Use `samtools flagstat` to count SAM flags from each output
####

module load dibig_tools samtools bamtools
source /apps/dibig_tools/1.0/lib/sh/utils.sh

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

### Inputs

TRIMMEDFQ=$PROJ/qc_trim_output

TRIMR1=good-${NAME}_R1.trimmed.paired.fastq.gz
TRIMR2=good-${NAME}_R2.trimmed.paired.fastq.gz

### Output

OUTPUT=$PROJ/BSmap_output
   if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

### (1) Run BSMap

BSMAPPROG=/apps/dibig_tools/dmap/bsmap
REF=/blue/data/reference/icbr/TAIR10/TAIR10_fa/TAIR10.fna

$$BSMAPPROG -d $REF -a $TRIMMEDFQ/$TRIMR1 -b $TRIMMEDFQ/$TRIMR2 -o $OUTPUT/${NAME}.bam -p 4



### (2) SAMtools flagstat

FLAGSTATS=$OUTPUT/samtools_stats_by_chr
 if  [ ! -e $FLAGSTATS ]; then mkdir -p $FLAGSTATS; fi

samtools view $OUTPUT/${NAME}.bam 
samtools flagstat $OUTPUT/${NAME}.bam > $FLAGSTATS/${NAME}_samflag_stats.txt




