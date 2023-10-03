#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=cscall
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.cscall.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --array=1-12

### Load modules

module samtools dibig_tools

PROJ=/ufrc/concannon/share/jnewman/mingqi_arab


# Design file
    DESIGN_FILE=$PROJ/design_files/sample_design_file_nofq.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

        GROUP=${ARRAY[0]}
        SAMPLE=${ARRAY[1]}

INPUT=$PROJ/BSmap_output/merged_bams

bisconv.py split ${INPUT}/${SAMPLE}.duprg.bam ${INPUT}/${SAMPLE}.C.bam ${INPUT}/${SAMPLE}.G.bam

for SITE in CG CHG CHH GC; do
METHREF=/ufrc/data/reference/icbr/TAIR10/TAIR10_fa/TAIR10-${SITE}.bin
REF=/ufrc/data/reference/icbr/TAIR10/TAIR10_fa/TAIR10.fna

OUTBED=$PROJ/BSmap_output/cscall_output
   if [ ! -e $OUTBED ]; then mkdir -p $OUTBED; fi
REPORTS=$OUTBED/reports
   if [ ! -e $REPORTS ]; then mkdir -p $REPORTS; fi

cscall -call -ad \
  -i ${METHREF} \
  -r $REF \
  -o $OUTBED/${SAMPLE}_${SITE}.bed \
  -report $REPORTS/${SAMPLE}-report.csv \
  -lc $REPORTS/${SAMPLE}-lc.csv \
  -hist $REPORTS/${SAMPLE}-hist.csv \
  -mat $REPORTS/${SAMPLE}-mat.csv  \
  -mind 1 \
  -mins 1 \
  -b ${INPUT}/${SAMPLE}.C.bam \
  -bt ${INPUT}/${SAMPLE}.G.bam 
done
