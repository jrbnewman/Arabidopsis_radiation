#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --job-name=star_refs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/out.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --array=1-36

#36


#module load R python/3.6

#export PYTHONPATH=/blue/concannon/share/jnewman/software/TEtranscripts/lib/python3.8/site-packages:$PYTHONPATH
export PATH=/blue/concannon/share/jnewman/software/TEtranscripts_py36/bin:$PATH

PROJ=/blue/concannon/share/jnewman/mingqi_arab

GENEGTF=$PROJ/references/Arabidopsis_thaliana.TAIR10.37.gtf
TEGTF=$PROJ/references/TAIR10_TE.gtf


DESIGN_FILE=$PROJ/design_files/arabidopsis_samples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}

BAM=$PROJ/star_aln/${SAMPLEID}_Aligned.out.bam

OUTPUT1=$PROJ/TEcount_output_uniq/$SAMPLEID
OUTPUT2=$PROJ/TEcount_output_multi/$SAMPLEID

mkdir -p $OUTPUT2
mkdir -p $OUTPUT1



# Run TEcount to generate counts table

TEcount -b $BAM --GTF $GENEGTF --TE $TEGTF --project $OUTPUT1 --mode uniq
TEcount -b $BAM --GTF $GENEGTF --TE $TEGTF --project $OUTPUT2 --mode multi


