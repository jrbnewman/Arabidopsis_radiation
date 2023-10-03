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

module load star

PROJ=/blue/concannon/share/jnewman/mingqi_arab
NUMPROC=4

REFDIR=$PROJ/references/TAIR_STAR
if [ ! -e $REFDIR ]; then mkdir -p $REFDIR; fi

FASTA=$REFDIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
GTF=$REFDIR/Arabidopsis_thaliana.TAIR10.37.gtf
JUNCLENGTH=100


STAR --runThreadN ${NUMPROC} \
     --runMode genomeGenerate \
     --genomeDir ${REFDIR} \
     --genomeFastaFiles ${FASTA} \
     --sjdbGTFfile ${GTF} \
     --sjdbOverhang ${JUNCLENGTH}


