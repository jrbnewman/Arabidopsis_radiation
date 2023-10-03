#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --job-name=star_refs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/out.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#36


export PATH=/blue/concannon/share/jnewman/software/TEtranscripts_py36/bin:$PATH

PROJ=/blue/concannon/share/jnewman/mingqi_arab

GENEGTF=$PROJ/references/Arabidopsis_thaliana.TAIR10.37.gtf
TEGTF=$PROJ/references/TAIR10_TE.gtf


OUTPUT=$PROJ/TEtranscripts_output
mkdir -p $OUTPUT

TRT1=$PROJ/star_aln/7-1-1-1_Aligned.out.bam
TRT2=$PROJ/star_aln/8-1-1-2_Aligned.out.bam
TRT3=$PROJ/star_aln/9-1-1-3_Aligned.out.bam
CTL1=$PROJ/star_aln/1-M-1-1_Aligned.out.bam
CTL2=$PROJ/star_aln/2-M-1-2_Aligned.out.bam
CTL3=$PROJ/star_aln/3-M-1-3_Aligned.out.bam
CONTRAST=100_v_0cgy_1h


TEtranscripts -t ${TRT1} ${TRT2} ${TRT3} \
              -c ${CTL1} ${CTL2} ${CTL3} \
              --GTF $GENEGTF \
              --TE $TEGTF \
              --project $OUTPUT/${CONTRAST}_uniq \
              --mode uniq \
              --minread 1 \
              --padj 0.05 \
              --foldchange 1

TEtranscripts -t ${TRT1} ${TRT2} ${TRT3} \
              -c ${CTL1} ${CTL2} ${CTL3} \
              --GTF $GENEGTF \
              --TE $TEGTF \
              --project $OUTPUT/${CONTRAST}_multi \
              --mode multi \
              --minread 1 \
              --padj 0.05 \
              --foldchange 1

