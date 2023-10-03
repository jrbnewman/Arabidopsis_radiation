#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=freebayes
#SBATCH --mail-type=FAIL
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.freebayes.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --array=1-7


###1-25



mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
    module load samtools vcflib freebayes/1.3.1 bcftools bamaddrg vcftools


# Set paths
    PROJ=/blue/concannon/share/jnewman/mingqi_arab
    OUTPUT=$PROJ/bam_out

    FASTA=/blue/concannon/share/jnewman/references/tair10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

BCFDIR=$PROJ/vcf_files_freebayes
mkdir -p $BCFDIR


if [ ${SLURM_ARRAY_TASK_ID} == 7 ]; then
   CHROM=Pt
elif [ ${SLURM_ARRAY_TASK_ID} == 6 ]; then
  CHROM=Mt
else
  CHROM=${SLURM_ARRAY_TASK_ID}
fi


#bamaddrg -R ${CHROM} -b ${OUTPUT}/1-M-1-1.sorted.bam -s 0cGy_1h_rep1 -r 0cGy_1h_rep1 -b ${OUTPUT}/2-M-1-2.sorted.bam -s 0cGy_1h_rep2 -r 0cGy_1h_rep2 -b ${OUTPUT}/3-M-1-3.sorted.bam -s 0cGy_1h_rep3 -r 0cGy_1h_rep3 -b ${OUTPUT}/19-M-24-1.sorted.bam -s 0cGy_24h_rep1 -r 0cGy_24h_rep1 -b ${OUTPUT}/20-M-24-2.sorted.bam -s 0cGy_24h_rep2 -r 0cGy_24h_rep2 -b ${OUTPUT}/21-M-24-3.sorted.bam -s 0cGy_24h_rep3 -r 0cGy_24h_rep3 -b ${OUTPUT}/10-M-3-1.sorted.bam -s 0cGy_3h_rep1 -r 0cGy_3h_rep1 -b ${OUTPUT}/11-M-3-2.sorted.bam -s 0cGy_3h_rep2 -r 0cGy_3h_rep2 -b ${OUTPUT}/12-M-3-3.sorted.bam -s 0cGy_3h_rep3 -r 0cGy_3h_rep3 -b ${OUTPUT}/28-M-72-1.sorted.bam -s 0cGy_72h_rep1 -r 0cGy_72h_rep1 -b ${OUTPUT}/29-M-72-2.sorted.bam -s 0cGy_72h_rep2 -r 0cGy_72h_rep2 -b ${OUTPUT}/30-M-72-3.sorted.bam -s 0cGy_72h_rep3 -r 0cGy_72h_rep3 -b ${OUTPUT}/7-1-1-1.sorted.bam -s 100cGy_1h_rep1 -r 100cGy_1h_rep1 -b ${OUTPUT}/8-1-1-2.sorted.bam -s 100cGy_1h_rep2 -r 100cGy_1h_rep2 -b ${OUTPUT}/9-1-1-3.sorted.bam -s 100cGy_1h_rep3 -r 100cGy_1h_rep3 -b ${OUTPUT}/25-1-24-1.sorted.bam -s 100cGy_24h_rep1 -r 100cGy_24h_rep1 -b ${OUTPUT}/26-1-24-2.sorted.bam -s 100cGy_24h_rep2 -r 100cGy_24h_rep2 -b ${OUTPUT}/27-1-24-3.sorted.bam -s 100cGy_24h_rep3 -r 100cGy_24h_rep3 -b ${OUTPUT}/16-1-3-1.sorted.bam -s 100cGy_3h_rep1 -r 100cGy_3h_rep1 -b ${OUTPUT}/17-1-3-2.sorted.bam -s 100cGy_3h_rep2 -r 100cGy_3h_rep2 -b ${OUTPUT}/18-1-3-3.sorted.bam -s 100cGy_3h_rep3 -r 100cGy_3h_rep3 -b ${OUTPUT}/34-1-72-1.sorted.bam -s 100cGy_72h_rep1 -r 100cGy_72h_rep1 -b ${OUTPUT}/35-1-72-2.sorted.bam -s 100cGy_72h_rep2 -r 100cGy_72h_rep2 -b ${OUTPUT}/36-1-72-3.sorted.bam -s 100cGy_72h_rep3 -r 100cGy_72h_rep3 -b ${OUTPUT}/4-0-1-1-1.sorted.bam -s 10cGy_1h_rep1 -r 10cGy_1h_rep1 -b ${OUTPUT}/5-0-1-1-2.sorted.bam -s 10cGy_1h_rep2 -r 10cGy_1h_rep2 -b ${OUTPUT}/6-0-1-1-3.sorted.bam -s 10cGy_1h_rep3 -r 10cGy_1h_rep3 -b ${OUTPUT}/22-0-1-24-1.sorted.bam -s 10cGy_24h_rep1 -r 10cGy_24h_rep1 -b ${OUTPUT}/23-0-1-24-2.sorted.bam -s 10cGy_24h_rep2 -r 10cGy_24h_rep2 -b ${OUTPUT}/24-0-1-24-3.sorted.bam -s 10cGy_24h_rep3 -r 10cGy_24h_rep3 -b ${OUTPUT}/13-0-1-3-1.sorted.bam -s 10cGy_3h_rep1 -r 10cGy_3h_rep1 -b ${OUTPUT}/14-0-1-3-2.sorted.bam -s 10cGy_3h_rep2 -r 10cGy_3h_rep2 -b ${OUTPUT}/15-0-1-3-3.sorted.bam -s 10cGy_3h_rep3 -r 10cGy_3h_rep3 -b ${OUTPUT}/31-0-1-72-1.sorted.bam -s 10cGy_72h_rep1 -r 10cGy_72h_rep1 -b ${OUTPUT}/32-0-1-72-2.sorted.bam -s 10cGy_72h_rep2 -r 10cGy_72h_rep2 -b ${OUTPUT}/33-0-1-72-3.sorted.bam -s 10cGy_72h_rep3 -r 10cGy_72h_rep3 | freebayes-v1.3.1 --stdin -f $FASTA -p 2 > ${BCFDIR}/all_samples_arabidopsis_RNAseq.freeBayes.called.unfiltered.${CHROM}.vcf

#all_2021_clones.tetra.freeBayes.called.unfiltered.${CHROM}.vcf

vcfutils.pl varFilter -Q 20 -d 10  ${BCFDIR}/all_samples_arabidopsis_RNAseq.freeBayes.called.unfiltered.${CHROM}.vcf > ${BCFDIR}/all_samples_arabidopsis_RNAseq.freeBayes.called.minQ20_D10.${CHROM}.vcf


