#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=motif
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/dataPrep_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G


module load meme/5.4.1

PROJ=/blue/concannon/share/jnewman/mingqi_arab

INDIR=$PROJ/dreme_sequences
MOTIFDB=$PROJ/references/motif_databases/ARABD/ArabidopsisPBM_20140210.meme

OUTDIR1=$PROJ/streme_analysis
OUTDIR2=$PROJ/sea_analysis

for IN1 in At_DEG_de_01_1_promoter At_DEG_de_01_24_promoter At_DEG_de_01_3_promoter At_DEG_de_01_72_promoter At_DEG_de_1_1_promoter At_DEG_de_1_24_promoter At_DEG_de_1_3_promoter At_DEG_de_1_72_promoter At_DEG_de_1_promoter At_DEG_de_24_promoter At_DEG_de_3_promoter At_DEG_de_72_promoter At_DEG_de_any_promoter At_DEG_dn_01_1_promoter At_DEG_dn_01_24_promoter At_DEG_dn_01_3_promoter At_DEG_dn_01_72_promoter At_DEG_dn_1_1_promoter At_DEG_dn_1_24_promoter At_DEG_dn_1_3_promoter At_DEG_dn_1_72_promoter At_DEG_up_01_1_promoter At_DEG_up_01_24_promoter At_DEG_up_01_3_promoter At_DEG_up_01_72_promoter At_DEG_up_1_1_promoter At_DEG_up_1_24_promoter At_DEG_up_1_3_promoter At_DEG_up_1_72_promoter
do

mkdir -p ${OUTDIR1}/${IN1}
mkdir -p ${OUTDIR2}/${IN1}

streme --verbosity 1 --oc ${OUTDIR1}/${IN1} --dna --totallength 4000000 --time 14400 --minw 8 --maxw 15 --thresh 0.05 --align center --desc ${IN1} --p ${INDIR}/${IN1}.fasta
sea --verbosity 1 --oc ${OUTDIR2}/${IN1} --thresh 10.0 --align center --p ${INDIR}/${IN1}.fasta --m ${MOTIFDB}

done
