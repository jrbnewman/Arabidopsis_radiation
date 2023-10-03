#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=motif
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/dataPrep_%A.%a.out
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G


module load meme/5.4.1

PROJ=/blue/concannon/share/jnewman/mingqi_arab

INDIR=$PROJ/dreme_sequences
MOTIFDB=$PROJ/references/motif_databases/ARABD/ArabidopsisPBM_20140210.meme

OUTDIR1=$PROJ/streme_analysis
OUTDIR2=$PROJ/sea_analysis

IN1=dn_dar_01_72
IN2=dn_dar_1_72
IN3=dn_dar_72_shared2
IN4=dn_dmr_01_72_merged2
IN5=dn_dmr_1_72_merged2
IN6=dn_dmr_72_CG_shared2
IN7=dn_dmr_72_CHG_shared2
IN8=dn_dmr_72_CHH_shared2
IN9=dn_dmr_72_shared2
IN10=up_dar_01_72
IN11=up_dar_1_72
IN12=up_dar_72_shared2
IN13=up_dmr_01_72_merged2
IN14=up_dmr_1_72_merged2
IN15=up_dmr_72_CG_shared2
IN16=up_dmr_72_CHG_shared2
IN17=up_dmr_72_CHH_shared2
IN18=up_dmr_72_shared2


mkdir -p ${OUTDIR1}/${IN1}
mkdir -p ${OUTDIR2}/${IN1}

streme --verbosity 1 --oc ${OUTDIR1}/${IN1} --dna --totallength 4000000 --time 14400 --minw 8 --maxw 15 --thresh 0.05 --align center --desc ${IN1} --p ${INDIR}/${IN1}.fasta
sea --verbosity 1 --oc ${OUTDIR2}/${IN1} --thresh 10.0 --align center --p ${INDIR}/${IN1}.fasta --m ${MOTIFDB}


