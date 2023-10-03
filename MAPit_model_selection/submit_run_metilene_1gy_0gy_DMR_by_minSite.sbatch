#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=metilene_1gy
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/metilene_1gy_0gy_%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=12G



date;hostname


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# Load modules
module load metilene

PROJ=/blue/concannon/share/jnewman/mingqi_arab

GROUP1_PREFIX=C1gy
GROUP2_PREFIX=C0gy
## 2 sites

# CG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CG_minSites_2.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CG_100cGy_minSites_2.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# CHG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG_minSites_2.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHG_100cGy_minSites_2.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}


# CHH sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH_minSites_2.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHH_100cGy_minSites_2.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_minSites_2.DAR.output.txt

metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_2.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites (no negative)
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_noneg.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_no_negative_minSites_2.DAR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_2.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

## 3 sites

# CG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CG_minSites_3.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CG_100cGy_minSites_3.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# CHG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG_minSites_3.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHG_100cGy_minSites_3.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}


# CHH sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH_minSites_3.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHH_100cGy_minSites_3.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_minSites_3.DAR.output.txt

metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_3.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites (no negative)
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_noneg.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_no_negative_minSites_3.DAR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_3.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

## 5 sites

# CG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CG_minSites_5.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CG_100cGy_minSites_5.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# CHG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG_minSites_5.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHG_100cGy_minSites_5.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}


# CHH sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH_minSites_5.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHH_100cGy_minSites_5.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_minSites_5.DAR.output.txt

metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_5.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites (no negative)
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_noneg.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_no_negative_minSites_5.DAR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_5.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

## 10 sites

# CG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CG_minSites_10.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CG_100cGy_minSites_10.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# CHG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG_minSites_10.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHG_100cGy_minSites_10.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}


# CHH sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH_minSites_10.DMR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DMR_CHH_100cGy_minSites_10.bed -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_minSites_10.DAR.output.txt

metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_10.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites (no negative)
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_noneg.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_no_negative_minSites_10.DAR.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX}  -B $PROJ/DAR_GC_100cGy_minSites_10.bed  -f 2 -d 0.1 ${INPUT} > ${OUTPUT}


