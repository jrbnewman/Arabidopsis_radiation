#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=metilene_1gy
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/metilene_1gy_0gy.out
#SBATCH -t 1:00:00
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

# CG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CG.output.txt
#metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX} -f 3 -d 0.1 ${INPUT} > ${OUTPUT}

# CHG sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHG.output.txt
#metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX} -f 3 -d 0.1 ${INPUT} > ${OUTPUT}

# CHH sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_CHH.output.txt
#metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX} -f 3 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC.output.txt
#metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX} -f 3 -d 0.1 ${INPUT} > ${OUTPUT}

# GC sites (no negative)
INPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_noneg.txt
OUTPUT=$PROJ/at_rad_metilene_0gy_1gy_GC_no_negative.output.txt
metilene -a ${GROUP1_PREFIX} -b ${GROUP2_PREFIX} -f 3 -d 0.1 ${INPUT} > ${OUTPUT}

