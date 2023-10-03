#!/bin/bash
#SBATCH --job-name=makeRMdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --cpus-per-task=1
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mem=12gb
#SBATCH -t 12:00:00
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/repeatModeler_make_db.out


module load repeatmodeler

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export ROZ=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}


PROJ=/blue/concannon/share/jnewman/mingqi_arab
FASTA=/blue/concannon/share/jnewman/references/Brapa_15/Brapa_R500_V1.2.fasta
DBOUT=$PROJ/Brapa_R500

cd $PROJ

BuildDatabase -name $DBOUT $FASTA


