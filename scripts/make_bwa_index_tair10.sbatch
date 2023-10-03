#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=bwa_index
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS.out.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb
#

module load bwa/0.7.15

REF=/blue/concannon/share/jnewman/references/tair10

PREFIX=$REF/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
FASTA=$REF/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

bwa index -a bwtsw -p $PREFIX $FASTA

