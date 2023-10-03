#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=csfilter
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOGS/out.build_cs_ref.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb


### Load modules

module load dibig_tools 

### Set directories

PROJ=/ufrc/concannon/share/jnewman/mingqi_arab
FASTA=/ufrc/data/reference/icbr/TAIR10/TAIR10_fa/TAIR10.fna

cscall build -r $FASTA -s CG -o $PROJ/references/TAIR10-CG.bin 
cscall build -r $FASTA -s GC -o $PROJ/references/TAIR10-GC.bin 
cscall build -r $FASTA -s CHG -o $PROJ/references/TAIR10-CHG.bin 
cscall build -r $FASTA -s CHH -o $PROJ/references/TAIR10-CHH.bin 
