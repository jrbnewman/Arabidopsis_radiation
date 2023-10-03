#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=rsem_ref
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/mingqi_arab/scripts/SLURM_LOG_out.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem=16gb


PROJ=/ufrc/concannon/share/jnewman/mingqi_arab
REF=$PROJ/references

module load R/3.2.0
module load rsem/1.2.28
module load star/2.5.1b
module load cufflinks/2.2.1
module load bowtie/1.1.1

rsem-prepare-reference --gtf ${REF}/Arabidopsis_thaliana.TAIR10.37.gtf --star --star-sjdboverhang 100 ${REF}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ${REF}/TAIR10.rsem.tg

perl -ne 'if(/^#/){}else{if(/gene_id \"(.*)\"\; transcript_id "(.*)"\; gene_name "(.*)"; gene_source.*"\;/){print $1,"\t",$3,"\n"}}' ${REF}/Arabidopsis_thaliana.TAIR10.37.gtf | sort -u | sort -k1,1 -k2,2 > ${REF}/tair10.gid.gnm.txt
perl -ne 'if(/^#/){}else{if(/gene_id \"(.*)\"\; transcript_id "(.*)"\; exon_number .*\; gene_name "(.*)"; gene_source.*"\;/){print $1,"\t",$2,"\t",$3,"\n"}}' ${REF}/Arabidopsis_thaliana.TAIR10.37.gtf | sort -u | sort -k1,1 -k2,2 > ${REF}/tair10.gid.tid.gnm.txt

