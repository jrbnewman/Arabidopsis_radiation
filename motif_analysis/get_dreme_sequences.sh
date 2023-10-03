#!/bin/bash 

## extract FASTA sequence from BED files for motif analysis 


PROJ=$HOME/concannon/DTRA/at_rad_motif_analysis

INDIR=$PROJ/input_data
OUTDIR=$PROJ/dreme_sequences
mkdir -p $OUTDIR

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

FASTA=$HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

bedtools getfasta -bed ${INDIR}/${IN1}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN1}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN2}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN2}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN3}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN3}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN4}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN4}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN5}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN5}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN6}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN6}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN7}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN7}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN8}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN8}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN9}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN9}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN10}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN10}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN11}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN11}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN12}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN12}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN13}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN13}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN14}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN14}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN15}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN15}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN16}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN16}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN17}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN17}.fasta -name
bedtools getfasta -bed ${INDIR}/${IN18}.bed -fi ${FASTA} -fo ${OUTDIR}/${IN18}.fasta -name


