PROJ=/mnt/c/Users/Yamanto/Desktop/dtra_arabidopsis_manuscript/at_rad_motif_analysis

INDIR=$PROJ/dreme_sequences
MOTIFDB=/mnt/e/meme/motif_databases/ARABD/ArabidopsisPBM_20140210.meme

OUTDIR1=$PROJ/streme_analysis
#OUTDIR2=$PROJ/sea_analysis

#for IN1 in dn_dar_01_72 dn_dar_1_72 dn_dar_72_shared2 dn_dmr_01_72_merged2 dn_dmr_1_72_merged2 dn_dmr_72_CG_shared2 dn_dmr_72_CHG_shared2 dn_dmr_72_CHH_shared2 dn_dmr_72_shared2 up_dar_01_72 up_dar_1_72 up_dar_72_shared2 up_dmr_01_72_merged2 up_dmr_1_72_merged2 up_dmr_72_CG_shared2 up_dmr_72_CHG_shared2 up_dmr_72_CHH_shared2 up_dmr_72_shared2;

#for IN1 in dn_dar_01_72 dn_dar_1_72 dn_dar_72_shared2 dn_dmr_72_CG_shared2 dn_dmr_72_CHG_shared2 dn_dmr_72_CHH_shared2 dn_dmr_72_shared2 up_dar_01_72 up_dar_1_72 up_dar_72_shared2 up_dmr_72_CG_shared2 up_dmr_72_CHG_shared2 up_dmr_72_CHH_shared2 up_dmr_72_shared2; do

#for IN1 in dn_dmr_01_72_cg dn_dmr_01_72_chg dn_dmr_01_72_chh dn_dmr_1_72_cg dn_dmr_1_72_chg dn_dmr_1_72_chh up_dmr_01_72_cg up_dmr_01_72_chg up_dmr_01_72_chh up_dmr_1_72_cg up_dmr_1_72_chg up_dmr_1_72_chh dmr_01_72_cg dmr_1_72_cg dmr_01_72_chg dmr_1_72_chg dmr_01_72_chh dmr_1_72_chh

for IN1 in dmr_01_72_cg dmr_01_72_chg dmr_01_72_chh dmr_1_72_cg dmr_1_72_chg dmr_1_72_chh dn_dar_01_72 dn_dar_1_72 dn_dmr_01_72_cg dn_dmr_01_72_chg dn_dmr_01_72_chh dn_dmr_1_72_cg dn_dmr_1_72_chg dn_dmr_1_72_chh up_dar_01_72 up_dar_1_72 up_dmr_01_72_cg up_dmr_01_72_chg up_dmr_01_72_chh up_dmr_1_72_cg up_dmr_1_72_chg up_dmr_1_72_chh 
do
	echo ${IN1}

mkdir -p ${OUTDIR1}/${IN1}
#mkdir -p ${OUTDIR2}/${IN1}

streme --verbosity 1 --oc ${OUTDIR1}/${IN1} --dna --totallength 4000000 --time 14400 --minw 8 --maxw 15 --thresh 0.05 --align center --desc ${IN1} --p ${INDIR}/${IN1}.fasta.masked
#sea --verbosity 1 --oc ${OUTDIR2}/${IN1} --thresh 10.0 --align center --p ${INDIR}/${IN1}.fasta --m ${MOTIFDB}

done


