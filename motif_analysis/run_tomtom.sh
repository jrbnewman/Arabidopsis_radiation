PROJ=/mnt/c/Users/Yamanto/Desktop/dtra_arabidopsis_manuscript/at_rad_motif_analysis
MOTIF1=/mnt/e/meme/motif_databases/ARABD/ArabidopsisPBM_20140210.meme
MOTIF2=/mnt/e/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme


for IN1 in dmr_01_72_cg dmr_01_72_chg dmr_01_72_chh dmr_1_72_cg dmr_1_72_chg dmr_1_72_chh dn_dar_01_72 dn_dar_1_72 dn_dmr_01_72_cg dn_dmr_01_72_chg dn_dmr_01_72_chh dn_dmr_1_72_cg dn_dmr_1_72_chg dn_dmr_1_72_chh up_dar_01_72 up_dar_1_72 up_dmr_01_72_cg up_dmr_01_72_chg up_dmr_01_72_chh up_dmr_1_72_cg up_dmr_1_72_chg up_dmr_1_72_chh
do
	echo ${IN1}

	FASTAIN=$PROJ/dreme_sequences/${IN1}.fasta.masked
	STREMEOUT=$PROJ/streme_analysis/${IN1}/streme.txt

	TOMOUT1=$PROJ/tomtom_output_PBM/${IN1}
	mkdir -p $TOMOUT1

	TOMOUT2=$PROJ/tomtom_output_DAP/${IN1}
	mkdir -p $TOMOUT2

echo "${IN1} Tomtom"
ls $MOTIF1
ls $MOTIF2

tomtom --oc $TOMOUT1 $STREMEOUT $MOTIF1
tomtom --oc $TOMOUT2 $STREMEOUT $MOTIF2

done

