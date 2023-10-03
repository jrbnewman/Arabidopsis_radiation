#!/bin/bash

PROJ=/nfshome/jrbnewman/concannon/DTRA
PATH=$HOME/conda/.envs/lab/bin:$PATH

#for GENE in CMT3 CMT1 DRM3 MET2 CMT2 DRM2 MET1 new_ROS1_AT2G36490;
#do

GENE=new_ROS1_AT2G36490_gene
INPUT=$PROJ/at_rad_exp_${GENE}.csv
OUTPNG1=$PROJ/at_rad_exp_${GENE}_CPM.pdf
OUTPNG2=$PROJ/at_rad_exp_${GENE}_logCPM.pdf
python $PROJ/scripts/lineplot.py -i $INPUT -o $OUTPNG1 -x timepoint -y cpm -g group -t "${GENE}"
#python $PROJ/scripts/lineplot.py -i $INPUT -o $OUTPNG2 -x timepoint -y log_cpm -g group -t "${GENE}"

INPUT=$PROJ/at_rad_exp_${GENE}.csv
OUTPNG1=$PROJ/at_rad_exp_${GENE}_CPM_nosd.pdf
OUTPNG2=$PROJ/at_rad_exp_${GENE}_logCPM_nosd.pdf
#python $PROJ/scripts/lineplot_nosd.py -i $INPUT -o $OUTPNG1 -x timepoint -y cpm -g group -t "${GENE}"
#python $PROJ/scripts/lineplot_nosd.py -i $INPUT -o $OUTPNG2 -x timepoint -y log_cpm -g group -t "${GENE}"

#done
