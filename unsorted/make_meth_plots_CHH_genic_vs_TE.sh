#!/bin/bash


PATH=$HOME/conda/.envs/lab/bin:$PATH
PROJ=$PATCON/DTRA


INDIR=$PROJ
OUTDIR=$PROJ
mkdir -p $OUTDIR


for FILE in dmr_plot_chh_up_dmr_01_72_100_bins dmr_plot_chh_up_dmr_1_72_100_bins dmr_plot_chh_dn_dmr_01_72_100_bins dmr_plot_chh_dn_dmr_1_72_100_bins
do
python $PROJ/scripts/methylation_plot.py \
--input $INDIR/${FILE}.csv \
--output $OUTDIR/${FILE}.pdf \
--y-min 0 \
--y-max 100 \
--title "CHH methylation at DMRs (center +/- 1kb)" \
--y-label "Methylation (%)" \
--x-label "Position relative to center of DMR"

done

