#!/bin/bash


PATH=$HOME/conda/.envs/lab/bin:$PATH
PROJ=$PATCON/DTRA


INDIR=$PROJ
OUTDIR=$PROJ
mkdir -p $OUTDIR


for FILE in dar_plot_up_dar_01_72_100_bins dar_plot_up_dar_1_72_100_bins dar_plot_dn_dar_01_72_100_bins dar_plot_dn_dar_1_72_100_bins
do
python $PROJ/scripts/methylation_plot.py \
--input $INDIR/${FILE}.csv \
--output $OUTDIR/${FILE}.pdf \
--y-min 0 \
--y-max 100 \
--title "GC accessibility at DARs (center +/- 1kb)" \
--y-label "Methylation 100U-0U (%)" \
--x-label "Position relative to center of DAR"

done



