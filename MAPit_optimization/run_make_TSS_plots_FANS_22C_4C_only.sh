#!/bin/bash

PROJ=$PATCON/DTRA/arabidopsis_wgbs_cold

## Make TSS plots

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/FANS_only_GC_1kb_TSS_10X_sites_nobin_HOMER.csv \
	--output $PROJ/analysis_output/FANS_only_GC_TSS_plot_10X_nobin_202009240.pdf \
	--y-min 0 \
	--y-max 70 \
	--title "GC methylation - FANS" \
	--y-label "GC methylation (%)" \
	--x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
        --input $PROJ/analysis_output/FANS_only_GC_1kb_TSS_10X_sites_binned_HOMER.csv \
        --output $PROJ/analysis_output/FANS_only_GC_TSS_plot_10X_binned_202009240.pdf \
        --y-min 0 \
	--y-max 70 \
	--title "GC methylation - FANS" \
        --y-label "GC methylation (%)" \
        --x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
        --input $PROJ/analysis_output/GC_22C_4C_100U_0U_1kb_TSS_10X_sites_nobin_HOMER.csv \
        --output $PROJ/analysis_output/GC_22C_4C_100U_0U_TSS_plot_10X_nobin_202009240.pdf \
        --y-min 0 \
        --y-max 70 \
        --title "GC methylation - 22C and 4C" \
        --y-label "GC methylation (%)" \
        --x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
        --input $PROJ/analysis_output/GC_22C_4C_100U_0U_1kb_TSS_10X_sites_binned_HOMER.csv \
        --output $PROJ/analysis_output/GC_22C_4C_100U_0U_TSS_plot_10X_binned_202009240.pdf \
        --y-min 0 \
        --y-max 70 \
        --title "GC methylation - 22C and 4C" \
        --y-label "GC methylation (%)" \
        --x-label "Position relative to TSS"


