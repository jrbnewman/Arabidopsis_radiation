#!/bin/bash

PROJ=$PATCON/DTRA/arabidopsis_wgbs_cold

## Make TSS plots

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/FANS_GC_1kb_TSS_10X_sites_nobin_HOMER.csv \
	--output $PROJ/analysis_output/FANS_GC_TSS_plot_10X_nobin_20200530.pdf \
	--y-min 0.1 \
	--y-max 0.7 \
	--title "GC methylation - FANS vs 22C 100U" \
	--y-label "GC methylation (%)" \
	--x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/FANS_GC_1kb_TSS_10X_sites_binned_HOMER.csv \
	--output $PROJ/analysis_output/FANS_GC_TSS_plot_10X_binned_20200530.pdf \
	--y-min 0.1 \
	--y-max 0.7 \
	--title "GC methylation - FANS vs 22C 100U" \
	--y-label "GC methylation (%)" \
	--x-label "Position relative to TSS"

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/GC_diff_1kb_TSS_10X_sites_nobin_HOMER.csv \
	--output $PROJ/analysis_output/GC_diff_TSS_plot_10X_nobin_20200530.pdf \
	--y-min 0.3 \
	--y-max 0.6 \
	--title "Accessibility - 22C vs 4C" \
	--y-label "GC methylation (100U - 0U) (%)" \
	--x-label "Position relative to TSS"

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/GC_diff_1kb_TSS_10X_sites_binned_HOMER.csv \
	--output $PROJ/analysis_output/GC_diff_TSS_plot_10X_binned_20200530.pdf \
	--y-min 0.3 \
	--y-max 0.6 \
	--title "Accessibility - 22C vs 4C" \
	--y-label "GC methylation (100U - 0U) (%)" \
	--x-label "Position relative to TSS"

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/CG_diff_1kb_TSS_10X_sites_nobin_HOMER.csv \
	--output $PROJ/analysis_output/CG_diff_TSS_plot_10X_nobin_20200530.pdf \
	--y-min 0 \
	--y-max 0.35 \
	--title "CG methylation - 22C vs 4C" \
	--y-label "Methylation (%)" \
	--x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/CG_diff_1kb_TSS_10X_sites_binned_HOMER.csv \
	--output $PROJ/analysis_output/CG_diff_TSS_plot_10X_binned_20200530.pdf \
	--y-min 0 \
	--y-max 0.35 \
	--title "CG methylation - 22C vs 4C" \
	--y-label "Methylation (%)" \
	--x-label "Position relative to TSS"

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/CHG_diff_1kb_TSS_10X_sites_nobin_HOMER.csv \
	--output $PROJ/analysis_output/CHG_diff_TSS_plot_10X_nobin_20200530.pdf \
	--y-min 0 \
	--y-max 0.15 \
	--title "CHG methylation - 22C vs 4C" \
	--y-label "Methylation (%)" \
	--x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/CHG_diff_1kb_TSS_10X_sites_binned_HOMER.csv \
	--output $PROJ/analysis_output/CHG_diff_TSS_plot_10X_binned_20200530.pdf \
	--y-min 0 \
	--y-max 0.15 \
	--title "CHG methylation - 22C vs 4C" \
	--y-label "Methylation (%)" \
	--x-label "Position relative to TSS"

python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/CHH_diff_1kb_TSS_10X_sites_nobin_HOMER.csv \
	--output $PROJ/analysis_output/CHH_diff_TSS_plot_10X_nobin_20200530.pdf \
	--y-min 0 \
	--y-max 0.05 \
	--title "CHH methylation - 22C vs 4C" \
	--y-label "Methylation (%)" \
	--x-label "Position relative to TSS"


python $PROJ/scripts/methylation_plot.py \
	--input $PROJ/analysis_output/CHH_diff_1kb_TSS_10X_sites_binned_HOMER.csv \
	--output $PROJ/analysis_output/CHH_diff_TSS_plot_10X_binned_20200530.pdf \
	--y-min 0 \
	--y-max 0.05 \
	--title "CHH methylation - 22C vs 4C" \
	--y-label "Methylation (%)" \
	--x-label "Position relative to TSS"


