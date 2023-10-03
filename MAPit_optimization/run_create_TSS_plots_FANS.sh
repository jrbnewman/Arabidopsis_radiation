#!/bin/sh

python tss_accessibility_plots_FANS_02jrbn.py --input ../analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_binned_02jrbn.csv  --output ../analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_binned_plot_02jrbn
python tss_accessibility_plots_FANS_02jrbn.py --input ../analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_nobin_02jrbn.csv --output ../analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_nobin_plot_02jrbn
python tss_accessibility_plots_FANS_02jrbn.py --input ../analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_binned_02jrbn.csv --output ../analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_binned_plot_02jrbn
python tss_accessibility_plots_FANS_02jrbn.py --input ../analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_nobin_02jrbn.csv --output ../analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_nobin_plot_02jrbn


