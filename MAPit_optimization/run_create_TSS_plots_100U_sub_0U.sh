#!/bin/sh
python tss_accessibility_plots_02jrbn.py --input ../analysis_output/GC_accessibility_1kb_TSS_10X_sites_binned.csv  --output ../analysis_output/GC_accessibility_1kb_TSS_10X_sites_binned_plot_zoomin
python tss_accessibility_plots_02jrbn.py --input ../analysis_output/GC_accessibility_1kb_TSS_10X_sites_nobin.csv --output ../analysis_output/GC_accessibility_1kb_TSS_10X_sites_nobin_plot_zoomin
python tss_accessibility_plots_02jrbn.py --input ../analysis_output/GC_accessibility_1kb_TSS_all_sites_binned.csv --output ../analysis_output/GC_accessibility_1kb_TSS_all_sites_binned_plot_zoomin
python tss_accessibility_plots_02jrbn.py --input ../analysis_output/GC_accessibility_1kb_TSS_all_sites_nobin.csv --output ../analysis_output/GC_accessibility_1kb_TSS_all_sites_nobin_plot_zoomin



