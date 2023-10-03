#!/bin/bash

## Create GC accessibility plots using new annotations (transcripts from GTF, and HOMER)

PROJ=$PATCON/arabidopsis_wgbs_cold
PYSCRIPT=$PROJ/scripts/tss_accessibility_plots_FANS_03jrbn.py

##

for CONDIT in 22C_100U 22C_100U_0U 4C_100U_0U 22C_100U_0U_common 4C_100U_0U_common FANS_0p5U FANS_1p5U FANS_25U FANS_5U; do
    echo "Making plots for " ${CONDIT}
    ## TSS, all sites, no bin GTF annotation
    INPUT1=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_nobin_transcript_TSS.csv
    OUTPUT1=$PROJ/analysis_output/GC_accessibility_1kb_TSS_all_sites_nobin_transcript_TSS_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT1} --condition ${CONDIT} --output ${OUTPUT1}

    ## TSS, 10X sites, no bin GTF annotation
    INPUT2=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_nobin_transcript_TSS.csv
    OUTPUT2=$PROJ/analysis_output/GC_accessibility_1kb_TSS_10X_sites_nobin_transcript_TSS_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT2} --condition ${CONDIT} --output ${OUTPUT2}

    ## TSS, all sites, binned GTF annotation
    INPUT3=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_binned_transcript_TSS.csv
    OUTPUT3=$PROJ/analysis_output/GC_accessibility_1kb_TSS_all_sites_binned_transcript_TSS_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT3} --condition ${CONDIT} --output ${OUTPUT3}

    ## TSS, 10X sites, binned GTF annotation
    INPUT4=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_binned_transcript_TSS.csv
    OUTPUT4=$PROJ/analysis_output/GC_accessibility_1kb_TSS_10X_sites_binned_transcript_TSS_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT4} --condition ${CONDIT} --output ${OUTPUT4}

    ## TSS, all sites, no bin HOMER annotation
    INPUT5=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_nobin_HOMER.csv
    OUTPUT5=$PROJ/analysis_output/GC_accessibility_1kb_TSS_all_sites_nobin_HOMER_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT5} --condition ${CONDIT} --output ${OUTPUT5}

    ## TSS, 10X sites, no bin HOMER annotation
    INPUT6=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_nobin_HOMER.csv
    OUTPUT6=$PROJ/analysis_output/GC_accessibility_1kb_TSS_10X_sites_nobin_HOMER_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT6} --condition ${CONDIT} --output ${OUTPUT6}

    ## TSS, all sites, binned HOMER annotation
    INPUT7=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_binned_HOMER.csv
    OUTPUT7=$PROJ/analysis_output/GC_accessibility_1kb_TSS_all_sites_binned_HOMER_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT7} --condition ${CONDIT} --output ${OUTPUT7}

    ## TSS, 10X sites, binned HOMER annotation
    INPUT8=$PROJ/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_binned_HOMER.csv
    OUTPUT8=$PROJ/analysis_output/GC_accessibility_1kb_TSS_10X_sites_binned_HOMER_${CONDIT}_20200212.pdf
    python ${PYSCRIPT} --input ${INPUT8} --condition ${CONDIT} --output ${OUTPUT8}
    echo "Done for " ${CONDIT} "!"
done
