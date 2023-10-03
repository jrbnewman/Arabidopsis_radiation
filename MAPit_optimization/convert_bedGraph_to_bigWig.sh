#!/bin/sh

# Convert bedGraphs into bigWigs for Mike

INDIR=$PATCON/arabidopsis_wgbs_cold/analysis_output/bedGraphs
OUTDIR=${INDIR}/bigWigs

if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

BG2BW=$PATCON/ucsc_tools/bedGraphToBigWig

CHROMSIZE=$PATCON/useful_arabidopsis_data/TAIR10/downloaded_files/chrom.sizes

cd $INDIR

for FILE in *.bedGraph
do
   BASENAME=$(basename $FILE .bedGraph)
   echo $BASENAME
${BG2BW} ${FILE} ${CHROMSIZE} ${OUTDIR}/${BASENAME}.bw
done

