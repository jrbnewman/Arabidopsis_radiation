#!/bin/bash
PATH=$HOME/conda/.envs/lab/bin:$PATH

PROJ=$HOME/concannon/DTRA

for TYPE in CHG CHH
do
IN=at_mean_${TYPE}_methylation_across_gene
python $PROJ/scripts/lineplot_nosd.py -i $PROJ/${IN}.csv -min 0 -max 0.2 -o $PROJ/${IN}_lineplot.pdf -g dose -y methylation -x pos -t "Mean ${TYPE} methylation across genic regions"
done
