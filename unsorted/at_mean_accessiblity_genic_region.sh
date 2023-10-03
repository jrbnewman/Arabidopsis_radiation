#!/bin/bash
PATH=$HOME/conda/.envs/lab/bin:$PATH

PROJ=$HOME/concannon/DTRA

IN=at_mean_accessiblity_across_gene

python $PROJ/scripts/lineplot_nosd.py -i $PROJ/${IN}.csv -min 0 -max 0.5 -o $PROJ/${IN}_lineplot.pdf -g dose -y accessibility -x pos -t "Mean accessibility across genic regions"

