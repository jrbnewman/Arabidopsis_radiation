#!/bin/bash
PATH=$HOME/conda/.envs/lab/bin:$PATH

PROJ=$HOME/concannon/DTRA

IN=at_mean_accessiblity_per_gene

python $PROJ/scripts/boxplot_violinplot.py -i $PROJ/${IN}.csv -o $PROJ/${IN}_boxplot -g dose -y accessibility -x dose -t "mean genic accessibility by exposure"
