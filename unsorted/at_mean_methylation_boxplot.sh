#!/bin/bash
PATH=$HOME/conda/.envs/lab/bin:$PATH

PROJ=$HOME/concannon/DTRA

for TYPE in CG CHG CHH
do
IN=at_mean_${TYPE}_methylation_per_gene

python $PROJ/scripts/boxplot_violinplot.py -i $PROJ/${IN}.csv -o $PROJ/${IN}_boxplot -g dose -y methylation -x dose -t "mean genic ${TYPE} methylation by exposure"
done

