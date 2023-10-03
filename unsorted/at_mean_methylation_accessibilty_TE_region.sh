#!/bin/bash
PATH=$HOME/conda/.envs/lab/bin:$PATH

PROJ=$HOME/concannon/DTRA

IN1=at_mean_accessiblity_across_TEs
IN2=at_mean_CG_methylation_across_TEs
IN3=at_mean_CHG_methylation_across_TEs
IN4=at_mean_CHH_methylation_across_TEs

python $PROJ/scripts/lineplot_nosd.py -i $PROJ/${IN1}.csv -min 0 -max 0.2 -o $PROJ/${IN1}_lineplot.pdf -g dose -y accessibility -x pos -t "Mean GC accessibility across genic regions"

python $PROJ/scripts/lineplot_nosd.py -i $PROJ/${IN2}.csv -min 0 -max 0.2 -o $PROJ/${IN2}_lineplot.pdf -g dose -y methylation -x pos -t "Mean CG methylation across genic regions"

python $PROJ/scripts/lineplot_nosd.py -i $PROJ/${IN3}.csv -min 0 -max 0.2 -o $PROJ/${IN3}_lineplot.pdf -g dose -y methylation -x pos -t "Mean CHG methylation across genic regions"

python $PROJ/scripts/lineplot_nosd.py -i $PROJ/${IN4}.csv -min 0 -max 0.2 -o $PROJ/${IN4}_lineplot.pdf -g dose -y methylation -x pos -t "Mean CHH methylation across genic regions"


