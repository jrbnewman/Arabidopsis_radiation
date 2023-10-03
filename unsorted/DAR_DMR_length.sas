libname arabMAP "/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs";

/* Calc distribution of DAR length */

ods listing;
ods html close;
proc datasets lib=work kill noprint;
run;
quit;

data dar_length;
  set arabMAP.results_by_dar_annot;
  dar_length = region_stop - region_start;
  if mean_methyl_diff_TRT_CTL > 0 then dar_direction=1;
  else if mean_methyl_diff_TRT_CTL < 0 then dar_direction=-1;
  else dar_direction=0;
run;

proc sort data=dar_length;
  by comparison dar_direction;
proc means data=dar_length noprint;
  by comparison dar_direction;
  var dar_length;
  output out=dar_length_distrib mean=mean stddev=stddev min=min p5=p5 p10=p10 q1=q1 median=medain q3=q3 p90=p90 p95=p95 max=max;
run;

proc print data=dar_length_distrib;
run;

data dar_length_gt100;
  set dar_length;
  if dar_length >= 100;
run;

proc sort data=dar_length_gt100;
  by comparison dar_direction;
proc means data=dar_length_gt100 noprint;
  by comparison dar_direction;
  var dar_length;
  output out=dar_length_gt100_distrib mean=mean stddev=stddev min=min p5=p5 p10=p10 q1=q1 median=medain q3=q3 p90=p90 p95=p95  max=max;
run;

proc print data=dar_length_gt100_distrib;
run;


data dmr_length;
  set arabMAP.results_by_dmr_annot;
  if (num_sites_FDR05_diff_10perc>=2 and abs(mean_methyl_diff) > 0) or (flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1);
  dmr_length = region_stop - region_start;
  if mean_methyl_diff > 0 then dmr_direction=1;
  else if mean_methyl_diff < 0 then dmr_direction=-1;
  else dmr_direction=0;
run;

proc sort data=dmr_length;
  by site_type comparison dmr_direction;
proc means data=dmr_length noprint;
  by site_type comparison dmr_direction;
  var dmr_length;
  output out=dmr_length_distrib mean=mean stddev=stddev min=min p5=p5 p10=p10 q1=q1 median=medain q3=q3 p90=p90 p95=p95  max=max;
run;

proc print data=dmr_length_distrib;
run;



