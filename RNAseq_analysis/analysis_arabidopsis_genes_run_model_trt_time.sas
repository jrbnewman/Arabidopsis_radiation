/* Running model -- by fusion, time: treatment, plate as random */

ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";


/* Get list of "on" fusions */

data on_gene;
   set rs.arab_flag_gene_on_cpm_gt0;
   if flag_gene_on_mock_cpm0=1 and flag_gene_on_01gy_cpm0=1 and flag_gene_on_1gy_cpm0=1 ;
   keep gene_id;
run;

/* Merge in with count data */

proc sort data=on_gene;
   by gene_id;
proc sort data=rs.cpm_norm_counts_by_gene;
   by gene_id;
run;

data gene_counts_on;
   merge rs.cpm_norm_counts_by_gene (in=in1) on_gene (in=in2);
   by gene_id;
   if in1 and in2;
   log_cpm=log(cpm + 1);
run;

/* Merge in design file
   Want variables: treatment, genotype, time, plate */
   
/* Run model */

* Need to sort by fusion and time!;

proc sort data=gene_counts_on;
   by gene_id time;
   run;


%macro runModels(measure,outname);
ods listing close;
proc glimmix data=gene_counts_on;
  by gene_id time;
  class treatment;
  model &measure. = treatment / ddfm=kr;
  output out=resid resid=resid pred=pred student=stu;
ods output tests3=anova;
run;
quit;


/* Flag residuals */

proc univariate data = resid normal noprint;
  by gene_id time;
  var Resid;
  output out = normtest probn=pnorm;
  run;

data flag_resids;
  set normtest;
  if pnorm = . then flag_fail_norm = .;
        else if pnorm le 0.05 then flag_fail_norm = 1;
        else flag_fail_norm = 0;
  run;

proc freq data = flag_resids noprint;
  tables flag_fail_norm / out=fusions_flag_fail_norm;
  run;

ods listing;
proc print data=fusions_flag_fail_norm;
run;

/* Make permenant */

data rs.arab_gene_resids_trt_&outname.;
  set flag_resids;
  run;
data rs.arab_gene_anova_gt_trt_&outname.;
  set anova ;
  run ;
data rs.arab_gene_norm_gt_trt_&outname.;
  set fusions_flag_fail_norm ;
  run ;

%mend;

%runModels(log_cpm,lcpm);
%runModels(cpm,cpm);


