/* Running model -- by fusion, time: genotype*time*treatment, plate as random */

ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";



/* Get list of "on" fusions */
* as we are now doing contrasts that involve ethanol-only or control-only comparisons, I am changing;
* this from "on in both control and ETOH" to "on in either control and ETOH";
* Can always drop fusions/events later, but before FDR;


data on_fusions;
   set rs.arab_flag_fusion_on_gt0;
   if flag_fusion_on_mock_apn0=1 or flag_fusion_on_01gy_apn0=1 or flag_fusion_on_1gy_apn0=1 ;
   keep fusion_id;
run;

/* Merge in with count data */

proc sort data=on_fusions;
   by fusion_id;
proc sort data=rs.arab_fusion_counts_q3_norm_cent;
   by fusion_id;
run;

data fusion_counts_on;
   merge rs.arab_fusion_counts_q3_norm_cent (in=in1) on_fusions (in=in2);
   by fusion_id;
   if in1 and in2;
run;

/* Merge in design file
   Want variables: treatment, time */
   
data sample_key;
  set rs.arab_design_file;
  keep sample_number treatment time;
run;

proc sort data=sample_key nodup;
   by sample_number;
proc sort data=fusion_counts_on;
   by sample_number;
run;

data fusions_on_w_key;
   merge sample_key (in=in1) fusion_counts_on (in=in2);
   by sample_number;
   if in1 and in2;
   run;

/* Run model */

* Need to sort by fusion and time!;

proc sort data=fusions_on_w_key;
   by fusion_id time;
   run;

ods listing close;

%macro runModels(measure,outname);
   
ods listing close;
proc glimmix data=fusions_on_w_key;
  by fusion_id;
  class time treatment ;
  model &measure. = time|treatment / htype=1;
  output out=resid resid=resid pred=pred student=stu;
ods output tests1=anova;
run;
quit;

/* Flag residuals */

proc univariate data = resid normal noprint;
  by fusion_id;
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

/* Make permenant */

data rs.arab_resid_main_trt_tm_&outname.;
  set flag_resids;
  run;
data rs.arab_anova_main_trt_tm_&outname. ;
  set anova ;
  run ;    
data rs.arab_norm_main_trt_tm_&outname. ;
  set fusions_flag_fail_norm ;
  run ;

%mend;

%runModels(log_q3_q3_apn,non);
*%runModels(mean_log_q3_center,mean);
*%runModels(median_log_q3_center,med);
*%runModels(q3_log_q3_center,q3);
   


