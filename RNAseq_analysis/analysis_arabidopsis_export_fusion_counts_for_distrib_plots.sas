ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Prep data for exporting so I can plot the distribution of coverage. We can use this to decide if
   there are low-coverage samples that need dropping

   I want to include: sampleID, fusionID, log_apn, plate, condition (trt*time) */

data sample_key;
  set rs.arab_design_file;
  length condition $20.;
  condition=catx("_",treatment,time);
  keep sample_number condition treatment time;
run;


proc sort data=sample_key nodup;
  by sample_number;
proc sort data=rs.arab_fusion_counts_w_key;
  by sample_number;
run;

data counts_w_key;
   merge sample_key (in=in1) rs.arab_fusion_counts_w_key (in=in2);
   by sample_number;
   if in1 and in2;
run;

data calc_log;
  set counts_w_key;
  log_apn=log(apn+1);
run;

proc sort data=calc_log;
  by sample_number fusion_id;
run;


data drop_zero;
  set calc_log;
  where apn > 0;
run;

/* Export for plotting coverage distribution */

proc export data=calc_log
     outfile="!PATCON/arabidopsis/analysis_output/distrib_logapn_fusions_arabidopsis.csv"
     dbms=csv replace;
run;

proc export data=drop_zero
     outfile="!PATCON/arabidopsis/analysis_output/distrib_logapn_fusions_arabidopsis_gt0.csv"
     dbms=csv replace;
run;

