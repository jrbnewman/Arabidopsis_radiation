ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";

/* Calculating the means for treatment*time. */

/* Get coverage counts. Going to calculate means on the non-log, UQ-normalized APN */

data coverage;
   set rs.arab_fusion_counts_q3_norm_cent;
   keep sample_number fusion_id q3_q3_apn;
run;

/* Get timepoints and treatment for samples */
   
data sample_key;
  set rs.arab_design_file;
  keep sample_number treatment time;
run;

proc sort data=sample_key nodup;
   by sample_number;
proc sort data=coverage;
   by sample_number;
run;

data fusions_counts_w_key;
   merge sample_key (in=in1) coverage (in=in2);
   by sample_number;
   if in1 and in2;
   run;

/* Calc means for treatment*time by fusion */

proc sort data=fusions_counts_w_key;
   by fusion_id treatment time;
proc means data=fusions_counts_w_key noprint;
   by fusion_id treatment time;
   var q3_q3_apn;
   output out=mean_apn_by_trt_time mean=;
run;


/* Transpose */

data mean_apn_by_trt_time_2;
  set mean_apn_by_trt_time;
  trt_time=catx('_',treatment,time);
run;

proc sort data=mean_apn_by_trt_time_2;
  by fusion_id trt_time;
proc transpose data=mean_apn_by_trt_time_2 out=mean_apn_sbys;
   by fusion_id;
   var q3_q3_apn;
   id trt_time;
run;

/* Make permenant */

data rs.arab_fus_mean_apn_by_trt_time;
   set mean_apn_sbys;
   drop _NAME_;
   rename _0_1gy_1=mean_q3apn_01gy_1h
          _0_1gy_3=mean_q3apn_01gy_3h
          _0_1gy_24=mean_q3apn_01gy_24h
          _0_1gy_72=mean_q3apn_01gy_72h
          _1gy_1=mean_q3apn_1gy_1h
          _1gy_3=mean_q3apn_1gy_3h
          _1gy_24=mean_q3apn_1gy_24h
          _1gy_72=mean_q3apn_1gy_72h
          Mock_1=mean_q3apn_mock_1h
          Mock_3=mean_q3apn_mock_3h
          Mock_24=mean_q3apn_mock_24h
          Mock_72=mean_q3apn_mock_72h ;
run;


