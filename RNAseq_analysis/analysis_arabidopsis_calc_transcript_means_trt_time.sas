ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";

/* Calculating the means for treatment*time. */

/* Get coverage counts. Going to calculate means on the non-log, UQ-normalized APN */

data coverage;
   set rs.arab_xscript_counts_q3_norm_cent;
   keep sample_id transcript_id q3_q3_tpm;
run;

/* Get timepoints and treatment for samples */
   
data sample_key;
  set rs.arab_design_file;
   length sample_id2 $50.;
   if grays=0 then sample_id2=catx("-",sample_number,"M",time,replicate);
   else if grays=0.1 then sample_id2=catx("-",sample_number,"0-1",time,replicate);
   else if grays=1 then sample_id2=catx("-",sample_number,"1",time,replicate);
  keep sample_id2 treatment time;
  rename sample_id2=sample_id;
run;

proc sort data=sample_key nodup;
   by sample_id;
proc sort data=coverage;
   by sample_id;
run;

data xscript_counts_w_key;
   merge sample_key (in=in1) coverage (in=in2);
   by sample_id;
   if in1 and in2;
   run;

/* Calc means for treatment*time by fusion */

proc sort data=xscript_counts_w_key;
   by transcript_id treatment time;
proc means data=xscript_counts_w_key noprint;
   by transcript_id treatment time;
   var q3_q3_tpm;
   output out=mean_tpm_by_trt_time mean=;
run;


/* Transpose */

data mean_tpm_by_trt_time_2;
  set mean_tpm_by_trt_time;
  trt_time=catx('_',treatment,time);
run;

proc sort data=mean_tpm_by_trt_time_2;
  by transcript_id trt_time;
proc transpose data=mean_tpm_by_trt_time_2 out=mean_tpm_sbys;
   by transcript_id;
   var q3_q3_tpm;
   id trt_time;
run;

/* Make permenant */

data rs.arab_xs_mean_tpm_by_trt_time;
   set mean_tpm_sbys;
   drop _NAME_;
   rename _0_1gy_1=mean_q3tpm_01gy_1h
          _0_1gy_3=mean_q3tpm_01gy_3h
          _0_1gy_24=mean_q3tpm_01gy_24h
          _0_1gy_72=mean_q3tpm_01gy_72h
          _1gy_1=mean_q3tpm_1gy_1h
          _1gy_3=mean_q3tpm_1gy_3h
          _1gy_24=mean_q3tpm_1gy_24h
          _1gy_72=mean_q3tpm_1gy_72h
          Mock_1=mean_q3tpm_mock_1h
          Mock_3=mean_q3tpm_mock_3h
          Mock_24=mean_q3tpm_mock_24h
          Mock_72=mean_q3tpm_mock_72h ;
run;


