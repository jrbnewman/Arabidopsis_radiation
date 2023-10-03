ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Flag fusions on/off at APN>0 and APN>5 */

data counts_w_key;
  set rs.counts_by_isoform;
  if tpm > 0 then flag_tpm_gt0=1;
  else flag_tpm_gt0=0;

  if tpm ge 5 then flag_tpm_ge5=1;
  else flag_tpm_ge5=0;
run;

/* Flag on per treatment */

proc sort data=counts_w_key;
   by treatment transcript_id;
proc means data=counts_w_key noprint;
   by treatment transcript_id;
   var flag_tpm_gt0 flag_tpm_ge5;
   output out=mean_on_by_trt mean=;
run;

data flag_iso_on_by_trt;
  set mean_on_by_trt;

  /* For TPM > 0 */
  if flag_tpm_gt0 ge 0.5 then flag_xscript_on_gt0=1;
  else if flag_tpm_gt0 = 0 then flag_xscript_on_gt0=0;
  else flag_xscript_on_gt0=.;

  /* For TPM ge 5 */
  if flag_tpm_ge5 ge 0.5 then flag_xscript_on_ge5=1;
  else if flag_tpm_ge5 = 0 then flag_xscript_on_ge5=0;
  else flag_xscript_on_ge5=.;

  keep transcript_id treatment flag_xscript_on_gt0 flag_xscript_on_ge5;
run;

proc sort data=flag_iso_on_by_trt;
   by transcript_id treatment;
proc transpose data=flag_iso_on_by_trt out=flag_iso_trt_sbys_tpm0;
   by transcript_id;
   id treatment;
   var flag_xscript_on_gt0;
run;

proc transpose data=flag_iso_on_by_trt out=flag_iso_trt_sbys_tpm5;
   by transcript_id;
   id treatment;
   var flag_xscript_on_ge5;
run;

data treat_tpm0;
   set flag_iso_trt_sbys_tpm0;
   drop _NAME_;
   rename _0_1gy=flag_xscript_on_01gy_apn0
          _1gy=flag_xscript_on_1gy_apn0
          Mock=flag_xscript_on_mock_apn0;
run;

data treat_tpm5;
   set flag_iso_trt_sbys_tpm5;
   drop _NAME_;
   rename _0_1gy=flag_xscript_on_01gy_apn5
          _1gy=flag_xscript_on_1gy_apn5
          Mock=flag_xscript_on_mock_apn5;
run;

/* Flag on per time */

proc sort data=counts_w_key;
   by time transcript_id;
proc means data=counts_w_key noprint;
   by time transcript_id;
   var flag_tpm_gt0 flag_tpm_ge5;
   output out=mean_on_by_time mean=;
run;

data flag_iso_on_by_time;
  set mean_on_by_time;

  /* For APN > 0 */
  if flag_tpm_gt0 ge 0.5 then flag_xscript_on_gt0=1;
  else if flag_tpm_gt0 = 0 then flag_xscript_on_gt0=0;
  else flag_xscript_on_gt0=.;

  /* For APN ge 5 */
  if flag_tpm_ge5 ge 0.5 then flag_xscript_on_ge5=1;
  else if flag_tpm_ge5 = 0 then flag_xscript_on_ge5=0;
  else flag_xscript_on_ge5=.;

  keep transcript_id time flag_xscript_on_gt0 flag_xscript_on_ge5;
run;

proc sort data=flag_iso_on_by_time;
   by transcript_id time;
proc transpose data=flag_iso_on_by_time out=flag_iso_time_sbys_tpm0;
   by transcript_id;
   id time;
   var flag_xscript_on_gt0;
run;

proc transpose data=flag_iso_on_by_time out=flag_iso_time_sbys_tpm5;
   by transcript_id;
   id time;
   var flag_xscript_on_ge5;
run;

data time_tpm0;
   set flag_iso_time_sbys_tpm0;
   drop _NAME_;
   rename _1=flag_xscript_on_1hr_tpm0
          _3=flag_xscript_on_3hr_tpm0
          _24=flag_xscript_on_24hr_tpm0
          _72=flag_xscript_on_72hr_tpm0;
run;

data time_tpm5;
   set flag_iso_time_sbys_tpm5;
   drop _NAME_;
   rename _1=flag_xscript_on_1hr_tpm5
          _3=flag_xscript_on_3hr_tpm5
          _24=flag_xscript_on_24hr_tpm5
          _72=flag_xscript_on_72hr_tpm5;
run;

/* Flag on per treatment-by-time */

proc sort data=counts_w_key;
   by treatment time transcript_id;
proc means data=counts_w_key noprint;
   by treatment time transcript_id;
   var flag_tpm_gt0 flag_tpm_ge5;
   output out=mean_on_by_trttime mean=;
run;

data flag_iso_on_by_trttime;
  length condition $10.;
  set mean_on_by_trttime;
  condition=catx("_",treatment,time);

  /* For TPM > 0 */
  if flag_tpm_gt0 ge 0.5 then flag_xscript_on_gt0=1;
  else if flag_tpm_gt0 = 0 then flag_xscript_on_gt0=0;
  else flag_xscript_on_gt0=.;

  /* For TPM ge 5 */
  if flag_tpm_ge5 ge 0.5 then flag_xscript_on_ge5=1;
  else if flag_tpm_ge5 = 0 then flag_xscript_on_ge5=0;
  else flag_xscript_on_ge5=.;

  keep transcript_id condition flag_xscript_on_gt0 flag_xscript_on_ge5;
run;

proc sort data=flag_iso_on_by_trttime;
   by transcript_id condition;
proc transpose data=flag_iso_on_by_trttime out=flag_iso_trttime_sbys_tpm0;
   by transcript_id;
   id condition;
   var flag_xscript_on_gt0;
run;

proc transpose data=flag_iso_on_by_trttime out=flag_iso_trttime_sbys_tpm5;
   by transcript_id;
   id condition;
   var flag_xscript_on_ge5;
run;

data trttime_tpm0;
   set flag_iso_trttime_sbys_tpm0;
   drop _NAME_;
   rename _0_1gy_1=flag_xscript_on_01gy_1hr_tpm0
          _0_1gy_3=flag_xscript_on_01gy_3hr_tpm0
          _0_1gy_24=flag_xscript_on_01gy_24hr_tpm0
          _0_1gy_72=flag_xscript_on_01gy_72hr_tpm0
          _1gy_1=flag_xscript_on_1gy_1hr_tpm0
          _1gy_3=flag_xscript_on_1gy_3hr_tpm0
          _1gy_24=flag_xscript_on_1gy_24hr_tpm0
          _1gy_72=flag_xscript_on_1gy_72hr_tpm0
          Mock_1=flag_xscript_on_mock_1hr_tpm0
          Mock_3=flag_xscript_on_mock_3hr_tpm0
          Mock_24=flag_xscript_on_mock_24hr_tpm0
          Mock_72=flag_xscript_on_mock_72hr_tpm0;
run;

data trttime_tpm5;
   set flag_iso_trttime_sbys_tpm5;
   drop _NAME_;
   rename _0_1gy_1=flag_xscript_on_01gy_1hr_tpm5
          _0_1gy_3=flag_xscript_on_01gy_3hr_tpm5
          _0_1gy_24=flag_xscript_on_01gy_24hr_tpm5
          _0_1gy_72=flag_xscript_on_01gy_72hr_tpm5
          _1gy_1=flag_xscript_on_1gy_1hr_tpm5
          _1gy_3=flag_xscript_on_1gy_3hr_tpm5
          _1gy_24=flag_xscript_on_1gy_24hr_tpm5
          _1gy_72=flag_xscript_on_1gy_72hr_tpm5
          Mock_1=flag_xscript_on_mock_1hr_tpm5
          Mock_3=flag_xscript_on_mock_3hr_tpm5
          Mock_24=flag_xscript_on_mock_24hr_tpm5
          Mock_72=flag_xscript_on_mock_72hr_tpm5;
run;

/* Merge and make permenant */

proc sort data=treat_tpm0;
   by transcript_id;
proc sort data=time_tpm0;
   by transcript_id;
proc sort data=trttime_tpm0;
   by transcript_id;
proc sort data=treat_tpm5;
   by transcript_id;
proc sort data=time_tpm5;
   by transcript_id;
proc sort data=trttime_tpm5;
   by transcript_id;
run;

data rs.arab_flag_transcript_on_gt0;
  merge treat_tpm0 (in=in1) time_tpm0 (in=in2) trttime_tpm0 (in=in3);
  by transcript_id;
  if in1 and in2 and in3;
run;

data rs.arab_flag_transcript_on_ge5;
  merge treat_tpm5 (in=in1) time_tpm5 (in=in2) trttime_tpm5 (in=in3);
  by transcript_id;
  if in1 and in2 and in3;
run;



