ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Flag fusions on/off at APN>0 and APN>5 */

data counts_w_key;
  set rs.arab_fusion_counts_w_key;
  if apn > 0 then flag_apn_gt0=1;
  else flag_apn_gt0=0;

  if apn ge 5 then flag_apn_ge5=1;
  else flag_apn_ge5=0;
run;

/* Flag on per treatment */

proc sort data=counts_w_key;
   by treatment fusion_id;
proc means data=counts_w_key noprint;
   by treatment fusion_id;
   var flag_apn_gt0 flag_apn_ge5;
   output out=mean_on_by_trt mean=;
run;

data flag_fus_on_by_trt;
  set mean_on_by_trt;

  /* For APN > 0 */
  if flag_apn_gt0 ge 0.5 then flag_fusion_on_gt0=1;
  else if flag_apn_gt0 = 0 then flag_fusion_on_gt0=0;
  else flag_fusion_on_gt0=.;

  /* For APN ge 5 */
  if flag_apn_ge5 ge 0.5 then flag_fusion_on_ge5=1;
  else if flag_apn_ge5 = 0 then flag_fusion_on_ge5=0;
  else flag_fusion_on_ge5=.;

  keep fusion_id treatment flag_fusion_on_gt0 flag_fusion_on_ge5;
run;

proc sort data=flag_fus_on_by_trt;
   by fusion_id treatment;
proc transpose data=flag_fus_on_by_trt out=flag_fus_trt_sbys_apn0;
   by fusion_id;
   id treatment;
   var flag_fusion_on_gt0;
run;

proc transpose data=flag_fus_on_by_trt out=flag_fus_trt_sbys_apn5;
   by fusion_id;
   id treatment;
   var flag_fusion_on_ge5;
run;

data treat_apn0;
   set flag_fus_trt_sbys_apn0;
   drop _NAME_;
   rename _0_1gy=flag_fusion_on_01gy_apn0
          _1gy=flag_fusion_on_1gy_apn0
          Mock=flag_fusion_on_mock_apn0;
run;

data treat_apn5;
   set flag_fus_trt_sbys_apn5;
   drop _NAME_;
   rename _0_1gy=flag_fusion_on_01gy_apn5
          _1gy=flag_fusion_on_1gy_apn5
          Mock=flag_fusion_on_mock_apn5;
run;

/* Flag on per time */

proc sort data=counts_w_key;
   by time fusion_id;
proc means data=counts_w_key noprint;
   by time fusion_id;
   var flag_apn_gt0 flag_apn_ge5;
   output out=mean_on_by_time mean=;
run;

data flag_fus_on_by_time;
  set mean_on_by_time;

  /* For APN > 0 */
  if flag_apn_gt0 ge 0.5 then flag_fusion_on_gt0=1;
  else if flag_apn_gt0 = 0 then flag_fusion_on_gt0=0;
  else flag_fusion_on_gt0=.;

  /* For APN ge 5 */
  if flag_apn_ge5 ge 0.5 then flag_fusion_on_ge5=1;
  else if flag_apn_ge5 = 0 then flag_fusion_on_ge5=0;
  else flag_fusion_on_ge5=.;

  keep fusion_id time flag_fusion_on_gt0 flag_fusion_on_ge5;
run;

proc sort data=flag_fus_on_by_time;
   by fusion_id time;
proc transpose data=flag_fus_on_by_time out=flag_fus_time_sbys_apn0;
   by fusion_id;
   id time;
   var flag_fusion_on_gt0;
run;

proc transpose data=flag_fus_on_by_time out=flag_fus_time_sbys_apn5;
   by fusion_id;
   id time;
   var flag_fusion_on_ge5;
run;

data time_apn0;
   set flag_fus_time_sbys_apn0;
   drop _NAME_;
   rename _1=flag_fusion_on_1hr_apn0
          _3=flag_fusion_on_3hr_apn0
          _24=flag_fusion_on_34hr_apn0
          _72=flag_fusion_on_72hr_apn0;
run;

data time_apn5;
   set flag_fus_time_sbys_apn5;
   drop _NAME_;
   rename _1=flag_fusion_on_1hr_apn5
          _3=flag_fusion_on_3hr_apn5
          _24=flag_fusion_on_34hr_apn5
          _72=flag_fusion_on_72hr_apn5;
run;

/* Flag on per treatment-by-time */

proc sort data=counts_w_key;
   by treatment time fusion_id;
proc means data=counts_w_key noprint;
   by treatment time fusion_id;
   var flag_apn_gt0 flag_apn_ge5;
   output out=mean_on_by_trttime mean=;
run;

data flag_fus_on_by_trttime;
  length condition $10.;
  set mean_on_by_trttime;
  condition=catx("_",treatment,time);

  /* For APN > 0 */
  if flag_apn_gt0 ge 0.5 then flag_fusion_on_gt0=1;
  else if flag_apn_gt0 = 0 then flag_fusion_on_gt0=0;
  else flag_fusion_on_gt0=.;

  /* For APN ge 5 */
  if flag_apn_ge5 ge 0.5 then flag_fusion_on_ge5=1;
  else if flag_apn_ge5 = 0 then flag_fusion_on_ge5=0;
  else flag_fusion_on_ge5=.;

  keep fusion_id condition flag_fusion_on_gt0 flag_fusion_on_ge5;
run;

proc sort data=flag_fus_on_by_trttime;
   by fusion_id condition;
proc transpose data=flag_fus_on_by_trttime out=flag_fus_trttime_sbys_apn0;
   by fusion_id;
   id condition;
   var flag_fusion_on_gt0;
run;

proc transpose data=flag_fus_on_by_trttime out=flag_fus_trttime_sbys_apn5;
   by fusion_id;
   id condition;
   var flag_fusion_on_ge5;
run;

data trttime_apn0;
   set flag_fus_trttime_sbys_apn0;
   drop _NAME_;
   rename _0_1gy_1=flag_fusion_on_01gy_1hr_apn0
          _0_1gy_3=flag_fusion_on_01gy_3hr_apn0
          _0_1gy_24=flag_fusion_on_01gy_24hr_apn0
          _0_1gy_72=flag_fusion_on_01gy_72hr_apn0
          _1gy_1=flag_fusion_on_1gy_1hr_apn0
          _1gy_3=flag_fusion_on_1gy_3hr_apn0
          _1gy_24=flag_fusion_on_1gy_24hr_apn0
          _1gy_72=flag_fusion_on_1gy_72hr_apn0
          Mock_1=flag_fusion_on_mock_1hr_apn0
          Mock_3=flag_fusion_on_mock_3hr_apn0
          Mock_24=flag_fusion_on_mock_24hr_apn0
          Mock_72=flag_fusion_on_mock_72hr_apn0;
run;

data trttime_apn5;
   set flag_fus_trttime_sbys_apn5;
   drop _NAME_;
   rename _0_1gy_1=flag_fusion_on_01gy_1hr_apn5
          _0_1gy_3=flag_fusion_on_01gy_3hr_apn5
          _0_1gy_24=flag_fusion_on_01gy_24hr_apn5
          _0_1gy_72=flag_fusion_on_01gy_72hr_apn5
          _1gy_1=flag_fusion_on_1gy_1hr_apn5
          _1gy_3=flag_fusion_on_1gy_3hr_apn5
          _1gy_24=flag_fusion_on_1gy_24hr_apn5
          _1gy_72=flag_fusion_on_1gy_72hr_apn5
          Mock_1=flag_fusion_on_mock_1hr_apn5
          Mock_3=flag_fusion_on_mock_3hr_apn5
          Mock_24=flag_fusion_on_mock_24hr_apn5
          Mock_72=flag_fusion_on_mock_72hr_apn5;
run;

/* Merge and make permenant */

proc sort data=treat_apn0;
   by fusion_id;
proc sort data=time_apn0;
   by fusion_id;
proc sort data=trttime_apn0;
   by fusion_id;
proc sort data=treat_apn5;
   by fusion_id;
proc sort data=time_apn5;
   by fusion_id;
proc sort data=trttime_apn5;
   by fusion_id;
run;

data rs.arab_flag_fusion_on_gt0;
  merge treat_apn0 (in=in1) time_apn0 (in=in2) trttime_apn0 (in=in3);
  by fusion_id;
  if in1 and in2 and in3;
run;

data rs.arab_flag_fusion_on_ge5;
  merge treat_apn5 (in=in1) time_apn5 (in=in2) trttime_apn5 (in=in3);
  by fusion_id;
  if in1 and in2 and in3;
run;



