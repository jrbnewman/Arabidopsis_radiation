ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";

/* Flag genes on/off, CPM>0 in at least 50% of samples per treatment */

data counts_w_key;
  set rs.cpm_norm_counts_by_gene;
  if cpm > 0 then flag_cpm_gt0=1;
  else flag_cpm_gt0=0;

  if cpm ge 5 then flag_cpm_ge5=1;
  else flag_cpm_ge5=0;
run;


/* Flag on per treatment */

proc sort data=counts_w_key;
   by treatment gene_id;
proc means data=counts_w_key noprint;
   by treatment gene_id;
   var flag_cpm_gt0 flag_cpm_ge5;
   output out=mean_on_by_trt mean=;
run;

data flag_iso_on_by_trt;
  set mean_on_by_trt;

  /* For CPM > 0 */
  if flag_cpm_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_cpm_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  /* For CPM ge 5 */
  if flag_cpm_ge5 ge 0.5 then flag_gene_on_ge5=1;
  else if flag_cpm_ge5 = 0 then flag_gene_on_ge5=0;
  else flag_gene_on_ge5=.;

  keep gene_id treatment flag_gene_on_gt0 flag_gene_on_ge5;
run;

proc sort data=flag_iso_on_by_trt;
   by gene_id treatment;
proc transpose data=flag_iso_on_by_trt out=flag_iso_trt_sbys_cpm0;
   by gene_id;
   id treatment;
   var flag_gene_on_gt0;
run;

proc transpose data=flag_iso_on_by_trt out=flag_iso_trt_sbys_cpm5;
   by gene_id;
   id treatment;
   var flag_gene_on_ge5;
run;

data treat_cpm0;
   set flag_iso_trt_sbys_cpm0;
   drop _NAME_;
   rename _0_1gy=flag_gene_on_01gy_cpm0
          _1gy=flag_gene_on_1gy_cpm0
          Mock=flag_gene_on_mock_cpm0;
run;

data treat_cpm5;
   set flag_iso_trt_sbys_cpm5;
   drop _NAME_;
   rename _0_1gy=flag_gene_on_01gy_cpm5
          _1gy=flag_gene_on_1gy_cpm5
          Mock=flag_gene_on_mock_cpm5;
run;

/* Flag on per time */

proc sort data=counts_w_key;
   by time gene_id;
proc means data=counts_w_key noprint;
   by time gene_id;
   var flag_cpm_gt0 flag_cpm_ge5;
   output out=mean_on_by_time mean=;
run;

data flag_iso_on_by_time;
  set mean_on_by_time;

  /* For APN > 0 */
  if flag_cpm_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_cpm_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  /* For APN ge 5 */
  if flag_cpm_ge5 ge 0.5 then flag_gene_on_ge5=1;
  else if flag_cpm_ge5 = 0 then flag_gene_on_ge5=0;
  else flag_gene_on_ge5=.;

  keep gene_id time flag_gene_on_gt0 flag_gene_on_ge5;
run;

proc sort data=flag_iso_on_by_time;
   by gene_id time;
proc transpose data=flag_iso_on_by_time out=flag_iso_time_sbys_cpm0;
   by gene_id;
   id time;
   var flag_gene_on_gt0;
run;

proc transpose data=flag_iso_on_by_time out=flag_iso_time_sbys_cpm5;
   by gene_id;
   id time;
   var flag_gene_on_ge5;
run;

data time_cpm0;
   set flag_iso_time_sbys_cpm0;
   drop _NAME_;
   rename _1=flag_gene_on_1hr_cpm0
          _3=flag_gene_on_3hr_cpm0
          _24=flag_gene_on_24hr_cpm0
          _72=flag_gene_on_72hr_cpm0;
run;

data time_cpm5;
   set flag_iso_time_sbys_cpm5;
   drop _NAME_;
   rename _1=flag_gene_on_1hr_cpm5
          _3=flag_gene_on_3hr_cpm5
          _24=flag_gene_on_24hr_cpm5
          _72=flag_gene_on_72hr_cpm5;
run;

/* Flag on per treatment-by-time */

proc sort data=counts_w_key;
   by treatment time gene_id;
proc means data=counts_w_key noprint;
   by treatment time gene_id;
   var flag_cpm_gt0 flag_cpm_ge5;
   output out=mean_on_by_trttime mean=;
run;

data flag_iso_on_by_trttime;
  length condition $10.;
  set mean_on_by_trttime;
  condition=catx("_",treatment,time);

  /* For TPM > 0 */
  if flag_cpm_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_cpm_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  /* For TPM ge 5 */
  if flag_cpm_ge5 ge 0.5 then flag_gene_on_ge5=1;
  else if flag_cpm_ge5 = 0 then flag_gene_on_ge5=0;
  else flag_gene_on_ge5=.;

  keep gene_id condition flag_gene_on_gt0 flag_gene_on_ge5;
run;

proc sort data=flag_iso_on_by_trttime;
   by gene_id condition;
proc transpose data=flag_iso_on_by_trttime out=flag_iso_trttime_sbys_cpm0;
   by gene_id;
   id condition;
   var flag_gene_on_gt0;
run;

proc transpose data=flag_iso_on_by_trttime out=flag_iso_trttime_sbys_cpm5;
   by gene_id;
   id condition;
   var flag_gene_on_ge5;
run;

data trttime_cpm0;
   set flag_iso_trttime_sbys_cpm0;
   drop _NAME_;
   rename _0_1gy_1=flag_gene_on_01gy_1hr_cpm0
          _0_1gy_3=flag_gene_on_01gy_3hr_cpm0
          _0_1gy_24=flag_gene_on_01gy_24hr_cpm0
          _0_1gy_72=flag_gene_on_01gy_72hr_cpm0
          _1gy_1=flag_gene_on_1
data rs.arab_flag_gene_on_cpm_gt0;
  merge treat_cpm0 (in=in1) time_cpm0 (in=in2) trttime_cpm0 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

data rs.arab_flag_gene_on_cpm_ge5;
  merge treat_cpm5 (in=in1) time_cpm5 (in=in2) trttime_cpm5 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;
gy_1hr_cpm0
          _1gy_3=flag_gene_on_1gy_3hr_cpm0
          _1gy_24=flag_gene_on_1gy_24hr_cpm0
          _1gy_72=flag_gene_on_1gy_72hr_cpm0
          Mock_1=flag_gene_on_mock_1hr_cpm0
          Mock_3=flag_gene_on_mock_3hr_cpm0
          Mock_24=flag_gene_on_mock_24hr_cpm0
          Mock_72=flag_gene_on_mock_72hr_cpm0;
run;

data trttime_cpm5;
   set flag_iso_trttime_sbys_cpm5;
   drop _NAME_;
   rename _0_1gy_1=flag_gene_on_01gy_1hr_cpm5
          _0_1gy_3=flag_gene_on_01gy_3hr_cpm5
          _0_1gy_24=flag_gene_on_01gy_24hr_cpm5
          _0_1gy_72=flag_gene_on_01gy_72hr_cpm5
          _1gy_1=flag_gene_on_1gy_1hr_cpm5
          _1gy_3=flag_gene_on_1gy_3hr_cpm5
          _1gy_24=flag_gene_on_1gy_24hr_cpm5
          _1gy_72=flag_gene_on_1gy_72hr_cpm5
          Mock_1=flag_gene_on_mock_1hr_cpm5
          Mock_3=flag_gene_on_mock_3hr_cpm5
          Mock_24=flag_gene_on_mock_24hr_cpm5
          Mock_72=flag_gene_on_mock_72hr_cpm5;
run;

/* Merge and make permenant */

proc sort data=treat_cpm0;
   by gene_id;
proc sort data=time_cpm0;
   by gene_id;
proc sort data=trttime_cpm0;
   by gene_id;
proc sort data=treat_cpm5;
   by gene_id;
proc sort data=time_cpm5;
   by gene_id;
proc sort data=trttime_cpm5;
   by gene_id;
run;

data rs.arab_flag_gene_on_cpm_gt0;
  merge treat_cpm0 (in=in1) time_cpm0 (in=in2) trttime_cpm0 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

data rs.arab_flag_gene_on_cpm_ge5;
  merge treat_cpm5 (in=in1) time_cpm5 (in=in2) trttime_cpm5 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

/* Count CPM > 0 */

ods listing;
proc freq data=rs.arab_flag_gene_on_cpm_gt0;
   tables flag_gene_on_01gy_cpm0 flag_gene_on_1gy_cpm0 flag_gene_on_Mock_cpm0;
run;

/*
                          The FREQ Procedure

   flag_gene_on_                             Cumulative    Cumulative
       01gy_cpm0    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        8651       27.88          8651        27.88
               1       22378       72.12         31029       100.00

                       Frequency Missing = 2948


   flag_gene_on_                             Cumulative    Cumulative
        1gy_cpm0    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        8520       27.46          8520        27.46
               1       22506       72.54         31026       100.00

                       Frequency Missing = 2951


   flag_gene_on_                             Cumulative    Cumulative
       mock_cpm0    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        8682       27.94          8682        27.94
               1       22387       72.06         31069       100.00

                       Frequency Missing = 2908

*/


