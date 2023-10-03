/* Formatting fusion data for the 13 ethanol response genes for SVN and SS */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

/* Get fusion annotations */

data symbol2gene;
  set tair.gene2symbol_ens;
  if symbol="" then symbol=gene_id;
  keep symbol gene_id;
run;

/* Get LSmeans */

data lsmeans;
  set rs.arab_gene_cntrs_lsmeans_lcpm;
  if time_trt="1_0.1gy" then group="LSmeans_01gy_1h";
  if time_trt="1_1gy" then group="LSmeans_1gy_1h";
  if time_trt="1_Mock" then group="LSmeans_Mock_1h";

  if time_trt="3_0.1gy" then group="LSmeans_01gy_3h";
  if time_trt="3_1gy" then group="LSmeans_1gy_3h";
  if time_trt="3_Mock" then group="LSmeans_Mock_3h";

  if time_trt="24_0.1gy" then group="LSmeans_01gy_24h";
  if time_trt="24_1gy" then group="LSmeans_1gy_24h";
  if time_trt="24_Mock" then group="LSmeans_Mock_24h";

  if time_trt="72_0.1gy" then group="LSmeans_01gy_72h";
  if time_trt="72_1gy" then group="LSmeans_1gy_72h";
  if time_trt="72_Mock" then group="LSmeans_Mock_72h";
  keep gene_id group estimate;
run;

proc sort data=lsmeans;
   by gene_id group;
proc transpose data=lsmeans out=lsmeans_sbys;
   by gene_id;
   id group;
   var estimate;
run;


data lsmeans_sbys2;
   set lsmeans_sbys;
   drop _NAME_ _LABEL_;
run;


/* Get estimates */

data estimates;
  set rs.arab_gene_cntrs_estim_lcpm;
  length est_label $32.;
  if label="0.1gy-Mock: 1h" then est_label="Estimate_01gy_v_Mock_1h";
  else if label= "0.1gy-Mock: 3h" then est_label="Estimate_01gy_v_Mock_3h";
  else if label= "0.1gy-Mock: 24h" then est_label="Estimate_01gy_v_Mock_24h";
  else if label= "0.1gy-Mock: 72h" then est_label="Estimate_01gy_v_Mock_72h";

  else if label= "1gy-Mock: 1h" then est_label="Estimate_1gy_v_Mock_1h";
  else if label= "1gy-Mock: 3h" then est_label="Estimate_1gy_v_Mock_3h";
  else if label= "1gy-Mock: 24h" then est_label="Estimate_1gy_v_Mock_24h";
  else if label= "1gy-Mock: 72h" then est_label="Estimate_1gy_v_Mock_72h";

  else delete;
  keep gene_id est_label estimate;
run;

proc sort data=estimates;
   by gene_id est_label;
proc transpose data=estimates out=estimates_sbys;
   by gene_id;
   id est_label;
   var estimate;
run;


data estimates_sbys2;
   set estimates_sbys;
   drop _NAME_ _LABEL_;
run;


/* detection flags */

data flag_dtct;
   set rs.arab_flag_gene_on_cpm_gt0;
   keep gene_id    flag_gene_on_01gy_cpm0    flag_gene_on_1gy_cpm0
   flag_gene_on_Mock_cpm0    flag_gene_on_1hr_cpm0   flag_gene_on_3hr_cpm0
   flag_gene_on_24hr_cpm0   flag_gene_on_72hr_cpm0   flag_gene_on_01gy_1hr_cpm0
   flag_gene_on_01gy_3hr_cpm0   flag_gene_on_01gy_24hr_cpm0
   flag_gene_on_01gy_72hr_cpm0   flag_gene_on_1gy_1hr_cpm0
   flag_gene_on_1gy_3hr_cpm0   flag_gene_on_1gy_24hr_cpm0
   flag_gene_on_1gy_72hr_cpm0   flag_gene_on_Mock_1hr_cpm0
   flag_gene_on_Mock_3hr_cpm0   flag_gene_on_Mock_24hr_cpm0
   flag_gene_on_Mock_72hr_cpm0;

   rename flag_gene_on_01gy_cpm0=flag_01gy_on
   flag_gene_on_1gy_cpm0=flag_1gy_on
   flag_gene_on_Mock_cpm0=flag_Mock_on
   flag_gene_on_1hr_cpm0=flag_1hr_on
   flag_gene_on_3hr_cpm0=flag_3hr_on
   flag_gene_on_24hr_cpm0=flag_24hr_on
   flag_gene_on_72hr_cpm0=flag_72gy_on
   flag_gene_on_01gy_1hr_cpm0=flag_01gy_1hr_on
   flag_gene_on_01gy_3hr_cpm0=flag_01gy_3hr_on
   flag_gene_on_01gy_24hr_cpm0=flag_01gy_24hr_on
   flag_gene_on_01gy_72hr_cpm0=flag_01gy_72hr_on
   flag_gene_on_1gy_1hr_cpm0=flag_1gy_1hr_on
   flag_gene_on_1gy_3hr_cpm0=flag_1gy_3hr_on
   flag_gene_on_1gy_24hr_cpm0=flag_1gy_24hr_on
   flag_gene_on_1gy_72hr_cpm0=flag_1gy_72hr_on
   flag_gene_on_Mock_1hr_cpm0=flag_Mock_1hr_on
   flag_gene_on_Mock_3hr_cpm0=flag_Mock_3hr_on
   flag_gene_on_Mock_24hr_cpm0=flag_Mock_24hr_on
   flag_gene_on_Mock_72hr_cpm0=flag_Mock_72hr_on;

run;


/* Get gene means */

data means;
  set rs.arab_gene_mean_cpm_by_trt_time;
run;


/* Get Up/Down signs */

data sign;
  set rs.arab_sign_by_cntrst_gene_cpm_fdr;
  drop flag_gene_on_: ;
run;

/* Get FDR corrected P values */

data fdr;
   set rs.fdr_by_gene_log_cpm;
run;


proc sort data=symbol2gene;
   by gene_id;
proc sort data=lsmeans_sbys2;
   by gene_id;
proc sort data=estimates_sbys2;
   by gene_id;
proc sort data=sign;
   by gene_id;
proc sort data=fdr;
   by gene_id;
proc sort data=flag_dtct;
   by gene_id;
proc sort data=means;
   by gene_id;
run;

data gene_summary gene_skipped oops;
   merge symbol2gene (in=in1) flag_dtct (in=in2) fdr lsmeans_sbys2 estimates_sbys2 sign means;
   by gene_id;
   if in2 then output gene_summary;
   else if in1 then output gene_skipped;
run;

/* Add in log2 FC */

data gene_summary2;
  set gene_summary;
   if symbol="" then symbol=gene_id;
      if fdr_01gy_v_Mock_1h ne . then log2_fc_01gy_v_Mock_lh=log2(mean_cpm_01gy_1h/mean_cpm_Mock_1h);
      if fdr_01gy_v_Mock_3h ne . then log2_fc_01gy_v_Mock_3h=log2(mean_cpm_01gy_3h/mean_cpm_Mock_3h);
      if fdr_01gy_v_Mock_24h ne . then log2_fc_01gy_v_Mock_24h=log2(mean_cpm_01gy_24h/mean_cpm_Mock_24h);
      if fdr_01gy_v_Mock_72h ne . then log2_fc_01gy_v_Mock_72h=log2(mean_cpm_01gy_72h/mean_cpm_Mock_72h);

      if fdr_1gy_v_Mock_1h ne . then log2_fc_1gy_v_Mock_lh=log2(mean_cpm_1gy_1h/mean_cpm_Mock_1h);
      if fdr_1gy_v_Mock_3h ne . then log2_fc_1gy_v_Mock_3h=log2(mean_cpm_1gy_3h/mean_cpm_Mock_3h);
      if fdr_1gy_v_Mock_24h ne . then log2_fc_1gy_v_Mock_24h=log2(mean_cpm_1gy_24h/mean_cpm_Mock_24h);
      if fdr_1gy_v_Mock_72h ne . then log2_fc_1gy_v_Mock_72h=log2(mean_cpm_1gy_72h/mean_cpm_Mock_72h);
run;


data check;
  set gene_summary2;
  if log2_fc_01gy_v_Mock_24h = . then fc=0;
  else if log2_fc_01gy_v_Mock_24h > 0 then fc=1;
  else if log2_fc_01gy_v_Mock_24h < 0 then fc=-1;
  else fc=.;
run;

proc freq data=check noprint;
   tables fc*sign_01gy_v_Mock_24h_fdr20 / out=check2;
proc print data=check2;
run;

proc freq data=check noprint;
   tables flag_Mock_on*flag_01gy_on*flag_1gy_on / out=check2;
proc print data=check2;
run;

data check3;
  set check;
  if flag_mock_on=1 then delete;
  if flag_01gy_on=1 then delete;
  if flag_1gy_on=1 then delete;
run;

data check4;
  set gene_summary2;
  if sign_01gy_v_Mock_24h_fdr20="";
run;

/* Reorder variables */

proc contents data=gene_summary2;
run;

data gene_summary3;
retain
gene_id symbol transcript_id_s_
flag_Mock_on flag_01gy_on flag_1gy_on
flag_1hr_on flag_3hr_on flag_24hr_on flag_72gy_on
flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
P_time fdr_time flag_time_fdr05 flag_time_fdr10 flag_time_fdr20
P_treatment fdr_treatment flag_treatment_fdr05 flag_treatment_fdr10 flag_treatment_fdr20
P_trt_by_time fdr_trt_by_time flag_trt_by_time_fdr05 flag_trt_by_time_fdr10 flag_trt_by_time_fdr20

mean_cpm_mock_1h mean_cpm_mock_3h mean_cpm_mock_24h mean_cpm_mock_72h
mean_cpm_01gy_1h mean_cpm_01gy_3h mean_cpm_01gy_24h mean_cpm_01gy_72h 
mean_cpm_1gy_1h mean_cpm_1gy_3h mean_cpm_1gy_24h mean_cpm_1gy_72h

LSmeans_Mock_1h LSmeans_Mock_3h LSmeans_Mock_24 LSmeans_Mock_72
LSmeans_01gy_1h LSmeans_01gy_3h LSmeans_01gy_24 LSmeans_01gy_72
LSmeans_1gy_1h LSmeans_1gy_3h LSmeans_1gy_24h LSmeans_1gy_72h

flag_gene_01gy_v_Mock_1h_dd
Estimate_01gy_v_Mock_1h log2_fc_01gy_v_Mock_lh p_01gy_v_Mock_1h fdr_01gy_v_Mock_1h
flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_1h_fdr10 flag_01gy_v_Mock_1h_fdr20
sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_1h_fdr10 sign_01gy_v_Mock_1h_fdr20

flag_gene_01gy_v_Mock_3h_dd
Estimate_01gy_v_Mock_3h log2_fc_01gy_v_Mock_3h p_01gy_v_Mock_3h fdr_01gy_v_Mock_3h
flag_01gy_v_Mock_3h_fdr05 flag_01gy_v_Mock_3h_fdr10 flag_01gy_v_Mock_3h_fdr20
sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_3h_fdr10 sign_01gy_v_Mock_3h_fdr20

flag_gene_01gy_v_Mock_24h_dd
Estimate_01gy_v_Mock_24h log2_fc_01gy_v_Mock_24h p_01gy_v_Mock_24h fdr_01gy_v_Mock_24h
flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_24h_fdr10 flag_01gy_v_Mock_24h_fdr20
sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_24h_fdr10 sign_01gy_v_Mock_24h_fdr20

flag_gene_01gy_v_Mock_72h_dd
Estimate_01gy_v_Mock_72h log2_fc_01gy_v_Mock_72h p_01gy_v_Mock_72h fdr_01gy_v_Mock_72h
flag_01gy_v_Mock_72h_fdr05 flag_01gy_v_Mock_72h_fdr10 flag_01gy_v_Mock_72h_fdr20 
sign_01gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_72h_fdr10 sign_01gy_v_Mock_72h_fdr20

flag_gene_1gy_v_Mock_1h_dd
Estimate_1gy_v_Mock_1h log2_fc_1gy_v_Mock_lh p_1gy_v_Mock_1h fdr_1gy_v_Mock_1h
flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_1h_fdr10 flag_1gy_v_Mock_1h_fdr20
sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_1h_fdr10 sign_1gy_v_Mock_1h_fdr20

flag_gene_1gy_v_Mock_3h_dd
Estimate_1gy_v_Mock_3h log2_fc_1gy_v_Mock_3h p_1gy_v_Mock_3h fdr_1gy_v_Mock_3h
flag_1gy_v_Mock_3h_fdr05 flag_1gy_v_Mock_3h_fdr10 flag_1gy_v_Mock_3h_fdr20
sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_3h_fdr10 sign_1gy_v_Mock_3h_fdr20

flag_gene_1gy_v_Mock_24h_dd
Estimate_1gy_v_Mock_24h log2_fc_1gy_v_Mock_24h p_1gy_v_Mock_24h fdr_1gy_v_Mock_24h
flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_24h_fdr10 flag_1gy_v_Mock_24h_fdr20
sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_24h_fdr10 sign_1gy_v_Mock_24h_fdr20

flag_gene_1gy_v_Mock_72h_dd
Estimate_1gy_v_Mock_72h log2_fc_1gy_v_Mock_72h p_1gy_v_Mock_72h fdr_1gy_v_Mock_72h
flag_1gy_v_Mock_72h_fdr05 flag_1gy_v_Mock_72h_fdr10 flag_1gy_v_Mock_72h_fdr20
sign_1gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_72h_fdr10 sign_1gy_v_Mock_72h_fdr20

;
  
         set gene_summary2;


if flag_Mock_on=. then flag_Mock_on=0;
if flag_01gy_on=. then flag_01gy_on=0;
if flag_1gy_on=. then flag_1gy_on=0;
if flag_1hr_on=. then flag_1hr_on=0;
if flag_3hr_on=. then flag_3hr_on=0;
if flag_24hr_on=. then flag_24hr_on=0;
if flag_72gy_on=. then flag_72gy_on=0;

  if flag_Mock_1hr_on=. then flag_Mock_1hr_on=0;
  if flag_Mock_3hr_on=. then flag_Mock_3hr_on=0;
  if flag_Mock_24hr_on=. then flag_Mock_24hr_on=0;
  if flag_Mock_72hr_on=. then flag_Mock_72hr_on=0;
  if flag_01gy_1hr_on=. then flag_01gy_1hr_on=0;
  if flag_01gy_3hr_on=. then flag_01gy_3hr_on=0;
  if flag_01gy_24hr_on=. then flag_01gy_24hr_on=0;
  if flag_01gy_72hr_on=. then flag_01gy_72hr_on=0;
  if flag_1gy_1hr_on=. then flag_1gy_1hr_on=0;
  if flag_1gy_3hr_on=. then flag_1gy_3hr_on=0;
  if flag_1gy_24hr_on=. then flag_1gy_24hr_on=0;
  if flag_1gy_72hr_on=. then flag_1gy_72hr_on=0;

  keep
gene_id symbol transcript_id_s_
flag_Mock_on flag_01gy_on flag_1gy_on
flag_1hr_on flag_3hr_on flag_24hr_on flag_72gy_on
flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
P_time fdr_time flag_time_fdr05 flag_time_fdr10 flag_time_fdr20
P_treatment fdr_treatment flag_treatment_fdr05 flag_treatment_fdr10 flag_treatment_fdr20
P_trt_by_time fdr_trt_by_time flag_trt_by_time_fdr05 flag_trt_by_time_fdr10 flag_trt_by_time_fdr20

mean_cpm_mock_1h mean_cpm_mock_3h mean_cpm_mock_24h mean_cpm_mock_72h
mean_cpm_01gy_1h mean_cpm_01gy_3h mean_cpm_01gy_24h mean_cpm_01gy_72h 
mean_cpm_1gy_1h mean_cpm_1gy_3h mean_cpm_1gy_24h mean_cpm_1gy_72h

LSmeans_Mock_1h LSmeans_Mock_3h LSmeans_Mock_24 LSmeans_Mock_72
LSmeans_01gy_1h LSmeans_01gy_3h LSmeans_01gy_24 LSmeans_01gy_72
LSmeans_1gy_1h LSmeans_1gy_3h LSmeans_1gy_24h LSmeans_1gy_72h

flag_gene_01gy_v_Mock_1h_dd
Estimate_01gy_v_Mock_1h log2_fc_01gy_v_Mock_lh p_01gy_v_Mock_1h fdr_01gy_v_Mock_1h
flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_1h_fdr10 flag_01gy_v_Mock_1h_fdr20
sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_1h_fdr10 sign_01gy_v_Mock_1h_fdr20

flag_gene_01gy_v_Mock_3h_dd
Estimate_01gy_v_Mock_3h log2_fc_01gy_v_Mock_3h p_01gy_v_Mock_3h fdr_01gy_v_Mock_3h
flag_01gy_v_Mock_3h_fdr05 flag_01gy_v_Mock_3h_fdr10 flag_01gy_v_Mock_3h_fdr20
sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_3h_fdr10 sign_01gy_v_Mock_3h_fdr20

flag_gene_01gy_v_Mock_24h_dd
Estimate_01gy_v_Mock_24h log2_fc_01gy_v_Mock_24h p_01gy_v_Mock_24h fdr_01gy_v_Mock_24h
flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_24h_fdr10 flag_01gy_v_Mock_24h_fdr20
sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_24h_fdr10 sign_01gy_v_Mock_24h_fdr20

flag_gene_01gy_v_Mock_72h_dd
Estimate_01gy_v_Mock_72h log2_fc_01gy_v_Mock_72h p_01gy_v_Mock_72h fdr_01gy_v_Mock_72h
flag_01gy_v_Mock_72h_fdr05 flag_01gy_v_Mock_72h_fdr10 flag_01gy_v_Mock_72h_fdr20 
sign_01gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_72h_fdr10 sign_01gy_v_Mock_72h_fdr20

flag_gene_1gy_v_Mock_1h_dd
Estimate_1gy_v_Mock_1h log2_fc_1gy_v_Mock_lh p_1gy_v_Mock_1h fdr_1gy_v_Mock_1h
flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_1h_fdr10 flag_1gy_v_Mock_1h_fdr20
sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_1h_fdr10 sign_1gy_v_Mock_1h_fdr20

flag_gene_1gy_v_Mock_3h_dd
Estimate_1gy_v_Mock_3h log2_fc_1gy_v_Mock_3h p_1gy_v_Mock_3h fdr_1gy_v_Mock_3h
flag_1gy_v_Mock_3h_fdr05 flag_1gy_v_Mock_3h_fdr10 flag_1gy_v_Mock_3h_fdr20
sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_3h_fdr10 sign_1gy_v_Mock_3h_fdr20

flag_gene_1gy_v_Mock_24h_dd
Estimate_1gy_v_Mock_24h log2_fc_1gy_v_Mock_24h p_1gy_v_Mock_24h fdr_1gy_v_Mock_24h
flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_24h_fdr10 flag_1gy_v_Mock_24h_fdr20
sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_24h_fdr10 sign_1gy_v_Mock_24h_fdr20

flag_gene_1gy_v_Mock_72h_dd
Estimate_1gy_v_Mock_72h log2_fc_1gy_v_Mock_72h p_1gy_v_Mock_72h fdr_1gy_v_Mock_72h
flag_1gy_v_Mock_72h_fdr05 flag_1gy_v_Mock_72h_fdr10 flag_1gy_v_Mock_72h_fdr20
sign_1gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_72h_fdr10 sign_1gy_v_Mock_72h_fdr20
;
  

if sign_01gy_v_Mock_1h_fdr05="" then sign_01gy_v_Mock_1h_fdr05="Off";
if sign_01gy_v_Mock_1h_fdr10="" then sign_01gy_v_Mock_1h_fdr10="Off";
if sign_01gy_v_Mock_1h_fdr20="" then sign_01gy_v_Mock_1h_fdr20="Off";

if sign_01gy_v_Mock_3h_fdr05="" then sign_01gy_v_Mock_3h_fdr05="Off";
if sign_01gy_v_Mock_3h_fdr10="" then sign_01gy_v_Mock_3h_fdr10="Off";
if sign_01gy_v_Mock_3h_fdr20="" then sign_01gy_v_Mock_3h_fdr20="Off";

if sign_01gy_v_Mock_24h_fdr05="" then sign_01gy_v_Mock_24h_fdr05="Off";
if sign_01gy_v_Mock_24h_fdr10="" then sign_01gy_v_Mock_24h_fdr10="Off";
if sign_01gy_v_Mock_24h_fdr20="" then sign_01gy_v_Mock_24h_fdr20="Off";

if sign_01gy_v_Mock_72h_fdr05="" then sign_01gy_v_Mock_72h_fdr05="Off";
if sign_01gy_v_Mock_72h_fdr10="" then sign_01gy_v_Mock_72h_fdr10="Off";
if sign_01gy_v_Mock_72h_fdr20="" then sign_01gy_v_Mock_72h_fdr20="Off";

if sign_1gy_v_Mock_1h_fdr05="" then sign_1gy_v_Mock_1h_fdr05="Off";
if sign_1gy_v_Mock_1h_fdr10="" then sign_1gy_v_Mock_1h_fdr10="Off";
if sign_1gy_v_Mock_1h_fdr20="" then sign_1gy_v_Mock_1h_fdr20="Off";

if sign_1gy_v_Mock_3h_fdr05="" then sign_1gy_v_Mock_3h_fdr05="Off";
if sign_1gy_v_Mock_3h_fdr10="" then sign_1gy_v_Mock_3h_fdr10="Off";
if sign_1gy_v_Mock_3h_fdr20="" then sign_1gy_v_Mock_3h_fdr20="Off";

if sign_1gy_v_Mock_24h_fdr05="" then sign_1gy_v_Mock_24h_fdr05="Off";
if sign_1gy_v_Mock_24h_fdr10="" then sign_1gy_v_Mock_24h_fdr10="Off";
if sign_1gy_v_Mock_24h_fdr20="" then sign_1gy_v_Mock_24h_fdr20="Off";

if sign_1gy_v_Mock_72h_fdr05="" then sign_1gy_v_Mock_72h_fdr05="Off";
if sign_1gy_v_Mock_72h_fdr10="" then sign_1gy_v_Mock_72h_fdr10="Off";
if sign_1gy_v_Mock_72h_fdr20="" then sign_1gy_v_Mock_72h_fdr20="Off";


;
  rename flag_72gy_on=flag_72hr_on
         LSmeans_Mock_24=LSmeans_Mock_24h
         LSmeans_Mock_72=LSmeans_Mock_72h
         LSmeans_01gy_24=LSmeans_01gy_24h
         LSmeans_01gy_72=LSmeans_01gy_72h
 ;

run;

/* Export 5%, 10%, 20% FDR separately : full data plus reduced set (gene, symbol, logFC, FDR P, sign), for analyzed only */

data fdr05_full;
   set gene_summary3;
   drop flag_time_fdr10 flag_time_fdr20
flag_treatment_fdr10 flag_treatment_fdr20
flag_trt_by_time_fdr10 flag_trt_by_time_fdr20
flag_01gy_v_Mock_1h_fdr10 flag_01gy_v_Mock_1h_fdr20
sign_01gy_v_Mock_1h_fdr10 sign_01gy_v_Mock_1h_fdr20
flag_01gy_v_Mock_3h_fdr10 flag_01gy_v_Mock_3h_fdr20
sign_01gy_v_Mock_3h_fdr10 sign_01gy_v_Mock_3h_fdr20
flag_01gy_v_Mock_24h_fdr10 flag_01gy_v_Mock_24h_fdr20
sign_01gy_v_Mock_24h_fdr10 sign_01gy_v_Mock_24h_fdr20
flag_01gy_v_Mock_72h_fdr10 flag_01gy_v_Mock_72h_fdr20 
sign_01gy_v_Mock_72h_fdr10 sign_01gy_v_Mock_72h_fdr20
flag_1gy_v_Mock_1h_fdr10 flag_1gy_v_Mock_1h_fdr20
sign_1gy_v_Mock_1h_fdr10 sign_1gy_v_Mock_1h_fdr20
flag_1gy_v_Mock_3h_fdr10 flag_1gy_v_Mock_3h_fdr20
sign_1gy_v_Mock_3h_fdr10 sign_1gy_v_Mock_3h_fdr20
flag_1gy_v_Mock_24h_fdr10 flag_1gy_v_Mock_24h_fdr20
sign_1gy_v_Mock_24h_fdr10 sign_1gy_v_Mock_24h_fdr20
flag_1gy_v_Mock_72h_fdr10 flag_1gy_v_Mock_72h_fdr20
sign_1gy_v_Mock_72h_fdr10 sign_1gy_v_Mock_72h_fdr20
;
run;

data fdr05_reduced;
  set gene_summary3;
  where p_treatment ne . ;
  keep gene_id symbol transcript_id_s_ log2_fc_01gy_v_Mock_lh fdr_01gy_v_Mock_1h flag_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_1h_fdr05
                      log2_fc_01gy_v_Mock_3h fdr_01gy_v_Mock_3h flag_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_3h_fdr05
                      log2_fc_01gy_v_Mock_24h fdr_01gy_v_Mock_24h flag_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_24h_fdr05
                      log2_fc_01gy_v_Mock_72h fdr_01gy_v_Mock_72h flag_01gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_72h_fdr05
                      log2_fc_1gy_v_Mock_lh fdr_1gy_v_Mock_1h flag_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_1h_fdr05
                      log2_fc_1gy_v_Mock_3h fdr_1gy_v_Mock_3h flag_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_3h_fdr05
                      log2_fc_1gy_v_Mock_24h fdr_1gy_v_Mock_24h flag_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_24h_fdr05
                      log2_fc_1gy_v_Mock_72h fdr_1gy_v_Mock_72h flag_1gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_72h_fdr05
                      flag_gene_01gy_v_Mock_1h_dd flag_gene_01gy_v_Mock_3h_dd flag_gene_01gy_v_Mock_24h_dd
                      flag_gene_01gy_v_Mock_72h_dd
                      flag_gene_1gy_v_Mock_1h_dd flag_gene_1gy_v_Mock_3h_dd flag_gene_1gy_v_Mock_24h_dd
                      flag_gene_1gy_v_Mock_72h_dd
;
run;


data fdr10_full;
   set gene_summary3;
   drop flag_time_fdr05 flag_time_fdr20
flag_treatment_fdr05 flag_treatment_fdr20
flag_trt_by_time_fdr05 flag_trt_by_time_fdr20
flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_1h_fdr20
sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_1h_fdr20
flag_01gy_v_Mock_3h_fdr05 flag_01gy_v_Mock_3h_fdr20
sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_3h_fdr20
flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_24h_fdr20
sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_24h_fdr20
flag_01gy_v_Mock_72h_fdr05 flag_01gy_v_Mock_72h_fdr20 
sign_01gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_72h_fdr20
flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_1h_fdr20
sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_1h_fdr20
flag_1gy_v_Mock_3h_fdr05 flag_1gy_v_Mock_3h_fdr20
sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_3h_fdr20
flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_24h_fdr20
sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_24h_fdr20
flag_1gy_v_Mock_72h_fdr05 flag_1gy_v_Mock_72h_fdr20
sign_1gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_72h_fdr20
;
run;

data fdr10_reduced;
  set gene_summary3;
  where p_treatment ne . ;
  keep gene_id symbol transcript_id_s_ log2_fc_01gy_v_Mock_lh fdr_01gy_v_Mock_1h flag_01gy_v_Mock_1h_fdr10 sign_01gy_v_Mock_1h_fdr10
                      log2_fc_01gy_v_Mock_3h fdr_01gy_v_Mock_3h flag_01gy_v_Mock_3h_fdr10 sign_01gy_v_Mock_3h_fdr10
                      log2_fc_01gy_v_Mock_24h fdr_01gy_v_Mock_24h flag_01gy_v_Mock_24h_fdr10 sign_01gy_v_Mock_24h_fdr10
                      log2_fc_01gy_v_Mock_72h fdr_01gy_v_Mock_72h flag_01gy_v_Mock_72h_fdr10 sign_01gy_v_Mock_72h_fdr10
                      log2_fc_1gy_v_Mock_lh fdr_1gy_v_Mock_1h flag_1gy_v_Mock_1h_fdr10 sign_1gy_v_Mock_1h_fdr10
                      log2_fc_1gy_v_Mock_3h fdr_1gy_v_Mock_3h flag_1gy_v_Mock_3h_fdr10 sign_1gy_v_Mock_3h_fdr10
                      log2_fc_1gy_v_Mock_24h fdr_1gy_v_Mock_24h flag_1gy_v_Mock_24h_fdr10 sign_1gy_v_Mock_24h_fdr10
                      log2_fc_1gy_v_Mock_72h fdr_1gy_v_Mock_72h flag_1gy_v_Mock_72h_fdr10 sign_1gy_v_Mock_72h_fdr10
                      flag_gene_01gy_v_Mock_1h_dd flag_gene_01gy_v_Mock_3h_dd flag_gene_01gy_v_Mock_24h_dd
                      flag_gene_01gy_v_Mock_72h_dd
                      flag_gene_1gy_v_Mock_1h_dd flag_gene_1gy_v_Mock_3h_dd flag_gene_1gy_v_Mock_24h_dd
                      flag_gene_1gy_v_Mock_72h_dd
;

run;



data fdr20_full;
   set gene_summary3;
   drop flag_time_fdr05 flag_time_fdr10
flag_treatment_fdr05 flag_treatment_fdr10
flag_trt_by_time_fdr05 flag_trt_by_time_fdr10
flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_1h_fdr10
sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_1h_fdr10
flag_01gy_v_Mock_3h_fdr05 flag_01gy_v_Mock_3h_fdr10
sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_3h_fdr10
flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_24h_fdr10
sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_24h_fdr10
flag_01gy_v_Mock_72h_fdr05 flag_01gy_v_Mock_72h_fdr10 
sign_01gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_72h_fdr10
flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_1h_fdr10
sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_1h_fdr10
flag_1gy_v_Mock_3h_fdr05 flag_1gy_v_Mock_3h_fdr10
sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_3h_fdr10
flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_24h_fdr10
sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_24h_fdr10
flag_1gy_v_Mock_72h_fdr05 flag_1gy_v_Mock_72h_fdr10
sign_1gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_72h_fdr10
;
run;

data fdr20_reduced;
  set gene_summary3;
  where p_treatment ne . ;
  keep gene_id symbol transcript_id_s_ log2_fc_01gy_v_Mock_lh fdr_01gy_v_Mock_1h flag_01gy_v_Mock_1h_fdr20 sign_01gy_v_Mock_1h_fdr20
                      log2_fc_01gy_v_Mock_3h fdr_01gy_v_Mock_3h flag_01gy_v_Mock_3h_fdr20 sign_01gy_v_Mock_3h_fdr20
                      log2_fc_01gy_v_Mock_24h fdr_01gy_v_Mock_24h flag_01gy_v_Mock_24h_fdr20 sign_01gy_v_Mock_24h_fdr20
                      log2_fc_01gy_v_Mock_72h fdr_01gy_v_Mock_72h flag_01gy_v_Mock_72h_fdr20 sign_01gy_v_Mock_72h_fdr20
                      log2_fc_1gy_v_Mock_lh fdr_1gy_v_Mock_1h flag_1gy_v_Mock_1h_fdr20 sign_1gy_v_Mock_1h_fdr20
                      log2_fc_1gy_v_Mock_3h fdr_1gy_v_Mock_3h flag_1gy_v_Mock_3h_fdr20 sign_1gy_v_Mock_3h_fdr20
                      log2_fc_1gy_v_Mock_24h fdr_1gy_v_Mock_24h flag_1gy_v_Mock_24h_fdr20 sign_1gy_v_Mock_24h_fdr20
                      log2_fc_1gy_v_Mock_72h fdr_1gy_v_Mock_72h flag_1gy_v_Mock_72h_fdr20 sign_1gy_v_Mock_72h_fdr20
                      flag_gene_01gy_v_Mock_1h_dd flag_gene_01gy_v_Mock_3h_dd flag_gene_01gy_v_Mock_24h_dd
                      flag_gene_01gy_v_Mock_72h_dd
                      flag_gene_1gy_v_Mock_1h_dd flag_gene_1gy_v_Mock_3h_dd flag_gene_1gy_v_Mock_24h_dd
                      flag_gene_1gy_v_Mock_72h_dd
;

run;





/* Export */

proc export data=fdr05_full outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_cpm_tair10_fdr05_full_results.csv"
dbms=csv replace;
run;
proc export data=fdr05_reduced outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_cpm_tair10_fdr05_simple_results.csv"
dbms=csv replace;
run;


proc export data=fdr10_full outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_cpm_tair10_fdr10_full_results.csv"
dbms=csv replace;
run;
proc export data=fdr10_reduced outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_cpm_tair10_fdr10_simple_results.csv"
dbms=csv replace;
run;


proc export data=fdr20_full outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_cpm_tair10_fdr20_full_results.csv"
dbms=csv replace;
run;
proc export data=fdr20_reduced outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_cpm_tair10_fdr20_simple_results.csv"
dbms=csv replace;
run;


/* Make permanent */

data rs.arab_results_by_gene;
  set gene_summary3;
run;

