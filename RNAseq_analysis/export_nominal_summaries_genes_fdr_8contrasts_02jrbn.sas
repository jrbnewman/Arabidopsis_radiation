/* Formatting fusion data for the 13 ethanol response genes for SVN and SS */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

/* Get fusion annotations */

data symbol2gene;
  set tair.gene2symbol_ens;
  keep symbol gene_id;
run;

data gene_info;
  set tair.tair20_exons;
  keep gene_id start chrom stop strand exon_id;
run;

proc sort data=gene_info nodup;
    by gene_id chrom strand start stop exon_id;
proc means data=gene_info noprint;
    by gene_id chrom strand;
    var start stop;
    output out=gene_coord min(start)=gene_start max(stop)=gene_stop;
run;

/* check that there aren't duplicate genes */

proc freq data=gene_coord noprint;
   tables gene_id /out=gene_check;
proc sort data=gene_check;
   by descending count;
run;
*good, only 1 obs per gene;

proc sort data=symbol2gene;
  by gene_id;
proc sort data=gene_coord;
  by gene_id;
run;

data gene_coord2;
  merge gene_coord (in=in1) symbol2gene (in=in2);
  by gene_id;
  if in1 ;
  if symbol="" then symbol=gene_id;
run;


/* Get LSmeans */

data lsmeans;
  set rs.arab_gene_cntrs_lsmeans;
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
  set rs.arab_gene_cntrs_estim;
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
   set rs.arab_flag_gene_on_gt0;
   keep gene_id    flag_gene_on_01gy_apn0    flag_gene_on_1gy_apn0
   flag_gene_on_Mock_apn0    flag_gene_on_1hr_apn0   flag_gene_on_3hr_apn0
   flag_gene_on_34hr_apn0   flag_gene_on_72hr_apn0   flag_gene_on_01gy_1hr_apn0
   flag_gene_on_01gy_3hr_apn0   flag_gene_on_01gy_24hr_apn0
   flag_gene_on_01gy_72hr_apn0   flag_gene_on_1gy_1hr_apn0
   flag_gene_on_1gy_3hr_apn0   flag_gene_on_1gy_24hr_apn0
   flag_gene_on_1gy_72hr_apn0   flag_gene_on_Mock_1hr_apn0
   flag_gene_on_Mock_3hr_apn0   flag_gene_on_Mock_24hr_apn0
   flag_gene_on_Mock_72hr_apn0;

   rename flag_gene_on_01gy_apn0=flag_01gy_on
   flag_gene_on_1gy_apn0=flag_1gy_on
   flag_gene_on_Mock_apn0=flag_Mock_on
   flag_gene_on_1hr_apn0=flag_1hr_on
   flag_gene_on_3hr_apn0=flag_3hr_on
   flag_gene_on_34hr_apn0=flag_24hr_on
   flag_gene_on_72hr_apn0=flag_72gy_on
   flag_gene_on_01gy_1hr_apn0=flag_01gy_1hr_on
   flag_gene_on_01gy_3hr_apn0=flag_01gy_3hr_on
   flag_gene_on_01gy_24hr_apn0=flag_01gy_24hr_on
   flag_gene_on_01gy_72hr_apn0=flag_01gy_72hr_on
   flag_gene_on_1gy_1hr_apn0=flag_1gy_1hr_on
   flag_gene_on_1gy_3hr_apn0=flag_1gy_3hr_on
   flag_gene_on_1gy_24hr_apn0=flag_1gy_24hr_on
   flag_gene_on_1gy_72hr_apn0=flag_1gy_72hr_on
   flag_gene_on_Mock_1hr_apn0=flag_Mock_1hr_on
   flag_gene_on_Mock_3hr_apn0=flag_Mock_3hr_on
   flag_gene_on_Mock_24hr_apn0=flag_Mock_24hr_on
   flag_gene_on_Mock_72hr_apn0=flag_Mock_72hr_on;

run;


/* Get gene means */

data means;
  set rs.arab_gene_mean_apn_by_trt_time;
run;


/* Get Up/Down signs */

data sign;
  set rs.arab_sign_by_contrast_gene_fdr2;
  drop flag_gene_on_: ;
run;

/* Get FDR corrected P values */

data fdr;
   set rs.fdr_by_gene_v2;
run;

/* Means */




proc sort data=gene_coord2;
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


data gene_summary gene_skipped;
   length gene_id $444.;
   format gene_id $444.;
   informat gene_id $444.;
   merge gene_coord2 (in=in1) flag_dtct (in=in2) fdr lsmeans_sbys2 estimates_sbys2 sign means;
   by gene_id;
   if in1 and in2 then output gene_summary;
   else if in1 then output gene_skipped;
run;

/* Add in log2 FC */

data gene_summary2;
  set gene_summary;
      if fdr_01gy_v_Mock_1h ne . then log2_fc_01gy_v_Mock_lh=log2(mean_q3apn_01gy_1h/mean_q3apn_Mock_1h);
      if fdr_01gy_v_Mock_3h ne . then log2_fc_01gy_v_Mock_3h=log2(mean_q3apn_01gy_3h/mean_q3apn_Mock_3h);
      if fdr_01gy_v_Mock_24h ne . then log2_fc_01gy_v_Mock_24h=log2(mean_q3apn_01gy_24h/mean_q3apn_Mock_24h);
      if fdr_01gy_v_Mock_72h ne . then log2_fc_01gy_v_Mock_72h=log2(mean_q3apn_01gy_72h/mean_q3apn_Mock_72h);

      if fdr_1gy_v_Mock_1h ne . then log2_fc_1gy_v_Mock_lh=log2(mean_q3apn_1gy_1h/mean_q3apn_Mock_1h);
      if fdr_1gy_v_Mock_3h ne . then log2_fc_1gy_v_Mock_3h=log2(mean_q3apn_1gy_3h/mean_q3apn_Mock_3h);
      if fdr_1gy_v_Mock_24h ne . then log2_fc_1gy_v_Mock_24h=log2(mean_q3apn_1gy_24h/mean_q3apn_Mock_24h);
      if fdr_1gy_v_Mock_72h ne . then log2_fc_1gy_v_Mock_72h=log2(mean_q3apn_1gy_72h/mean_q3apn_Mock_72h);
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




/* Reorder variables */

proc contents data=gene_summary2;
run;

data gene_summary3;
retain
gene_id symbol chrom gene_start gene_stop strand
flag_Mock_on flag_01gy_on flag_1gy_on
flag_1hr_on flag_3hr_on flag_24hr_on flag_72gy_on
flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
P_time fdr_time flag_time_fdr05 flag_time_fdr10 flag_time_fdr20
P_treatment fdr_treatment flag_treatment_fdr05 flag_treatment_fdr10 flag_treatment_fdr20
P_trt_by_time fdr_trt_by_time flag_trt_by_time_fdr05 flag_trt_by_time_fdr10 flag_trt_by_time_fdr20

mean_q3apn_mock_1h mean_q3apn_mock_3h mean_q3apn_mock_24h mean_q3apn_mock_72h
mean_q3apn_01gy_1h mean_q3apn_01gy_3h mean_q3apn_01gy_24h mean_q3apn_01gy_72h 
mean_q3apn_1gy_1h mean_q3apn_1gy_3h mean_q3apn_1gy_24h mean_q3apn_1gy_72h

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
gene_id symbol chrom gene_start gene_stop strand
flag_Mock_on flag_01gy_on flag_1gy_on
flag_1hr_on flag_3hr_on flag_24hr_on flag_72gy_on
flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
P_time fdr_time flag_time_fdr05 flag_time_fdr10 flag_time_fdr20
P_treatment fdr_treatment flag_treatment_fdr05 flag_treatment_fdr10 flag_treatment_fdr20
P_trt_by_time fdr_trt_by_time flag_trt_by_time_fdr05 flag_trt_by_time_fdr10 flag_trt_by_time_fdr20

mean_q3apn_mock_1h mean_q3apn_mock_3h mean_q3apn_mock_24h mean_q3apn_mock_72h
mean_q3apn_01gy_1h mean_q3apn_01gy_3h mean_q3apn_01gy_24h mean_q3apn_01gy_72h 
mean_q3apn_1gy_1h mean_q3apn_1gy_3h mean_q3apn_1gy_24h mean_q3apn_1gy_72h

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

 mean_q3apn_mock_1h=mean_norm_apn_mock_1h
 mean_q3apn_mock_3h=mean_norm_apn_mock_3h
 mean_q3apn_mock_24h=mean_norm_apn_mock_24h
 mean_q3apn_mock_72h=mean_norm_apn_mock_72h
 mean_q3apn_01gy_1h=mean_norm_apn_01gy_1h
 mean_q3apn_01gy_3h=mean_norm_apn_01gy_3h
 mean_q3apn_01gy_24h=mean_norm_apn_01gy_24h
 mean_q3apn_01gy_72h =mean_norm_apn_01gy_72h
 mean_q3apn_1gy_1h=mean_norm_apn_1gy_1h
 mean_q3apn_1gy_3h=mean_norm_apn_1gy_3h
 mean_q3apn_1gy_24h=mean_norm_apn_1gy_24h
 mean_q3apn_1gy_72h=mean_norm_apn_1gy_72h
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
  keep gene_id symbol log2_fc_01gy_v_Mock_lh fdr_01gy_v_Mock_1h flag_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_1h_fdr05
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
  keep gene_id symbol log2_fc_01gy_v_Mock_lh fdr_01gy_v_Mock_1h flag_01gy_v_Mock_1h_fdr10 sign_01gy_v_Mock_1h_fdr10
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
  keep gene_id symbol log2_fc_01gy_v_Mock_lh fdr_01gy_v_Mock_1h flag_01gy_v_Mock_1h_fdr20 sign_01gy_v_Mock_1h_fdr20
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

proc export data=fdr05_full outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes_fdr05_full_results.csv"
dbms=csv replace;
run;
proc export data=fdr05_reduced outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes_fdr05_simple_results.csv"
dbms=csv replace;
run;


proc export data=fdr10_full outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes_fdr10_full_results.csv"
dbms=csv replace;
run;
proc export data=fdr10_reduced outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes_fdr10_simple_results.csv"
dbms=csv replace;
run;


proc export data=fdr20_full outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes_fdr20_full_results.csv"
dbms=csv replace;
run;
proc export data=fdr20_reduced outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes_fdr20_simple_results.csv"
dbms=csv replace;
run;




/* New comparisons
Overlap between Treat vs. Mock for each timepoint*/

/* FDR 5% */
proc freq data=gene_summary3 noprint;
  tables flag_01gy_v_Mock_1h_fdr05*flag_01gy_v_Mock_3h_fdr05*flag_01gy_v_Mock_24h_fdr05*flag_01gy_v_Mock_72h_fdr05 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 flag_01gy_    flag_01gy_    flag_01gy_    flag_01gy_
 v_Mock_1h_    v_Mock_3h_      v_Mock_       v_Mock_
    fdr05         fdr05       24h_fdr05     72h_fdr05    COUNT

      .             .             .             .         7326
      0             0             0             0        16300
      0             0             0             1          987
      0             0             1             0          121
      0             0             1             1           43
      0             1             0             0           50
      0             1             0             1           19
      0             1             1             1            3
      1             0             0             0         3752
      1             0             0             1          285
      1             0             1             0           97
      1             0             1             1           28
      1             1             0             0           30
      1             1             0             1           17
      1             1             1             0            1
      1             1             1             1            5

*/
proc freq data=gene_summary3 noprint;
  tables sign_01gy_v_Mock_1h_fdr05*sign_01gy_v_Mock_3h_fdr05*sign_01gy_v_Mock_24h_fdr05*sign_01gy_v_Mock_72h_fdr05 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 sign_01gy_    sign_01gy_    sign_01gy_    sign_01gy_
 v_Mock_1h_    v_Mock_3h_     v_Mock_       v_Mock_
   fdr05         fdr05       24h_fdr05     72h_fdr05     COUNT

                                                          7326
     D             D             D             N             1
     D             D             N             D             5
     D             D             N             N            14
     D             D             N             U             3
     D             D             U             U             1
     D             N             D             D            12
     D             N             D             N            31
     D             N             N             D           108
     D             N             N             N          2010
     D             N             N             U            34
     D             N             U             N            18
     D             N             U             U             3
     D             U             N             D             2
     D             U             N             N             5
     N             D             D             D             1
     N             D             N             D             4
     N             D             N             N            21
     N             D             N             U            13
     N             D             U             U             1
     N             N             D             D            34
     N             N             D             N            66
     N             N             N             D           662
     N             N             N             N         16300
     N             N             N             U           325
     N             N             U             N            55
     N             N             U             U             9
     N             U             N             N            29
     N             U             N             U             2
     N             U             U             U             1
     U             D             N             D             1
     U             D             N             N             5
     U             D             N             U             5
     U             D             U             U             1
     U             N             D             D             4
     U             N             D             N            12
     U             N             N             D            35
     U             N             N             N          1742
     U             N             N             U           108
     U             N             U             N            36
     U             N             U             U             9
     U             U             N             N             6
     U             U             N             U             1
     U             U             U             U             3
*/

proc freq data=gene_summary3 noprint;
  tables flag_1gy_v_Mock_1h_fdr05*flag_1gy_v_Mock_3h_fdr05*flag_1gy_v_Mock_24h_fdr05*flag_1gy_v_Mock_72h_fdr05 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
  flag_1gy_     flag_1gy_    flag_1gy_    flag_1gy_
 v_Mock_1h_    v_Mock_3h_     v_Mock_      v_Mock_
    fdr05         fdr05      24h_fdr05    72h_fdr05    COUNT

      .             .            .            .         7326
      0             0            0            0        18321
      0             0            0            1           34
      0             0            1            0           26
      0             1            0            0          554
      0             1            0            1            8
      1             0            0            0         2589
      1             0            0            1           12
      1             0            1            0           13
      1             1            0            0          177
      1             1            0            1            3
      1             1            1            0            1
*/
proc freq data=gene_summary3 noprint;
  tables sign_1gy_v_Mock_1h_fdr05*sign_1gy_v_Mock_3h_fdr05*sign_1gy_v_Mock_24h_fdr05*sign_1gy_v_Mock_72h_fdr05 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 sign_1gy_     sign_1gy_     sign_1gy_    sign_1gy_
 v_Mock_1h_    v_Mock_3h_     v_Mock_      v_Mock_
   fdr05         fdr05       24h_fdr05    72h_fdr05    COUNT

                                                        7326
     D             D             N            D            1
     D             D             N            N           57
     D             N             D            N            1
     D             N             N            D            3
     D             N             N            N         1219
     D             N             N            U            1
     D             N             U            N            4
     D             U             N            N           27
     N             D             N            D            6
     N             D             N            N          336
     N             N             D            N           13
     N             N             N            D            5
     N             N             N            N        18321
     N             N             N            U           29
     N             N             U            N           13
     N             U             N            N          218
     N             U             N            U            2
     U             D             N            N           28
     U             D             U            N            1
     U             N             D            N            1
     U             N             N            D            1
     U             N             N            N         1370
     U             N             N            U            7
     U             N             U            N           7
     U             U             N            N          65
     U             U             N            U           2

*/



/* FDR 10% */

proc freq data=gene_summary3 noprint;
  tables flag_01gy_v_Mock_1h_fdr10*flag_01gy_v_Mock_3h_fdr10*flag_01gy_v_Mock_24h_fdr10*flag_01gy_v_Mock_72h_fdr10 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 flag_01gy_    flag_01gy_    flag_01gy_    flag_01gy_
 v_Mock_1h_    v_Mock_3h_      v_Mock_       v_Mock_
    fdr10         fdr10       24h_fdr10     72h_fdr10    COUNT

      .             .             .             .         7326
      0             0             0             0        14450
      0             0             0             1         1379
      0             0             1             0          206
      0             0             1             1           88
      0             1             0             0          126
      0             1             0             1           44
      0             1             1             0            3
      0             1             1             1            3
      1             0             0             0         4466
      1             0             0             1          565
      1             0             1             0          190
      1             0             1             1           78
      1             1             0             0           92
      1             1             0             1           35
      1             1             1             0            1
      1             1             1             1           12

*/
proc freq data=gene_summary3 noprint;
  tables sign_01gy_v_Mock_1h_fdr10*sign_01gy_v_Mock_3h_fdr10*sign_01gy_v_Mock_24h_fdr10*sign_01gy_v_Mock_72h_fdr10 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
  sign_01gy_    sign_01gy_    sign_01gy_    sign_01gy_
  v_Mock_1h_    v_Mock_3h_     v_Mock_       v_Mock_
    fdr10         fdr10       24h_fdr10     72h_fdr10     COUNT

                                                           7326
      D             D             D             D             3
      D             D             N             D            10
      D             D             N             N            38
      D             D             N             U             5
      D             D             U             U             1
      D             N             D             D            32
      D             N             D             N            70
      D             N             D             U             1
      D             N             N             D           207
      D             N             N             N          2404
      D             N             N             U            78
      D             N             U             N            35
      D             N             U             U             7
      D             U             D             N             1
      D             U             N             D             2
      D             U             N             N            22
      D             U             N             U             3
      N             D             N             D            10
      N             D             N             N            52
      N             D             N             U            30
      N             D             U             U             2
      N             N             D             D            63
      N             N             D             N           105
      N             N             N             D           896
      N             N             N             N         14450
      N             N             N             U           483
      N             N             U             D             1
      N             N             U             N           101
      N             N             U             U            24
      N             U             D             N             2
      N             U             N             D             1
      N             U             N             N            74
      N             U             N             U             3
      N             U             U             N             1
      N             U             U             U             1
      U             D             D             D             2
      U             D             N             D             1
      U             D             N             N            13
      U             D             N             U            11
      U             D             U             U             2
      U             N             D             D             6
      U             N             D             N            18
      U             N             N             D            78
      U             N             N             N          2062
      U             N             N             U           202
      U             N             U             D             1
      U             N             U             N            67
      U             N             U             U            31
      U             U             N             N            19
      U             U             N             U             3
      U             U             U             U             4
*/

proc freq data=gene_summary3 noprint;
  tables flag_1gy_v_Mock_1h_fdr10*flag_1gy_v_Mock_3h_fdr10*flag_1gy_v_Mock_24h_fdr10*flag_1gy_v_Mock_72h_fdr10 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 flag_1gy_     flag_1gy_    flag_1gy_    flag_1gy_
v_Mock_1h_    v_Mock_3h_     v_Mock_      v_Mock_
   fdr10         fdr10      24h_fdr10    72h_fdr10    COUNT

     .             .            .            .         7326
     0             0            0            0        16625
     0             0            0            1          163
     0             0            1            0           99
     0             0            1            1            5
     0             1            0            0         1064
     0             1            0            1           47
     0             1            1            0           14
     0             1            1            1            2
     1             0            0            0         3185
     1             0            0            1           56
     1             0            1            0           59
     1             0            1            1            2
     1             1            0            0          391
     1             1            0            1           20
     1             1            1            0            6
*/

proc freq data=gene_summary3 noprint;
  tables sign_1gy_v_Mock_1h_fdr10*sign_1gy_v_Mock_3h_fdr10*sign_1gy_v_Mock_24h_fdr10*sign_1gy_v_Mock_72h_fdr10 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 sign_1gy_     sign_1gy_     sign_1gy_    sign_1gy_
 v_Mock_1h_    v_Mock_3h_     v_Mock_      v_Mock_
   fdr10         fdr10       24h_fdr10    72h_fdr10    COUNT

                                                        7326
     D             D             N            D            5
     D             D             N            N          126
     D             N             D            D            1
     D             N             D            N            6
     D             N             N            D           17
     D             N             N            N         1519
     D             N             N            U            8
     D             N             U            N           14
     D             U             D            N            1
     D             U             N            D            1
     D             U             N            N           63
     D             U             N            U            2
     N             D             D            N            9
     N             D             N            D           23
     N             D             N            N          612
     N             D             N            U            1
     N             D             U            U            1
     N             N             D            D            1
     N             N             D            N           59
     N             N             N            D           36
     N             N             N            N        16625
     N             N             N            U          127
     N             N             U            N           40
     N             N             U            U            4
     N             U             D            D            1
     N             U             D            N            3
     N             U             N            D            1
     N             U             N            N          452
     N             U             N            U           22
     N             U             U            N            2
     U             D             D            N            1
     U             D             N            D            3
     U             D             N            N           55
     U             D             N            U            1
     U             D             U            N            1
     U             N             D            N            6
     U             N             N            D            2
     U             N             N            N         1666
     U             N             N            U           29
     U             N             U            N           33
     U             N             U            U            1
     U             U             N            N          147
     U             U             N            U            8
     U             U             U            N            3


*/


/* FDR 20% */
proc freq data=gene_summary3 noprint;
  tables flag_01gy_v_Mock_1h_fdr20*flag_01gy_v_Mock_3h_fdr20*flag_01gy_v_Mock_24h_fdr20*flag_01gy_v_Mock_72h_fdr20 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
 flag_01gy_    flag_01gy_    flag_01gy_    flag_01gy_
 v_Mock_1h_    v_Mock_3h_      v_Mock_       v_Mock_
    fdr20         fdr20       24h_fdr20     72h_fdr20    COUNT

      .             .             .             .         7326
      0             0             0             0        11545
      0             0             0             1         1674
      0             0             1             0          451
      0             0             1             1          254
      0             1             0             0          341
      0             1             0             1          118
      0             1             1             0           16
      0             1             1             1           11
      1             0             0             0         5159
      1             0             0             1         1091
      1             0             1             0          387
      1             0             1             1          234
      1             1             0             0          297
      1             1             0             1          111
      1             1             1             0           23
      1             1             1             1           26

*/

proc freq data=gene_summary3 noprint;
  tables sign_01gy_v_Mock_1h_fdr20*sign_01gy_v_Mock_3h_fdr20*sign_01gy_v_Mock_24h_fdr20*sign_01gy_v_Mock_72h_fdr20 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
sign_01gy_    sign_01gy_    sign_01gy_    sign_01gy_
v_Mock_1h_    v_Mock_3h_     v_Mock_       v_Mock_
  fdr20         fdr20       24h_fdr20     72h_fdr20     COUNT

                                                         7326
    D             D             D             D             4
    D             D             D             N             1
    D             D             N             D            45
    D             D             N             N           112
    D             D             N             U             9
    D             D             U             N             1
    D             D             U             U             3
    D             N             D             D           101
    D             N             D             N           127
    D             N             D             U             5
    D             N             N             D           354
    D             N             N             N          2523
    D             N             N             U           170
    D             N             U             D             2
    D             N             U             N            70
    D             N             U             U            25
    D             U             D             D             2
    D             U             D             N             9
    D             U             D             U             1
    D             U             N             D             1
    D             U             N             N            62
    D             U             N             U            10
    D             U             U             N             2
    N             D             D             D             4
    N             D             D             N             2
    N             D             N             D            32
    N             D             N             N           148
    N             D             N             U            74
    N             D             U             N             4
    N             D             U             U             6
    N             N             D             D           201
    N             N             D             N           203
    N             N             D             U             4
    N             N             N             D          1103
    N             N             N             N         12016
    N             N             N             U           673
    N             N             U             D             3
    N             N             U             N           272
    N             N             U             U            62
    N             U             D             N             4
    N             U             N             D             2
    N             U             N             N           212
    N             U             N             U            14
    N             U             U             N             7
    N             U             U             U             1
    U             D             D             D             3
    U             D             D             N             4
    U             D             N             D             4
    U             D             N             N            37
    U             D             N             U            28
    U             D             U             N             5
    U             D             U             U             5
    U             N             D             D            16
    U             N             D             N            25
    U             N             D             U             2
    U             N             N             D           134
    U             N             N             N          2165
    U             N             N             U           331
    U             N             U             D             1
    U             N             U             N           141
    U             N             U             U            66
    U             U             N             N            67
    U             U             N             U            10
    U             U             U             U             8


*/

proc freq data=gene_summary3 noprint;
  tables flag_1gy_v_Mock_1h_fdr20*flag_1gy_v_Mock_3h_fdr20*flag_1gy_v_Mock_24h_fdr20*flag_1gy_v_Mock_72h_fdr20 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
  flag_1gy_     flag_1gy_    flag_1gy_    flag_1gy_
 v_Mock_1h_    v_Mock_3h_     v_Mock_      v_Mock_
    fdr20         fdr20      24h_fdr20    72h_fdr20    COUNT

      .             .            .            .         7326
      0             0            0            0        13406
      0             0            0            1          704
      0             0            1            0          347
      0             0            1            1           34
      0             1            0            0         1739
      0             1            0            1          213
      0             1            1            0           86
      0             1            1            1           19
      1             0            0            0         3702
      1             0            0            1          248
      1             0            1            0          176
      1             0            1            1           32
      1             1            0            0          861
      1             1            0            1          117
      1             1            1            0           45
      1             1            1            1            9

*/

proc freq data=gene_summary3 noprint;
  tables sign_1gy_v_Mock_1h_fdr20*sign_1gy_v_Mock_3h_fdr20*sign_1gy_v_Mock_24h_fdr20*sign_1gy_v_Mock_72h_fdr20 / out=diff_count;
run;

proc print data=diff_count;
run;
/*
sign_1gy_     sign_1gy_     sign_1gy_    sign_1gy_
v_Mock_1h_    v_Mock_3h_     v_Mock_      v_Mock_
  fdr20         fdr20       24h_fdr20    72h_fdr20    COUNT

                                                       7326
    D             D             D            D            3
    D             D             D            N            3
    D             D             N            D           40
    D             D             N            N          268
    D             D             U            N           12
    D             N             D            D           10
    D             N             D            N           26
    D             N             N            D           76
    D             N             N            N         1813
    D             N             N            U           38
    D             N             U            N           51
    D             N             U            U            4
    D             U             D            N            4
    D             U             N            D            5
    D             U             N            N          154
    D             U             N            U           11
    D             U             U            N            6
    D             U             U            U            2
    N             D             D            D            9
    N             D             D            N           42
    N             D             D            U            1
    N             D             N            D          100
    N             D             N            N          891
    N             D             N            U            4
    N             D             U            N           20
    N             D             U            U            1
    N             N             D            D           25
    N             N             D            N          208
    N             N             N            D          159
    N             N             N            N        13406
    N             N             N            U          545
    N             N             U            N          139
    N             N             U            U            9
    N             U             D            D            2
    N             U             D            N            9
    N             U             N            D            4
    N             U             N            N          848
    N             U             N            U          105
    N             U             U            D            1
    N             U             U            N           15
    N             U             U            U            5
    U             D             D            N            2
    U             D             N            D            6
    U             D             N            N          129
    U             D             N            U            2
    U             D             U            N            6
    U             N             D            D            4
    U             N             D            N           15
    U             N             N            D           18
    U             N             N            N         1889
    U             N             N            U          116
    U             N             U            D            1
    U             N             U            N           84
    U             N             U            U           13
    U             U             D            N            2
    U             U             N            D            3
    U             U             N            N          310
    U             U             N            U           50
    U             U             U            N           10
    U             U             U            U            4
*/


