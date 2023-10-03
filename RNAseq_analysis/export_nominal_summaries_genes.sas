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


/* Get P-values for the following:
   treatment   time trt*time

 */

data main_p;
  set rs.arab_anova_main_trt_tm_genes;
  length effect_p $30.;
   format ProbF best32. ;
  if effect="treatment" then effect_p="P_treatment";
  if effect="time" then effect_p="P_time";
  if effect="time*treatment" then effect_p="P_trt_by_time";
  keep gene_id effect_p probf;
run;

proc sort data=main_p;
   by gene_id effect_p;
proc transpose data=main_p out=main_p_sbys;
  by gene_id;
  id effect_p;
  var probf;
run;

data main_p_sbys2;
   set main_p_sbys;
   drop _NAME_ _LABEL_;
run;

data contrast_p;
  set rs.arab_gene_cntrs_constr;
  length p_label $100.;
   format ProbF best32. ;

  if label="0.1gy-Mock: 1h" then p_label="p_01gy_v_Mock_1h";
  if label= "0.1gy-Mock: 3h" then p_label="p_01gy_v_Mock_3h";
  if label= "0.1gy-Mock: 24h" then p_label="p_01gy_v_Mock_24h";
  if label= "0.1gy-Mock: 72h" then p_label="p_01gy_v_Mock_72h";

  if label= "1gy-Mock: 1h" then p_label="p_1gy_v_Mock_1h";
  if label= "1gy-Mock: 3h" then p_label="p_1gy_v_Mock_3h";
  if label= "1gy-Mock: 24h" then p_label="p_1gy_v_Mock_24h";
  if label= "1gy-Mock: 72h" then p_label="p_1gy_v_Mock_72h";

  if label="1gy-0.1gy: 1h" then p_label="p_1gy_v_01gy_1h";
  if label="1gy-0.1gy: 3h" then p_label="p_1gy_v_01gy_3h";
  if label="1gy-0.1gy: 24h" then p_label="p_1gy_v_01gy_24h";
  if label="1gy-0.1gy: 72h" then p_label="p_1gy_v_01gy_72h";

  if label="Mock: 3h-1h" then p_label="p_mock_3h_v_1h";
  if label="Mock: 24h-1h" then p_label="p_mock_24h_v_1h";
  if label="Mock: 72h-1h" then p_label="p_mock_72h_v_1h";
  if label="Mock: 24h-3h" then p_label="p_mock_24h_v_3h";
  if label="Mock: 72h-3h" then p_label="p_mock_72h_v_3h";
  if label="Mock: 72h-24h" then p_label="p_mock_72h_v_24h";

  if label="0.1gy: 3h-1h" then p_label="p_01gy_3h_v_1h";
  if label="0.1gy: 24h-1h" then p_label="p_01gy_24h_v_1h";
  if label="0.1gy: 72h-1h" then p_label="p_01gy_72h_v_1h";
  if label="0.1gy: 24h-3h" then p_label="p_01gy_24h_v_3h";
  if label="0.1gy: 72h-3h" then p_label="p_01gy_72h_v_3h";
  if label="0.1gy: 72h-24h" then p_label="p_01gy_72h_v_24h";

  if label="1gy: 3h-1h" then p_label="p_1gy_3h_v_1h";
  if label="1gy: 24h-1h" then p_label="p_1gy_24h_v_1h";
  if label="1gy: 72h-1h" then p_label="p_1gy_72h_v_1h";
  if label="1gy: 24h-3h" then p_label="p_1gy_24h_v_3h";
  if label="1gy: 72h-3h" then p_label="p_1gy_72h_v_3h";
  if label="1gy: 72h-24h" then p_label="p_1gy_72h_v_24h";

  if label="0.1gy 3h-1h = Mock 3h-1h" then p_label="p_01gy_v_Mock_3h_sub_1h";
  if label="0.1gy 24h-1h = Mock 24h-1h" then p_label="p_01gy_v_Mock_24h_sub_1h";
  if label="0.1gy 72h-1h = Mock 72h-1h" then p_label="p_01gy_v_Mock_72h_sub_1h";
  if label="1gy 3h-1h = Mock 3h-1h" then p_label="p_1gy_v_Mock_3h_sub_1h";
  if label="1gy 24h-1h = Mock 24h-1h" then p_label="p_1gy_v_Mock_24h_sub_1h";
  if label="1gy 72h-1h = Mock 72h-1h" then p_label="p_1gy_v_Mock_72h_sub_1h";

  keep gene_id p_label probf;
run;

proc sort data=contrast_p;
   by gene_id p_label;
proc transpose data=contrast_p out=contrast_p_sbys;
   by gene_id;
   id p_label;
   var probf;
run;


data contrast_p_sbys2;
   set contrast_p_sbys;
   drop _NAME_ _LABEL_;
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
  if label= "0.1gy-Mock: 3h" then est_label="Estimate_01gy_v_Mock_3h";
  if label= "0.1gy-Mock: 24h" then est_label="Estimate_01gy_v_Mock_24h";
  if label= "0.1gy-Mock: 72h" then est_label="Estimate_01gy_v_Mock_72h";

  if label= "1gy-Mock: 1h" then est_label="Estimate_1gy_v_Mock_1h";
  if label= "1gy-Mock: 3h" then est_label="Estimate_1gy_v_Mock_3h";
  if label= "1gy-Mock: 24h" then est_label="Estimate_1gy_v_Mock_24h";
  if label= "1gy-Mock: 72h" then est_label="Estimate_1gy_v_Mock_72h";

  if label="1gy-0.1gy: 1h" then est_label="Estimate_1gy_v_01gy_1h";
  if label="1gy-0.1gy: 3h" then est_label="Estimate_1gy_v_01gy_3h";
  if label="1gy-0.1gy: 24h" then est_label="Estimate_1gy_v_01gy_24h";
  if label="1gy-0.1gy: 72h" then est_label="Estimate_1gy_v_01gy_72h";

  if label="Mock: 3h-1h" then est_label="Estimate_mock_3h_v_1h";
  if label="Mock: 24h-1h" then est_label="Estimate_mock_24h_v_1h";
  if label="Mock: 72h-1h" then est_label="Estimate_mock_72h_v_1h";
  if label="Mock: 24h-3h" then est_label="Estimate_mock_24h_v_3h";
  if label="Mock: 72h-3h" then est_label="Estimate_mock_72h_v_3h";
  if label="Mock: 72h-24h" then est_label="Estimate_mock_72h_v_24h";

  if label="0.1gy: 3h-1h" then est_label="Estimate_01gy_3h_v_1h";
  if label="0.1gy: 24h-1h" then est_label="Estimate_01gy_24h_v_1h";
  if label="0.1gy: 72h-1h" then est_label="Estimate_01gy_72h_v_1h";
  if label="0.1gy: 24h-3h" then est_label="Estimate_01gy_24h_v_3h";
  if label="0.1gy: 72h-3h" then est_label="Estimate_01gy_72h_v_3h";
  if label="0.1gy: 72h-24h" then est_label="Estimate_01gy_72h_v_24h";

  if label="1gy: 3h-1h" then est_label="Estimate_1gy_3h_v_1h";
  if label="1gy: 24h-1h" then est_label="Estimate_1gy_24h_v_1h";
  if label="1gy: 72h-1h" then est_label="Estimate_1gy_72h_v_1h";
  if label="1gy: 24h-3h" then est_label="Estimate_1gy_24h_v_3h";
  if label="1gy: 72h-3h" then est_label="Estimate_1gy_72h_v_3h";
  if label="1gy: 72h-24h" then est_label="Estimate_1gy_72h_v_24h";

  if label="0.1gy 3h-1h = Mock 3h-1h" then est_label="Estimate_01gy_v_Mock_3h_sub_1h";
  if label="0.1gy 24h-1h = Mock 24h-1h" then est_label="Estimate_01gy_v_Mock_24h_sub_1h";
  if label="0.1gy 72h-1h = Mock 72h-1h" then est_label="Estimate_01gy_v_Mock_72h_sub_1h";
  if label="1gy 3h-1h = Mock 3h-1h" then est_label="Estimate_1gy_v_Mock_3h_sub_1h";
  if label="1gy 24h-1h = Mock 24h-1h" then est_label="Estimate_1gy_v_Mock_24h_sub_1h";
  if label="1gy 72h-1h = Mock 72h-1h" then est_label="Estimate_1gy_v_Mock_72h_sub_1h";

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


/* Get Up/Down signs */

data sign;
  set rs.arab_sign_by_contrast_gene;
run;


proc sort data=gene_coord2;
   by gene_id;
proc sort data=main_p_sbys2;
   by gene_id;
proc sort data=contrast_p_sbys2;
   by gene_id;
proc sort data=lsmeans_sbys2;
   by gene_id;
proc sort data=estimates_sbys2;
   by gene_id;
proc sort data=sign;
   by gene_id;
proc sort data=flag_dtct;
   by gene_id;
run;


data gene_summary gene_skipped;
   length gene_id $444.;
   format gene_id $444.;
   informat gene_id $444.;
   merge gene_coord2 (in=in1) flag_dtct (in=in2) main_p_sbys2 contrast_p_sbys2 lsmeans_sbys2 estimates_sbys2 sign;
   by gene_id;
   if in1 and in2 then output gene_summary;
   else if in1 then output gene_skipped;
run;

/* Reorder variables */

proc contents data=gene_summary;
run;


data gene_summary2;
 retain gene_id symbol chrom gene_start gene_stop strand flag_Mock_on flag_01gy_on flag_1gy_on flag_1hr_on flag_3hr_on
         flag_24hr_on flag_72gy_on flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
         flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on flag_1gy_1hr_on flag_1gy_3hr_on
         flag_1gy_24hr_on flag_1gy_72hr_on LSmeans_Mock_1h LSmeans_Mock_3h LSmeans_Mock_24 LSmeans_Mock_72
         LSmeans_01gy_1h LSmeans_01gy_3h LSmeans_01gy_24 LSmeans_01gy_72 LSmeans_1gy_1h LSmeans_1gy_3h
         LSmeans_1gy_24h LSmeans_1gy_72h P_time P_treatment P_trt_by_time p_01gy_v_Mock_1h p_01gy_v_Mock_3h
         p_01gy_v_Mock_24h  p_01gy_v_Mock_72h p_1gy_v_Mock_1h p_1gy_v_Mock_3h p_1gy_v_Mock_24h p_1gy_v_Mock_72h
         p_1gy_v_01gy_1h p_1gy_v_01gy_3h p_1gy_v_01gy_24h p_1gy_v_01gy_72h p_mock_3h_v_1h p_mock_24h_v_1h
         p_mock_72h_v_1h p_mock_24h_v_3h p_mock_72h_v_3h p_mock_72h_v_24h p_01gy_3h_v_1h p_01gy_24h_v_1h
         p_01gy_72h_v_1h p_01gy_24h_v_3h p_01gy_72h_v_3h p_01gy_72h_v_24h p_1gy_3h_v_1h p_1gy_24h_v_1h p_1gy_72h_v_1h
         p_1gy_24h_v_3h p_1gy_72h_v_3h p_1gy_72h_v_24h

         p_01gy_v_Mock_3h_sub_1h
p_01gy_v_Mock_24h_sub_1h
p_01gy_v_Mock_72h_sub_1h
p_1gy_v_Mock_3h_sub_1h
p_1gy_v_Mock_24h_sub_1h
p_1gy_v_Mock_72h_sub_1h

         Estimate_01gy_v_Mock_1h  Estimate_01gy_v_Mock_3h
         Estimate_01gy_v_Mock_24h Estimate_01gy_v_Mock_72h Estimate_1gy_v_Mock_1h Estimate_1gy_v_Mock_3h
         Estimate_1gy_v_Mock_24h Estimate_1gy_v_Mock_72h Estimate_1gy_v_01gy_1h Estimate_1gy_v_01gy_3h
         Estimate_1gy_v_01gy_24h Estimate_1gy_v_01gy_72h Estimate_mock_3h_v_1h Estimate_mock_24h_v_1h
         Estimate_mock_72h_v_1h Estimate_mock_24h_v_3h Estimate_mock_72h_v_3h Estimate_mock_72h_v_24h
         Estimate_01gy_3h_v_1h Estimate_01gy_24h_v_1h Estimate_01gy_72h_v_1h Estimate_01gy_24h_v_3h
         Estimate_01gy_72h_v_3h Estimate_01gy_72h_v_24h Estimate_1gy_3h_v_1h Estimate_1gy_24h_v_1h
         Estimate_1gy_72h_v_1h Estimate_1gy_24h_v_3h Estimate_1gy_72h_v_3h Estimate_1gy_72h_v_24h

Estimate_01gy_v_Mock_3h_sub_1h
Estimate_01gy_v_Mock_24h_sub_1h
Estimate_01gy_v_Mock_72h_sub_1h
Estimate_1gy_v_Mock_3h_sub_1h
Estimate_1gy_v_Mock_24h_sub_1h
Estimate_1gy_v_Mock_72h_sub_1h


         sign_01gy_v_Mock_1h  sign_01gy_v_Mock_3h
         sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h
         sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h sign_1gy_v_01gy_1h sign_1gy_v_01gy_3h
         sign_1gy_v_01gy_24h sign_1gy_v_01gy_72h sign_mock_3h_v_1h sign_mock_24h_v_1h
         sign_mock_72h_v_1h sign_mock_24h_v_3h sign_mock_72h_v_3h sign_mock_72h_v_24h
         sign_01gy_3h_v_1h sign_01gy_24h_v_1h sign_01gy_72h_v_1h sign_01gy_24h_v_3h
         sign_01gy_72h_v_3h sign_01gy_72h_v_24h sign_1gy_3h_v_1h sign_1gy_24h_v_1h
         sign_1gy_72h_v_1h sign_1gy_24h_v_3h sign_1gy_72h_v_3h sign_1gy_72h_v_24h

sign_01gy_v_Mock_3h_sub_1h
sign_01gy_v_Mock_24h_sub_1h
sign_01gy_v_Mock_72h_sub_1h
sign_1gy_v_Mock_3h_sub_1h
sign_1gy_v_Mock_24h_sub_1h
sign_1gy_v_Mock_72h_sub_1h

;
  


         set gene_summary;


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

  drop _TYPE_ _FREQ_;
run;

/* Export */

proc export data=gene_summary2 outfile="!PATCON/arabidopsis/analysis_output/arab_rs_genes_tair10_all_genes.csv"
dbms=csv replace;
run;


/* Need to count number of sig fusions by effect/contrast, and number on per trt, time and trt*time */

proc freq data=gene_summary2 noprint;
   tables flag_Mock_on*flag_01gy_on*flag_1gy_on / out=count_on_by_trt;
run;

proc freq data=gene_summary2 noprint;
   tables flag_1hr_on*flag_3hr_on*flag_24hr_on*flag_72gy_on / out=count_on_by_time;
run;

proc freq data=gene_summary2 noprint;
   tables flag_Mock_1hr_on*flag_Mock_3hr_on*flag_Mock_24hr_on*flag_Mock_72hr_on*
          flag_01gy_1hr_on*flag_01gy_3hr_on*flag_01gy_24hr_on*flag_01gy_72hr_on*
          flag_1gy_1hr_on*flag_1gy_3hr_on*flag_1gy_24hr_on*flag_1gy_72hr_on / out=count_on_by_trt_time;
run;

proc export data=count_on_by_trt outfile="!PATCON/arabidopsis/analysis_output/genes_on_by_trt_apn0.csv"
    dbms=csv replace;
run;

proc export data=count_on_by_time outfile="!PATCON/arabidopsis/analysis_output/genes_on_by_time_apn0.csv"
    dbms=csv replace;
run;

proc export data=count_on_by_trt_time outfile="!PATCON/arabidopsis/analysis_output/genes_on_by_trt_time_apn0.csv"
    dbms=csv replace;
run;


/* Now flag and count P<0.05 for each main effect and contrast */

%macro flagp05(pvalue,flagname);
data flag_p05;
  set gene_summary2;
  if &pvalue. = . then &flagname.=.;
  else if &pvalue. <0.05 then &flagname.=1;
  else &flagname.=0;
run;

proc freq data=flag_p05;
  tables &flagname.;
run;
%mend;


%flagp05(P_time,flag_time_p05);
%flagp05(P_treatment,flag_trt_p05);
%flagp05(P_trt_by_time,flag_trt_time_p05);
%flagp05(p_01gy_v_Mock_1h,flag_01vMock_1_p05);
%flagp05(p_01gy_v_Mock_3h,flag_01vMock_3_p05);
%flagp05(p_01gy_v_Mock_24h,flag_01vMock_24_p05);
%flagp05(p_01gy_v_Mock_72h,flag_01vMock_72_p05);
%flagp05(p_1gy_v_Mock_1h,flag_1vMock_1_p05);
%flagp05(p_1gy_v_Mock_3h,flag_1vMock_3_p05);
%flagp05(p_1gy_v_Mock_24h,flag_1vMock_24_p05); 
%flagp05(p_1gy_v_Mock_72h,flag_1vMock_72_p05);
%flagp05(p_1gy_v_01gy_1h,flag_1v01gy_1_p05);
%flagp05(p_1gy_v_01gy_3h,flag_1v01gy_3_p05);
%flagp05(p_1gy_v_01gy_24h,flag_1v01gy_24_p05);
%flagp05(p_1gy_v_01gy_72h,flag_1v01gy_72_p05);
%flagp05(p_mock_3h_v_1h,flag_mock_3v1_p05);
%flagp05(p_mock_24h_v_1h,flag_mock_24v1_p05);
%flagp05(p_mock_72h_v_1h,flag_mock_72v1_p05);
%flagp05(p_mock_24h_v_3h,flag_mock_24v3_p05);
%flagp05(p_mock_72h_v_3h,flag_mock_72v3_p05);
%flagp05(p_mock_72h_v_24h,flag_mock_72v24_p05);
%flagp05(p_01gy_3h_v_1h,flag_01gy_3v1_p05);
%flagp05(p_01gy_24h_v_1h,flag_01gy_24v1_p05);
%flagp05(p_01gy_72h_v_1h,flag_01gy_72v1_p05);
%flagp05(p_01gy_24h_v_3h,flag_01gy_24v3_p05);
%flagp05(p_01gy_72h_v_3h,flag_01gy_72v3_p05);
%flagp05(p_01gy_72h_v_24h,flag_01gy_72v24_p05);
%flagp05(p_1gy_3h_v_1h,flag_1gy_3v1_p05);
%flagp05(p_1gy_24h_v_1h,flag_1gy_24v1_p05);
%flagp05(p_1gy_72h_v_1h,flag_1gy_72v1_p05);
%flagp05(p_1gy_24h_v_3h,flag_1gy_24v3_p05);
%flagp05(p_1gy_72h_v_3h,flag_1gy_72v3_p05);
%flagp05(p_1gy_72h_v_24h,flag_1gy_72v24_p05);
%flagp05(p_01gy_v_Mock_3h_sub_1h,flag_01gy_v_Mock_3h_sub_1h_p05);
%flagp05(p_01gy_v_Mock_24h_sub_1h,flag_01gy_v_Mock_24h_sub_1h_p05);
%flagp05(p_01gy_v_Mock_72h_sub_1h,flag_01gy_v_Mock_72h_sub_1h_p05);
%flagp05(p_1gy_v_Mock_3h_sub_1h,flag_1gy_v_Mock_3h_sub_1h_p05);
%flagp05(p_1gy_v_Mock_24h_sub_1h,flag_1gy_v_Mock_24h_sub_1h_p05);
%flagp05(p_1gy_v_Mock_72h_sub_1h,flag_1gy_v_Mock_72h_sub_1h_p05);


/*
                                           Cumulative    Cumulative
flag_time_p05    Frequency     Percent     Frequency      Percent
------------------------------------------------------------------
            0        4478       20.60          4478        20.60
            1       17260       79.40         21738       100.00

                     Frequency Missing = 7326

                                          Cumulative    Cumulative
 flag_trt_p05    Frequency     Percent     Frequency      Percent
 -----------------------------------------------------------------
            0       14714       67.69         14714        67.69
            1        7024       32.31         21738       100.00

                      Frequency Missing = 7326

      flag_trt_                             Cumulative    Cumulative
       time_p05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       15574       71.64         15574        71.64
              1        6164       28.36         21738       100.00

                      Frequency Missing = 7326

    flag_01vMock_                             Cumulative    Cumulative
            1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       15047       69.22         15047        69.22
                1        6691       30.78         21738       100.00

                        Frequency Missing = 7326

    flag_01vMock_                             Cumulative    Cumulative
            3_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18715       86.09         18715        86.09
                1        3023       13.91         21738       100.00

                        Frequency Missing = 7326
     flag_01vMock_                             Cumulative    Cumulative
            24_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       18698       86.02         18698        86.02
                 1        3040       13.98         21738       100.00

                         Frequency Missing = 7326
     flag_01vMock_                             Cumulative    Cumulative
            72_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       17446       80.26         17446        80.26
                 1        4292       19.74         21738       100.00

                         Frequency Missing = 7326

     flag_1vMock_                             Cumulative    Cumulative
            1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       16467       75.75         16467        75.75
                1        5271       24.25         21738       100.00

    flag_1vMock_                             Cumulative    Cumulative
           3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       17651       81.20         17651        81.20
               1        4087       18.80         21738       100.00

     flag_1vMock_                             Cumulative    Cumulative
           24_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       19084       87.79         19084        87.79
                1        2654       12.21         21738       100.00

    flag_1vMock_                             Cumulative    Cumulative
          72_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       18424       84.75         18424        84.75
               1        3314       15.25         21738       100.00


      flag_1v01gy_                             Cumulative    Cumulative
             1_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       19784       91.01         19784        91.01
                 1        1954        8.99         21738       100.00

      flag_1v01gy_                             Cumulative    Cumulative
             3_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       17916       82.42         17916        82.42
                 1        3822       17.58         21738       100.00

    flag_1v01gy_                             Cumulative    Cumulative
          24_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       19065       87.70         19065        87.70
               1        2673       12.30         21738       100.00

     flag_1v01gy_                             Cumulative    Cumulative
           72_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       19792       91.05         19792        91.05
                1        1946        8.95         21738       100.00

        flag_mock_                             Cumulative    Cumulative
           3v1_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       12527       57.63         12527        57.63
                 1        9211       42.37         21738       100.00

         flag_mock_                             Cumulative    Cumulative
           24v1_p05    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0       13290       61.14         13290        61.14
                  1        8448       38.86         21738       100.00

        flag_mock_                             Cumulative    Cumulative
          72v1_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       13987       64.34         13987        64.34
                 1        7751       35.66         21738       100.00


      flag_mock_                             Cumulative    Cumulative
        24v3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       10865       49.98         10865        49.98
               1       10873       50.02         21738       100.00

      flag_mock_                             Cumulative    Cumulative
        72v3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       15799       72.68         15799        72.68
               1        5939       27.32         21738       100.00

      flag_mock_                             Cumulative    Cumulative
       72v24_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       12607       58.00         12607        58.00
               1        9131       42.00         21738       100.00

        flag_01gy_                             Cumulative    Cumulative
           3v1_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       10285       47.31         10285        47.31
                 1       11453       52.69         21738       100.00

      flag_01gy_                             Cumulative    Cumulative
        24v1_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       10567       48.61         10567        48.61
               1       11171       51.39         21738       100.00


      flag_01gy_                             Cumulative    Cumulative
        72v1_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       11308       52.02         11308        52.02
               1       10430       47.98         21738       100.00

    flag_01gy_                             Cumulative    Cumulative
      24v3_p05    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0       10364       47.68         10364        47.68
             1       11374       52.32         21738       100.00

      flag_01gy_                             Cumulative    Cumulative
        72v3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       12855       59.14         12855        59.14
               1        8883       40.86         21738       100.00

     flag_01gy_                             Cumulative    Cumulative
      72v24_p05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       12522       57.60         12522        57.60
              1        9216       42.40         21738       100.00

                                                  Cumulative    Cumulative
     flag_1gy_3v1_p05    Frequency     Percent     Frequency      Percent
     ---------------------------------------------------------------------
                    0       10983       50.52         10983        50.52
                    1       10755       49.48         21738       100.00

       flag_1gy_                             Cumulative    Cumulative
        24v1_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       11431       52.59         11431        52.59
               1       10307       47.41         21738       100.00

        flag_1gy_                             Cumulative    Cumulative
         72v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       12259       56.39         12259        56.39
                1        9479       43.61         21738       100.00

       flag_1gy_                             Cumulative    Cumulative
        24v3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       10675       49.11         10675        49.11
               1       11063       50.89         21738       100.00

      flag_1gy_                             Cumulative    Cumulative
       72v3_p05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       14954       68.79         14954        68.79
              1        6784       31.21         21738       100.00


         flag_1gy_                             Cumulative    Cumulative
         72v24_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       12518       57.59         12518        57.59
                 1        9220       42.41         21738       100.00

    flag_01gy_v_
    Mock_3h_sub_                             Cumulative    Cumulative
          1h_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       16749       77.05         16749        77.05
               1        4989       22.95         21738       100.00

                       Frequency Missing = 7326

      flag_01gy_v_
     Mock_24h_sub_                             Cumulative    Cumulative
            1h_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       17130       78.80         17130        78.80
                 1        4608       21.20         21738       100.00

      flag_01gy_v_
     Mock_72h_sub_                             Cumulative    Cumulative
            1h_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       16518       75.99         16518        75.99
                 1        5220       24.01         21738       100.00

 flag_1gy_v_Mock_                             Cumulative    Cumulative
    3h_sub_1h_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       17653       81.21         17653        81.21
                1        4085       18.79         21738       100.00

flag_1gy_v_Mock_                             Cumulative    Cumulative
  24h_sub_1h_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       18056       83.06         18056        83.06
               1        3682       16.94         21738       100.00

 flag_1gy_v_Mock_                             Cumulative    Cumulative
   72h_sub_1h_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       17100       78.66         17100        78.66
                1        4638       21.34         21738       100.00
*/


/* Check: X-1 Y vs Mock: 0.1x1 for each time point, then within treatment */

data flag_p05_diffs;
   set gene_summary2;
   if p_01gy_v_Mock_3h_sub_1h = . then flag_01gy_3h_sub_1h=. ;
   else if p_01gy_v_Mock_3h_sub_1h < 0.05 then flag_01gy_3h_sub_1h=1 ;
   else  flag_01gy_3h_sub_1h=0 ;

   if p_01gy_v_Mock_24h_sub_1h = . then flag_01gy_24h_sub_1h=. ;
   else if p_01gy_v_Mock_24h_sub_1h < 0.05 then flag_01gy_24h_sub_1h=1 ;
   else  flag_01gy_24h_sub_1h=0 ;

   if p_01gy_v_Mock_72h_sub_1h = . then flag_01gy_72h_sub_1h=. ;
   else if p_01gy_v_Mock_72h_sub_1h < 0.05 then flag_01gy_72h_sub_1h=1 ;
   else  flag_01gy_72h_sub_1h=0 ;

   if p_1gy_v_Mock_3h_sub_1h = . then flag_1gy_3h_sub_1h=. ;
   else if p_1gy_v_Mock_3h_sub_1h < 0.05 then flag_1gy_3h_sub_1h=1 ;
   else  flag_1gy_3h_sub_1h=0 ;

   if p_1gy_v_Mock_24h_sub_1h = . then flag_1gy_24h_sub_1h=. ;
   else if p_1gy_v_Mock_24h_sub_1h < 0.05 then flag_1gy_24h_sub_1h=1 ;
   else  flag_1gy_24h_sub_1h=0 ;

   if p_1gy_v_Mock_72h_sub_1h = . then flag_1gy_72h_sub_1h=. ;
   else if p_1gy_v_Mock_72h_sub_1h < 0.05 then flag_1gy_72h_sub_1h=1 ;
   else  flag_1gy_72h_sub_1h=0 ;

run;

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_3h_sub_1h*flag_1gy_3h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*
        flag_01gy_    1gy_3h_
 Obs     3h_sub_1h     sub_1h    COUNT    PERCENT

  1          .           .        7326      .
  2          0           0       15863    72.9736
  3          0           1         886     4.0758
  4          1           0        1790     8.2344
  5          1           1        3199    14.7162


*/

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_24h_sub_1h*flag_1gy_24h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*

        flag_01gy_     flag_1gy_
 Obs    24h_sub_1h    24h_sub_1h    COUNT    PERCENT

  1          .             .         7326      .
  2          0             0        16259    74.7953
  3          0             1          871     4.0068
  4          1             0         1797     8.2666
  5          1             1         2811    12.9313


*/

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_72h_sub_1h*flag_1gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*

        flag_01gy_     flag_1gy_
 Obs    72h_sub_1h    72h_sub_1h    COUNT    PERCENT

  1          .             .         7326      .
  2          0             0        15442    71.0369
  3          0             1         1076     4.9499
  4          1             0         1658     7.6272
  5          1             1         3562    16.3861


*/
proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_3h_sub_1h*flag_01gy_24h_sub_1h*flag_01gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*

       flag_01gy_    flag_01gy_    flag_01gy_
Obs     3h_sub_1h    24h_sub_1h    72h_sub_1h    COUNT    PERCENT

 1          .             .             .         7326      .
 2          0             0             0        13267    61.0314
 3          0             0             1         1498     6.8912
 4          0             1             0         1150     5.2903
 5          0             1             1          834     3.8366
 6          1             0             0         1605     7.3834
 7          1             0             1          760     3.4962
 8          1             1             0          496     2.2817
 9          1             1             1         2128     9.7893



*/

proc freq data=flag_p05_diffs noprint;
   tables flag_1gy_3h_sub_1h*flag_1gy_24h_sub_1h*flag_1gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*          flag_
         1gy_3h_     flag_1gy_     flag_1gy_
  Obs     sub_1h    24h_sub_1h    72h_sub_1h    COUNT    PERCENT

   1        .            .             .         7326      .
   2        0            0             0        14297    65.7696
   3        0            0             1         1611     7.4110
   4        0            1             0          993     4.5680
   5        0            1             1          752     3.4594
   6        1            0             0         1429     6.5737
   7        1            0             1          719     3.3076
   8        1            1             0          381     1.7527
   9        1            1             1         1556     7.1580



*/

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_3h_sub_1h*flag_01gy_24h_sub_1h*flag_01gy_72h_sub_1h*flag_1gy_3h_sub_1h*flag_1gy_24h_sub_1h*flag_1gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;


/*
        flag_01gy_   flag_01gy_   flag_01gy_   1gy_3h_    flag_1gy_    flag_1gy_
  Obs    3h_sub_1h   24h_sub_1h   72h_sub_1h    sub_1h   24h_sub_1h   72h_sub_1h   COUNT   PERCENT

    1        .            .            .          .           .            .        7326     .
    2        0            0            0          0           0            0       11745   54.0298
    3        0            0            0          0           0            1         454    2.0885
    4        0            0            0          0           1            0         314    1.4445
    5        0            0            0          0           1            1          98    0.4508
    6        0            0            0          1           0            0         393    1.8079
    7        0            0            0          1           0            1         114    0.5244
    8        0            0            0          1           1            0          62    0.2852
    9        0            0            0          1           1            1          87    0.4002
   10        0            0            1          0           0            0         509    2.3415
   11        0            0            1          0           0            1         814    3.7446
   12        0            0            1          0           1            0           4    0.0184
   13        0            0            1          0           1            1          82    0.3772
   14        0            0            1          1           0            0          10    0.0460
   15        0            0            1          1           0            1          44    0.2024
   16        0            0            1          1           1            1          35    0.1610
   17        0            1            0          0           0            0         552    2.5393
   18        0            1            0          0           0            1          10    0.0460
   19        0            1            0          0           1            0         451    2.0747
   20        0            1            0          0           1            1          61    0.2806
   21        0            1            0          1           0            0          15    0.0690
   22        0            1            0          1           0            1           9    0.0414
   23        0            1            0          1           1            0          24    0.1104
   24        0            1            0          1           1            1          28    0.1288
  25        0            1            1          0           0            0        235    1.08106
  26        0            1            1          0           0            1        112    0.51523
  27        0            1            1          0           1            0         70    0.32202
  28        0            1            1          0           1            1        352    1.61928
  29        0            1            1          1           0            0         10    0.04600
  30        0            1            1          1           0            1          1    0.00460
  31        0            1            1          1           1            0          1    0.00460
  32        0            1            1          1           1            1         53    0.24381
  33        1            0            0          0           0            0        652    2.99936
  34        1            0            0          0           0            1         42    0.19321
  35        1            0            0          0           1            0         14    0.06440
  36        1            0            0          0           1            1         17    0.07820
  37        1            0            0          1           0            0        721    3.31677
  38        1            0            0          1           0            1         82    0.37722
  39        1            0            0          1           1            0         43    0.19781
  40        1            0            0          1           1            1         34    0.15641
  41        1            0            1          0           0            0        177    0.81424
  42        1            0            1          0           0            1         81    0.37262
  43        1            0            1          0           1            0          3    0.01380
  44        1            0            1          0           1            1         10    0.04600
  45        1            0            1          1           0            0        107    0.49223
  46        1            0            1          1           0            1        314    1.44448
  47        1            0            1          1           1            0          2    0.00920
  48        1            0            1          1           1            1         66    0.30362
  49        1            1            0          0           0            0         172   0.79124
  50        1            1            0          0           0            1           6   0.02760
  51        1            1            0          0           1            0          52   0.23921
  52        1            1            0          0           1            1           2   0.00920
  53        1            1            0          1           0            0          74   0.34042
  54        1            1            0          1           0            1           5   0.02300
  55        1            1            0          1           1            0         158   0.72684
  56        1            1            0          1           1            1          27   0.12421
  57        1            1            1          0           0            0         255   1.17306
  58        1            1            1          0           0            1          92   0.42322
  59        1            1            1          0           1            0          85   0.39102
  60        1            1            1          0           1            1         130   0.59803
  61        1            1            1          1           0            0          99   0.45542
  62        1            1            1          1           0            1         150   0.69004
  63        1            1            1          1           1            0          91   0.41862
  64        1            1            1          1           1            1        1226   5.63989

*/

/* New comparisons
Overlap between Treat vs. Mock for each timepoint*/


data flag_p05;
  set gene_summary2;
  if p_01gy_v_Mock_1h = . then flag_01vMock_1_p05=.;
  else if p_01gy_v_Mock_1h <0.05 then flag_01vMock_1_p05=1;
  else flag_01vMock_1_p05=0;

  if p_01gy_v_Mock_3h = . then flag_01vMock_3_p05=.;
  else if p_01gy_v_Mock_3h <0.05 then flag_01vMock_3_p05=1;
  else flag_01vMock_3_p05=0;

  if p_01gy_v_Mock_24h = . then flag_01vMock_24_p05=.;
  else if p_01gy_v_Mock_24h <0.05 then flag_01vMock_24_p05=1;
  else flag_01vMock_24_p05=0;

  if p_01gy_v_Mock_72h = . then flag_01vMock_72_p05=.;
  else if p_01gy_v_Mock_72h <0.05 then flag_01vMock_72_p05=1;
  else flag_01vMock_72_p05=0;


  if p_1gy_v_Mock_1h = . then flag_1vMock_1_p05=.;
  else if p_1gy_v_Mock_1h <0.05 then flag_1vMock_1_p05=1;
  else flag_1vMock_1_p05=0;

  if p_1gy_v_Mock_3h = . then flag_1vMock_3_p05=.;
  else if p_1gy_v_Mock_3h <0.05 then flag_1vMock_3_p05=1;
  else flag_1vMock_3_p05=0;

  if p_1gy_v_Mock_24h = . then flag_1vMock_24_p05=.;
  else if p_1gy_v_Mock_24h <0.05 then flag_1vMock_24_p05=1;
  else flag_1vMock_24_p05=0;

  if p_1gy_v_Mock_72h = . then flag_1vMock_72_p05=.;
  else if p_1gy_v_Mock_72h <0.05 then flag_1vMock_72_p05=1;
  else flag_1vMock_72_p05=0;

run;

proc freq data=flag_p05 noprint;
  tables flag_01vMock_1_p05*flag_01vMock_3_p05*flag_01vMock_24_p05*flag_01vMock_72_p05 / out=diff_count;
run;

proc print data=diff_count;
run;


/*
   flag_       flag_       flag_       flag_
 01vMock_    01vMock_    01vMock_    01vMock_
   1_p05       3_p05      24_p05      72_p05     COUNT

     .           .           .           .        7326
     0           0           0           0        9997
     0           0           0           1        1650
     0           0           1           0        1003
     0           0           1           1         558
     0           1           0           0        1217
     0           1           0           1         372
     0           1           1           0         157
     0           1           1           1          93
     1           0           0           0        3592
     1           0           0           1         899
     1           0           1           0         646
     1           0           1           1         370
     1           1           0           0         710
     1           1           0           1         261
     1           1           1           0         124
     1           1           1           1          89




*/

proc freq data=flag_p05 noprint;
  tables flag_1vMock_1_p05*flag_1vMock_3_p05*flag_1vMock_24_p05*flag_1vMock_72_p05 / out=diff_count;
run;

proc print data=diff_count;
run;

/*
   flag_      flag_      flag_1v     flag_1v
  1vMock_    1vMock_    Mock_24_    Mock_72_
   1_p05      3_p05        p05         p05      COUNT

     .          .           .           .        7326
     0          0           0           0       10845
     0          0           0           1        1512
     0          0           1           0        1103
     0          0           1           1         260
     0          1           0           0        1840
     0          1           0           1         466
     0          1           1           0         284
     0          1           1           1         157
     1          0           0           0        2826
     1          0           0           1         479
     1          0           1           0         479
     1          0           1           1         147
     1          1           0           0         883
     1          1           0           1         233
     1          1           1           0         164
     1          1           1           1          60


*/

proc freq data=flag_p05 noprint;
  tables sign_01gy_v_Mock_1h*sign_01gy_v_Mock_3h*sign_01gy_v_Mock_24h*sign_01gy_v_Mock_72h / out=diff_count;
run;

proc print data=diff_count;
run;

/*
        sign_01gy_    sign_01gy_    sign_01gy_    sign_01gy_
 Obs    v_Mock_1h     v_Mock_3h     v_Mock_24h    v_Mock_72h    COUNT    PERCENT

   1                                                             7326     .
   2        D             D             D             D            23    0.10581
   3        D             D             D             N            11    0.05060
   4        D             D             N             D            88    0.40482
   5        D             D             N             N           295    1.35707
   6        D             D             N             U            24    0.11041
   7        D             D             U             D             1    0.00460
   8        D             D             U             N            10    0.04600
   9        D             D             U             U             5    0.02300
  10        D             N             D             D           171    0.78664
  11        D             N             D             N           212    0.97525
  12        D             N             D             U             5    0.02300
  13        D             N             N             D           322    1.48128
  14        D             N             N             N          1901    8.74505
  15        D             N             N             U           158    0.72684
  16        D             N             U             D             6    0.02760
  17        D             N             U             N           136    0.62563
  18        D             N             U             U            46    0.21161
  19        D             U             D             D             4    0.01840
  20        D             U             D             N            30    0.13801
  21        D             U             D             U             2    0.00920
  22        D             U             N             D            10    0.04600
  23        D             U             N             N           140    0.64403
  24        D             U             N             U            26    0.11961
  25        D             U             U             N            10    0.04600
 26        D             U             U             U             3     0.0138
 27        N             D             D             D            40     0.1840
 28        N             D             D             N            30     0.1380
 29        N             D             D             U             1     0.0046
 30        N             D             N             D           155     0.7130
 31        N             D             N             N           567     2.6083
 32        N             D             N             U           137     0.6302
 33        N             D             U             D             2     0.0092
 34        N             D             U             N            52     0.2392
 35        N             D             U             U            32     0.1472
 36        N             N             D             D           371     1.7067
 37        N             N             D             N           421     1.9367
 38        N             N             D             U            10     0.0460
 39        N             N             N             D           987     4.5404
 40        N             N             N             N          9997    45.9886
 41        N             N             N             U           663     3.0500
 42        N             N             U             D            10     0.0460
 43        N             N             U             N           582     2.6773
 44        N             N             U             U           167     0.7682
 45        N             U             D             D             4     0.0184
 46        N             U             D             N            26     0.1196
 47        N             U             D             U             3     0.0138
 48        N             U             N             D            18     0.0828
 49        N             U             N             N           650     2.9902
 50        N             U             N             U            62     0.2852
   51        N             U             U             N            49    0.22541
   52        N             U             U             U            11    0.05060
   53        U             D             D             D             6    0.02760
   54        U             D             D             N            11    0.05060
   55        U             D             D             U             1    0.00460
   56        U             D             N             D            24    0.11041
   57        U             D             N             N           107    0.49223
   58        U             D             N             U            58    0.26681
   59        U             D             U             D             1    0.00460
   60        U             D             U             N            33    0.15181
   61        U             D             U             U            21    0.09661
   62        U             N             D             D            21    0.09661
   63        U             N             D             N            61    0.28061
   64        U             N             D             U             2    0.00920
   65        U             N             N             D           138    0.63483
   66        U             N             N             N          1691    7.77900
   67        U             N             N             U           281    1.29267
   68        U             N             U             D             4    0.01840
   69        U             N             U             N           237    1.09026
   70        U             N             U             U           115    0.52903
   71        U             U             D             N             4    0.01840
   72        U             U             N             D             5    0.02300
   73        U             U             N             N           168    0.77284
   74        U             U             N             U            26    0.11961
   75        U             U             U             N            15    0.06900
   76        U             U             U             U            22    0.10121
*/

proc freq data=flag_p05 noprint;
  tables sign_1gy_v_Mock_1h*sign_1gy_v_Mock_3h*sign_1gy_v_Mock_24h*sign_1gy_v_Mock_72h / out=diff_count;
run;

proc print data=diff_count;
run;

/*
         sign_      sign_
        1gy_v_     1gy_v_     sign_1gy_     sign_1gy_
 Obs    Mock_1h    Mock_3h    v_Mock_24h    v_Mock_72h    COUNT    PERCENT

   1                                                       7326     .
   2       D          D           D             D            24    0.11041
   3       D          D           D             N            27    0.12421
   4       D          D           N             D            75    0.34502
   5       D          D           N             N           248    1.14086
   6       D          D           N             U             2    0.00920
   7       D          D           U             D             3    0.01380
   8       D          D           U             N            40    0.18401
   9       D          N           D             D            48    0.22081
  10       D          N           D             N            87    0.40022
  11       D          N           D             U             2    0.00920
  12       D          N           N             D           156    0.71764
  13       D          N           N             N          1368    6.29313
  14       D          N           N             U            75    0.34502
  15       D          N           U             D             3    0.01380
  16       D          N           U             N           151    0.69464
  17       D          N           U             U            25    0.11501
  18       D          U           D             D             1    0.00460
  19       D          U           D             N             6    0.02760
  20       D          U           N             D             8    0.03680
  21       D          U           N             N           168    0.77284
  22       D          U           N             U            21    0.09661
  23       D          U           U             N            24    0.11041
  24       D          U           U             U             5    0.02300
  25       N          D           D             D           117     0.5382
  26       N          D           D             N           123     0.5658
  27       N          D           D             U             1     0.0046
  28       N          D           N             D           185     0.8510
  29       N          D           N             N           869     3.9976
  30       N          D           N             U            10     0.0460
  31       N          D           U             D             2     0.0092
  32       N          D           U             N            72     0.3312
  33       N          D           U             U             9     0.0414
  34       N          N           D             D           165     0.7590
  35       N          N           D             N           619     2.8475
  36       N          N           D             U             8     0.0368
  37       N          N           N             D           422     1.9413
  38       N          N           N             N         10845    49.8896
  39       N          N           N             U          1090     5.0143
  40       N          N           U             D             4     0.0184
  41       N          N           U             N           484     2.2265
  42       N          N           U             U            83     0.3818
  43       N          U           D             D             4     0.0184
  44       N          U           D             N            29     0.1334
  45       N          U           D             U             2     0.0092
  46       N          U           N             D            17     0.0782
  47       N          U           N             N           971     4.4668
  48       N          U           N             U           254     1.1685
 49       N          U           U             D             1    0.00460
 50       N          U           U             N            60    0.27601
 51       N          U           U             U            21    0.09661
 52       U          D           D             D             3    0.01380
 53       U          D           D             N             1    0.00460
 54       U          D           N             D            10    0.04600
 55       U          D           N             N           138    0.63483
 56       U          D           N             U             7    0.03220
 57       U          D           U             D             1    0.00460
 58       U          D           U             N            23    0.10581
 59       U          D           U             U             5    0.02300
 60       U          N           D             D             9    0.04140
 61       U          N           D             N            52    0.23921
 62       U          N           D             U             2    0.00920
 63       U          N           N             D            45    0.20701
 64       U          N           N             N          1458    6.70715
 65       U          N           N             U           203    0.93385
 66       U          N           U             D             6    0.02760
 67       U          N           U             N           189    0.86945
 68       U          N           U             U            52    0.23921
 69       U          U           D             D             3    0.01380
 70       U          U           D             N             4    0.01840
 71       U          U           D             U             1    0.00460
 72       U          U           N             D             4    0.01840
 73       U          U           N             N          329     1.51348
 74       U          U           N             U          106     0.48763
 75       U          U           U             N           39     0.17941
 76       U          U           U             U           14     0.06440

*/



proc freq data=flag_p05 ;
  tables sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h 
         sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h ;
run;

/*
  sign_
  01gy_v_                             Cumulative    Cumulative
  Mock_1h    Frequency     Percent     Frequency      Percent
  ------------------------------------------------------------
  D              3639       16.74          3639        16.74
  N             15047       69.22         18686        85.96
  U              3052       14.04         21738       100.00

                    Frequency Missing = 7326


  sign_
  01gy_v_                             Cumulative    Cumulative
  Mock_3h    Frequency     Percent     Frequency      Percent
  ------------------------------------------------------------
  D              1735        7.98          1735         7.98
  N             18715       86.09         20450        94.07
  U              1288        5.93         21738       100.00

  sign_01gy_                             Cumulative    Cumulative
  v_Mock_24h    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                 1470        6.76          1470         6.76
  N                18698       86.02         20168        92.78
  U                 1570        7.22         21738       100.00

                     Frequency Missing = 7326


  sign_01gy_                             Cumulative    Cumulative
  v_Mock_72h    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                 2411       11.09          2411        11.09
  N                17446       80.26         19857        91.35
  U                 1881        8.65         21738       100.00

   sign_
   1gy_v_
   Mock_                              Cumulative    Cumulative
   1h        Frequency     Percent     Frequency      Percent
   -----------------------------------------------------------
   D             2567       11.81          2567        11.81
   N            16467       75.75         19034        87.56
   U             2704       12.44         21738       100.00

                    Frequency Missing = 7326


   sign_
   1gy_v_
   Mock_                              Cumulative    Cumulative
   3h        Frequency     Percent     Frequency      Percent
   -----------------------------------------------------------
   D             1995        9.18          1995         9.18
   N            17651       81.20         19646        90.38
   U             2092        9.62         21738       100.00


 sign_
 1gy_v_
 Mock_                              Cumulative    Cumulative
 24h       Frequency     Percent     Frequency      Percent
 -----------------------------------------------------------
 D             1338        6.16          1338         6.16
 N            19084       87.79         20422        93.95
 U             1316        6.05         21738       100.00

                  Frequency Missing = 7326


 sign_
 1gy_v_
 Mock_                              Cumulative    Cumulative
 72h       Frequency     Percent     Frequency      Percent
 -----------------------------------------------------------
 D             1316        6.05          1316         6.05
 N            18424       84.75         19740        90.81
 U             1998        9.19         21738       100.00
*/

