/* Formatting fusion data for the 13 ethanol response genes for SVN and SS */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

/* Get fusion annotations */

data symbol2gene;
  set tair.gene2symbol_ens;
  keep symbol gene_id;
run;

data fus2gene;
  set tair.tair20_fusion_si_info;
  gene_id=compress(tranwrd(primary_fbgn,"gene:",""));
  keep fusion_id gene_id;
run;

proc sort data=symbol2gene nodup;
  by gene_id;
proc sort data=fus2gene nodup;
  by gene_id;
run;

data fus2keep;
  merge fus2gene (in=in1) symbol2gene (in=in2);
  by gene_id;
  if in1 and in2 then output;
  else if in1 then do;
      symbol=gene_id;
      output;
      end;
run;

proc sort data=fus2keep;
   by fusion_id gene_id;
proc freq data=fus2keep noprint;
  tables fusion_id/out=fus_count;
proc sort data=fus_count;
  by descending count;
run; *6 genes per fusion max;

data cat_fus2keep;
   array gene[20] $20.;
   array sym[20] $43.;
   retain gene1-gene20 sym1-sym20;
   set fus2keep;
   by fusion_id;
   if first.fusion_id then do;
      call missing(of gene1-gene20);
      call missing(of sym1-sym20);
      records = 0;
   end;
   records + 1;
   gene[records]=gene_id;
   sym[records]=symbol;
   if last.fusion_id then output;
run;

data cat_fus2keep2;
  set cat_fus2keep;
  length gene_id2 $450.;
  length symbol_cat $900.;
  gene_id2=catx("|", OF gene1-gene20);
  symbol_cat=catx("|", OF sym1-sym20);
  drop gene1-gene20 sym1-sym20 gene_id records symbol;
  rename gene_id2=gene_id;
run;

data fus_info;
  set tair.tair20_si_fusions_unique_flagged;
  keep chr fusion_start fusion_stop fusion_id strand flag_multigene flag_common flag_constitutive flag_alternative;
run;

proc sort data=cat_fus2keep2;
   by fusion_id;
proc sort data=fus_info;
   by fusion_id;
run;

data fus2keep_w_info;
  merge cat_fus2keep2 (in=in1) fus_info (in=in2);
  by fusion_id;
  if in1 and in2;
run;

/* Get P-values for the following:
   treatment   time trt*time

 */

data main_p;
  set rs.arab_anova_main_trt_tm_non;
  length effect_p $30.;
   format ProbF best32. ;
  if effect="treatment" then effect_p="P_treatment";
  if effect="time" then effect_p="P_time";
  if effect="time*treatment" then effect_p="P_trt_by_time";
  keep fusion_id effect_p probf;
run;

proc sort data=main_p;
   by fusion_id effect_p;
proc transpose data=main_p out=main_p_sbys;
  by fusion_id;
  id effect_p;
  var probf;
run;

data main_p_sbys2;
   set main_p_sbys;
   drop _NAME_ _LABEL_;
run;

data contrast_p;
  set rs.arab_fus_cntrs_constr_non;
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

  keep fusion_id p_label probf;
run;

proc sort data=contrast_p;
   by fusion_id p_label;
proc transpose data=contrast_p out=contrast_p_sbys;
   by fusion_id;
   id p_label;
   var probf;
run;


data contrast_p_sbys2;
   set contrast_p_sbys;
   drop _NAME_ _LABEL_;
run;


/* Get LSmeans */

data lsmeans;
  set rs.arab_fus_cntrs_lsmeans_non;
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
  keep fusion_id group estimate;
run;

proc sort data=lsmeans;
   by fusion_id group;
proc transpose data=lsmeans out=lsmeans_sbys;
   by fusion_id;
   id group;
   var estimate;
run;


data lsmeans_sbys2;
   set lsmeans_sbys;
   drop _NAME_ _LABEL_;
run;


/* Get estimates */

data estimates;
  set rs.arab_fus_cntrs_estim_non;
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

  keep fusion_id est_label estimate;
run;

proc sort data=estimates;
   by fusion_id est_label;
proc transpose data=estimates out=estimates_sbys;
   by fusion_id;
   id est_label;
   var estimate;
run;


data estimates_sbys2;
   set estimates_sbys;
   drop _NAME_ _LABEL_;
run;


/* detection flags */

data flag_dtct;
   set rs.arab_flag_fusion_on_gt0;
   keep fusion_id    flag_fusion_on_01gy_apn0    flag_fusion_on_1gy_apn0
   flag_fusion_on_Mock_apn0    flag_fusion_on_1hr_apn0   flag_fusion_on_3hr_apn0
   flag_fusion_on_34hr_apn0   flag_fusion_on_72hr_apn0   flag_fusion_on_01gy_1hr_apn0
   flag_fusion_on_01gy_3hr_apn0   flag_fusion_on_01gy_24hr_apn0
   flag_fusion_on_01gy_72hr_apn0   flag_fusion_on_1gy_1hr_apn0
   flag_fusion_on_1gy_3hr_apn0   flag_fusion_on_1gy_24hr_apn0
   flag_fusion_on_1gy_72hr_apn0   flag_fusion_on_Mock_1hr_apn0
   flag_fusion_on_Mock_3hr_apn0   flag_fusion_on_Mock_24hr_apn0
   flag_fusion_on_Mock_72hr_apn0;

   rename flag_fusion_on_01gy_apn0=flag_01gy_on
   flag_fusion_on_1gy_apn0=flag_1gy_on
   flag_fusion_on_Mock_apn0=flag_Mock_on
   flag_fusion_on_1hr_apn0=flag_1hr_on
   flag_fusion_on_3hr_apn0=flag_3hr_on
   flag_fusion_on_34hr_apn0=flag_24hr_on
   flag_fusion_on_72hr_apn0=flag_72gy_on
   flag_fusion_on_01gy_1hr_apn0=flag_01gy_1hr_on
   flag_fusion_on_01gy_3hr_apn0=flag_01gy_3hr_on
   flag_fusion_on_01gy_24hr_apn0=flag_01gy_24hr_on
   flag_fusion_on_01gy_72hr_apn0=flag_01gy_72hr_on
   flag_fusion_on_1gy_1hr_apn0=flag_1gy_1hr_on
   flag_fusion_on_1gy_3hr_apn0=flag_1gy_3hr_on
   flag_fusion_on_1gy_24hr_apn0=flag_1gy_24hr_on
   flag_fusion_on_1gy_72hr_apn0=flag_1gy_72hr_on
   flag_fusion_on_Mock_1hr_apn0=flag_Mock_1hr_on
   flag_fusion_on_Mock_3hr_apn0=flag_Mock_3hr_on
   flag_fusion_on_Mock_24hr_apn0=flag_Mock_24hr_on
   flag_fusion_on_Mock_72hr_apn0=flag_Mock_72hr_on;

run;


/* Get Up/Down signs */

data sign;
  set rs.arab_sign_by_contrast;
run;


proc sort data=fus2keep_w_info;
   by fusion_id;
proc sort data=main_p_sbys2;
   by fusion_id;
proc sort data=contrast_p_sbys2;
   by fusion_id;
proc sort data=lsmeans_sbys2;
   by fusion_id;
proc sort data=estimates_sbys2;
   by fusion_id;
proc sort data=sign;
   by fusion_id;
proc sort data=flag_dtct;
   by fusion_id;
run;


data fus_summary;
   merge fus2keep_w_info (in=in1) flag_dtct main_p_sbys2 contrast_p_sbys2 lsmeans_sbys2 estimates_sbys2 sign;
   by fusion_id;
   if in1;
run;

/* Reorder variables */

proc contents data=fus_summary;
run;


data fus_summary2;
  retain fusion_id chr fusion_start fusion_stop strand flag_multigene gene_id symbol_cat flag_Alternative
         flag_Common flag_Constitutive flag_Mock_on flag_01gy_on flag_1gy_on flag_1hr_on flag_3hr_on
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
         set fus_summary;
rename  fusion_id=exonic_region_id fusion_start=exonic_region_start fusion_stop=exonic_region_stop;
run;

/* Export */

proc export data=fus_summary2 outfile="!PATCON/arabidopsis/analysis_output/arab_rs_fusions_tair10_all_genes.csv"
dbms=csv replace;
run;


/* Need to count number of sig fusions by effect/contrast, and number on per trt, time and trt*time */

proc freq data=fus_summary2 noprint;
   tables flag_Mock_on*flag_01gy_on*flag_1gy_on / out=count_on_by_trt;
run;

proc freq data=fus_summary2 noprint;
   tables flag_1hr_on*flag_3hr_on*flag_24hr_on*flag_72gy_on / out=count_on_by_time;
run;

proc freq data=fus_summary2 noprint;
   tables flag_Mock_1hr_on*flag_Mock_3hr_on*flag_Mock_24hr_on*flag_Mock_72hr_on*
          flag_01gy_1hr_on*flag_01gy_3hr_on*flag_01gy_24hr_on*flag_01gy_72hr_on*
          flag_1gy_1hr_on*flag_1gy_3hr_on*flag_1gy_24hr_on*flag_1gy_72hr_on / out=count_on_by_trt_time;
run;

proc export data=count_on_by_trt outfile="!PATCON/arabidopsis/analysis_output/fusions_on_by_trt_apn0.csv"
    dbms=csv replace;
run;

proc export data=count_on_by_time outfile="!PATCON/arabidopsis/analysis_output/fusions_on_by_time_apn0.csv"
    dbms=csv replace;
run;

proc export data=count_on_by_trt_time outfile="!PATCON/arabidopsis/analysis_output/fusions_on_by_trt_time_apn0.csv"
    dbms=csv replace;
run;


/* Now flag and count P<0.05 for each main effect and contrast */

%macro flagp05(pvalue,flagname);
data flag_p05;
  set fus_summary2;
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
             0       31196       26.62         31196        26.62
             1       86000       73.38        117196       100.00

                      Frequency Missing = 21371

                                         Cumulative    Cumulative
flag_trt_p05    Frequency     Percent     Frequency      Percent
-----------------------------------------------------------------
           0       95812       81.75         95812        81.75
           1       21384       18.25        117196       100.00

                    Frequency Missing = 21371

        flag_trt_                             Cumulative    Cumulative
         time_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       97729       83.39         97729        83.39
                1       19467       16.61        117196       100.00

                       Frequency Missing = 21371

   flag_01vMock_                             Cumulative    Cumulative
           1_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       94945       81.01         94945        81.01
               1       22251       18.99        117196       100.00

   flag_01vMock_                             Cumulative    Cumulative
           3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      106556       90.92        106556        90.92
               1       10640        9.08        117196       100.00


   flag_01vMock_                             Cumulative    Cumulative
          24_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      105853       90.32        105853        90.32
               1       11343        9.68        117196       100.00


   flag_01vMock_                             Cumulative    Cumulative
          72_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      103391       88.22        103391        88.22
               1       13805       11.78        117196       100.00

     flag_1vMock_                             Cumulative    Cumulative
            1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      100170       85.47        100170        85.47
                1       17026       14.53        117196       100.00


    flag_1vMock_                             Cumulative    Cumulative
           3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      101698       86.78        101698        86.78
               1       15498       13.22        117196       100.00


     flag_1vMock_                             Cumulative    Cumulative
           24_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      108072       92.21        108072        92.21
                1        9124        7.79        117196       100.00

     flag_1vMock_                             Cumulative    Cumulative
           72_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      107642       91.85        107642        91.85
                1        9554        8.15        117196       100.00

     flag_1v01gy_                             Cumulative    Cumulative
            1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      110255       94.08        110255        94.08
                1        6941        5.92        117196       100.00

   flag_1v01gy_                             Cumulative    Cumulative
          3_p05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       99678       85.05         99678        85.05
              1       17518       14.95        117196       100.00


                   flag_1v01gy_                             Cumulative    Cumulative
                         24_p05    Frequency     Percent     Frequency      Percent
               ---------------------------------------------------------------------
                              0      106219       90.63        106219        90.63
                              1       10977        9.37        117196       100.00

                                     Frequency Missing = 21371

     flag_1v01gy_                             Cumulative    Cumulative
           72_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      110967       94.68        110967        94.68
                1        6229        5.32        117196       100.00


      flag_mock_                             Cumulative    Cumulative
         3v1_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       81387       69.45         81387        69.45
               1       35809       30.55        117196       100.00

       flag_mock_                             Cumulative    Cumulative
         24v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       81296       69.37         81296        69.37
                1       35900       30.63        117196       100.00

        flag_mock_                             Cumulative    Cumulative
          72v1_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       89228       76.14         89228        76.14
                 1       27968       23.86        117196       100.00

       flag_mock_                             Cumulative    Cumulative
         24v3_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       70174       59.88         70174        59.88
                1       47022       40.12        117196       100.00


       flag_mock_                             Cumulative    Cumulative
         72v3_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       94955       81.02         94955        81.02
                1       22241       18.98        117196       100.00

       flag_mock_                             Cumulative    Cumulative
        72v24_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       80394       68.60         80394        68.60
                1       36802       31.40        117196       100.00

      flag_01gy_                             Cumulative    Cumulative
         3v1_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       70007       59.73         70007        59.73
               1       47189       40.27        117196       100.00

       flag_01gy_                             Cumulative    Cumulative
         24v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       67617       57.70         67617        57.70
                1       49579       42.30        117196       100.00.

       flag_01gy_                             Cumulative    Cumulative
         72v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       75940       64.80         75940        64.80
                1       41256       35.20        117196       100.00

      flag_01gy_                             Cumulative    Cumulative
        24v3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       65053       55.51         65053        55.51
               1       52143       44.49        117196       100.00

       flag_01gy_                             Cumulative    Cumulative
         72v3_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       83271       71.05         83271        71.05
                1       33925       28.95        117196       100.00

      flag_01gy_                             Cumulative    Cumulative
       72v24_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       76435       65.22         76435        65.22
               1       40761       34.78        117196       100.00

                                              Cumulative    Cumulative
 flag_1gy_3v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       71766       61.24         71766        61.24
                1       45430       38.76        117196       100.00

        flag_1gy_                             Cumulative    Cumulative
         24v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       72439       61.81         72439        61.81
                1       44757       38.19        117196       100.00

        flag_1gy_                             Cumulative    Cumulative
         72v1_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       80961       69.08         80961        69.08
                1       36235       30.92        117196       100.00

       flag_1gy_                             Cumulative    Cumulative
        24v3_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       70788       60.40         70788        60.40
               1       46408       39.60        117196       100.00

        flag_1gy_                             Cumulative    Cumulative
         72v3_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       89413       76.29         89413        76.29
                1       27783       23.71        117196       100.00

      flag_1gy_                             Cumulative    Cumulative
      72v24_p05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       77960       66.52         77960        66.52
              1       39236       33.48        117196       100.00

       flag_01gy_v_
       Mock_3h_sub_                             Cumulative    Cumulative
             1h_p05    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0      100318       85.60        100318        85.60
                  1       16878       14.40        117196       100.00

      flag_01gy_v_
     Mock_24h_sub_                             Cumulative    Cumulative
            1h_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0      101831       86.89        101831        86.89
                 1       15365       13.11        117196       100.00

      flag_01gy_v_
     Mock_72h_sub_                             Cumulative    Cumulative
            1h_p05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0      100504       85.76        100504        85.76
                 1       16692       14.24        117196       100.00


flag_1gy_v_Mock_                             Cumulative    Cumulative
   3h_sub_1h_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      103683       88.47        103683        88.47
               1       13513       11.53        117196       100.00

flag_1gy_v_Mock_                             Cumulative    Cumulative
  24h_sub_1h_p05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0      104730       89.36        104730        89.36
               1       12466       10.64        117196       100.00

 flag_1gy_v_Mock_                             Cumulative    Cumulative
   72h_sub_1h_p05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0      102851       87.76        102851        87.76
                1       14345       12.24        117196       100.00
*/

/* Check: X-1 Y vs Mock: 0.1x1 for each time point, then within treatment */

data flag_p05_diffs;
   set fus_summary2;
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
               flag_
flag_01gy_    1gy_3h_
 3h_sub_1h     sub_1h    COUNT    PERCENT

     .           .       21371      .
     0           0       96534    82.3697
     0           1        3784     3.2288
     1           0        7149     6.1000
     1           1        9729     8.3015

*/

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_24h_sub_1h*flag_1gy_24h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*

flag_01gy_     flag_1gy_
24h_sub_1h    24h_sub_1h    COUNT    PERCENT

     .             .        21371      .
     0             0        97889    83.5259
     0             1         3942     3.3636
     1             0         6841     5.8372
     1             1         8524     7.2733

*/

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_72h_sub_1h*flag_1gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*
 flag_01gy_     flag_1gy_
 72h_sub_1h    72h_sub_1h    COUNT    PERCENT

      .             .        21371      .
      0             0        96056    81.9618
      0             1         4448     3.7954
      1             0         6795     5.7980
      1             1         9897     8.4448


*/
proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_3h_sub_1h*flag_01gy_24h_sub_1h*flag_01gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*
 flag_01gy_    flag_01gy_    flag_01gy_
  3h_sub_1h    24h_sub_1h    72h_sub_1h    COUNT    PERCENT

      .             .             .        21371      .
      0             0             0        86225    73.5733
      0             0             1         5899     5.0334
      0             1             0         5553     4.7382
      0             1             1         2641     2.2535
      1             0             0         6894     5.8825
      1             0             1         2813     2.4003
      1             1             0         1832     1.5632
      1             1             1         5339     4.5556



*/

proc freq data=flag_p05_diffs noprint;
   tables flag_1gy_3h_sub_1h*flag_1gy_24h_sub_1h*flag_1gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;

/*
  flag_
 1gy_3h_     flag_1gy_     flag_1gy_
  sub_1h    24h_sub_1h    72h_sub_1h    COUNT    PERCENT

    .            .             .        21371      .
    0            0             0        90450    77.1784
    0            0             1         5946     5.0736
    0            1             0         5018     4.2817
    0            1             1         2269     1.9361
    1            0             0         6048     5.1606
    1            0             1         2286     1.9506
    1            1             0         1335     1.1391
    1            1             1         3844     3.2800

*/

proc freq data=flag_p05_diffs noprint;
   tables flag_01gy_3h_sub_1h*flag_01gy_24h_sub_1h*flag_01gy_72h_sub_1h*flag_1gy_3h_sub_1h*flag_1gy_24h_sub_1h*flag_1gy_72h_sub_1h / out=p05_count;
proc print data=p05_count;
run;


/*
flag_01gy_   flag_01gy_   flag_01gy_   1gy_3h_    flag_1gy_    flag_1gy_
 3h_sub_1h   24h_sub_1h   72h_sub_1h    sub_1h   24h_sub_1h   72h_sub_1h   COUNT   PERCENT

     .            .            .          .           .            .       21371     .
     0            0            0          0           0            0       78716   67.1661
     0            0            0          0           0            1        2179    1.8593
     0            0            0          0           1            0        1830    1.5615
     0            0            0          0           1            1         474    0.4045
     0            0            0          1           0            0        1875    1.5999
     0            0            0          1           0            1         468    0.3993
     0            0            0          1           1            0         287    0.2449
     0            0            0          1           1            1         396    0.3379
     0            0            1          0           0            0        2504    2.1366
     0            0            1          0           0            1        2761    2.3559
     0            0            1          0           1            0          26    0.0222
     0            0            1          0           1            1         280    0.2389
     0            0            1          1           0            0          47    0.0401
     0            0            1          1           0            1         172    0.1468
     0            0            1          1           1            0           5    0.0043
     0            0            1          1           1            1         104    0.0887
     0            1            0          0           0            0        2743    2.3405
     0            1            0          0           0            1          34    0.0290
     0            1            0          0           1            0        2320    1.9796
     0            1            0          0           1            1         180    0.1536
     0            1            0          1           0            0          58    0.0495
     0            1            0          1           0            1          15    0.0128
     0            1            0          1           1            0         106    0.0904
     0            1            0          1           1            1          97   0.08277
     0            1            1          0           0            0         887   0.75685
     0            1            1          0           0            1         315   0.26878
     0            1            1          0           1            0         315   0.26878
     0            1            1          0           1            1         970   0.82767
     0            1            1          1           0            0          17   0.01451
     0            1            1          1           0            1           4   0.00341
     0            1            1          1           1            1         133   0.11349
     1            0            0          0           0            0        3164   2.69975
     1            0            0          0           0            1          88   0.07509
     1            0            0          0           1            0          67   0.05717
     1            0            0          0           1            1          42   0.03584
     1            0            0          1           0            0        3057   2.60845
     1            0            0          1           0            1         244   0.20820
     1            0            0          1           1            0         123   0.10495
     1            0            0          1           1            1         109   0.09301
     1            0            1          0           0            0         888   0.75771
     1            0            1          0           0            1         296   0.25257
     1            0            1          0           1            0          14   0.01195
     1            0            1          0           1            1          19   0.01621
     1            0            1          1           0            0         424   0.36179
     1            0            1          1           0            1        1006   0.85839
     1            0            1          1           1            0           1   0.00085
     1            0            1          1           1            1         165   0.14079
     1            1            0          0           0            0         694   0.59217
     1            1            0          0           0            1          19   0.01621
     1            1            0          0           1            0         193   0.16468
     1            1            0          0           1            1           3   0.00256
     1            1            0          1           0            0         290   0.24745
     1            1            0          1           0            1           4   0.00341
     1            1            0          1           1            0         533   0.45479
     1            1            0          1           1            1          96   0.08191
     1            1            1          0           0            0         854   0.72869
     1            1            1          0           0            1         254   0.21673
     1            1            1          0           1            0         253   0.21588
     1            1            1          0           1            1         301   0.25683
     1            1            1          1           0            0         280   0.23892
     1            1            1          1           0            1         373   0.31827
     1            1            1          1           1            0         280   0.23892
     1            1            1          1           1            1        2744   2.34138
*/
