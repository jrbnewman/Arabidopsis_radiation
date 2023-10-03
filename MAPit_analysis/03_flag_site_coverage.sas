libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Flag sites (by site type and condition) that have at least 10X coverage and evidence of methylation across both reps
   We will later only analyze those that have coverage in 0Gy and at least one other dose (or common ones, not sure) 
*/


* Native methylation;
data methyl_data;
    set wgbslocA.methylation_data_cg_chg_chh;
run;

proc sort data=methyl_data;
   by chr start_pos stop_pos site_type treatment units;
proc means data=methyl_data noprint;
   by chr start_pos stop_pos site_type treatment units;
   var total_C perc_methyl;
   output out=methyl_data_sum_reps sum(total_C)=total_C_summed mean(perc_methyl)=perc_methyl_avg;
run;

data flag_coverage;
   set methyl_data_sum_reps;
   if total_C_summed >= 10 then flag_10x=1; else flag_10x=0;
   if perc_methyl_avg > 0 then flag_methylated=1; else flag_methylated=0;
   if _FREQ_ = 2 then flag_both_reps=1; else flag_both_reps=0;
run;


proc sort data=flag_coverage;
   by chr start_pos stop_pos site_type treatment units;
run;

/* 10X coverage -- these will be the baseline set of sites to analyze */

proc transpose data=flag_coverage out=flag_coverage_sbys;
   by chr start_pos stop_pos site_type ;
   var flag_10x;
   id treatment units;
run;

data flag_coverage_sbys2;
  set flag_coverage_sbys;
  if _01Gy0U = . then _01Gy0U=0;
  if _1Gy0U = . then _1Gy0U=0;
  if _0Gy0U = . then _0Gy0U=0;
  rename _01Gy0U=flag_01Gy_coverage_ge_10x
         _1Gy0U=flag_1Gy_coverage_ge_10x
         _0Gy0U=flag_0Gy_coverage_ge_10x;
  drop _NAME_;
run;

/* methylation -- only do for sites with at least 10X coverage*/

proc transpose data=flag_coverage out=flag_meth_sbys;
   by chr start_pos stop_pos site_type ;
   var flag_methylated;
   id treatment units;
run;

data flag_meth_sbys2;
  set flag_meth_sbys;
  rename _01Gy0U=flag_01Gy_methyl_gt0
         _1Gy0U=flag_1Gy_methyl_gt0
         _0Gy0U=flag_0Gy_methyl_gt0;
  drop _NAME_;
run;



/* flag if from both reps */

proc transpose data=flag_coverage out=flag_reps_sbys;
   by chr start_pos stop_pos site_type ;
   var flag_both_reps;
   id treatment units;
run;

data flag_reps_sbys2;
  set flag_reps_sbys;
  rename _01Gy0U=flag_01Gy_both_reps
         _1Gy0U=flag_1Gy_both_reps
         _0Gy0U=flag_0Gy_both_reps;
  drop _NAME_;
run;

/* Make permanent */

data wgbslocA.flag_coverage_10x_cg_chg_chh;
  set flag_coverage_sbys2;
run;

data wgbslocA.flag_methylated_cg_chg_chh;
  set flag_meth_sbys2;
run;

data wgbslocA.flag_bothreps_cg_chg_chh;
  set flag_reps_sbys2;
run;



/* GC methylation */
data methyl_data;
    set wgbslocA.methylation_data_gc;
run;

proc sort data=methyl_data;
   by chr start_pos stop_pos site_type treatment units;
proc means data=methyl_data noprint;
   by chr start_pos stop_pos site_type treatment units;
   var total_C perc_methyl_norm;
   output out=methyl_data_sum_reps sum(total_C)=total_C_summed mean(perc_methyl)=perc_methyl_avg;
run;

data flag_coverage;
   set methyl_data_sum_reps;
   if total_C_summed >= 10 then flag_10x=1; else flag_10x=0;
   if perc_methyl_avg > 0 then flag_methylated=1; else flag_methylated=0;
   if _FREQ_ = 2 then flag_both_reps=1; else flag_both_reps=0;
run;


proc sort data=flag_coverage;
   by chr start_pos stop_pos site_type treatment units;
run;

/* 10X coverage -- these will be the baseline set of sites to analyze */

proc transpose data=flag_coverage out=flag_coverage_sbys;
   by chr start_pos stop_pos site_type ;
   var flag_10x;
   id treatment units;
run;

data flag_coverage_sbys2;
  set flag_coverage_sbys;

  if _01Gy0U = . then _01Gy0U=0;
  if _1Gy0U = . then _1Gy0U=0;
  if _0Gy0U = . then _0Gy0U=0;

  if _01Gy100U = . then _01Gy100U=0;
  if _1Gy100U = . then _1Gy100U=0;
  if _0Gy100U = . then _0Gy100U=0;

  rename _01Gy0U=flag_01Gy_0U_coverage_ge_10x
         _1Gy0U=flag_1Gy_0U_coverage_ge_10x
         _0Gy0U=flag_0Gy_0U_coverage_ge_10x
         _01Gy100U=flag_01Gy_100U_coverage_ge_10x
         _1Gy100U=flag_1Gy_100U_coverage_ge_10x
         _0Gy100U=flag_0Gy_100U_coverage_ge_10x;
  drop _NAME_;
run;

/* methylation -- only do for sites with at least 10X coverage*/

proc transpose data=flag_coverage out=flag_meth_sbys;
   by chr start_pos stop_pos site_type ;
   var flag_methylated;
   id treatment units;
run;

data flag_meth_sbys2;
  set flag_meth_sbys;
  rename _01Gy0U=flag_01Gy_0U_methyl_gt0
         _1Gy0U=flag_1Gy_0U_methyl_gt0
         _0Gy0U=flag_0Gy_0U_methyl_gt0
         _01Gy100U=flag_01Gy_100U_methyl_gt0
         _1Gy100U=flag_1Gy_100U_methyl_gt0
         _0Gy100U=flag_0Gy_100U_methyl_gt0;
  drop _NAME_;
run;



/* flag if from both reps */

proc transpose data=flag_coverage out=flag_reps_sbys;
   by chr start_pos stop_pos site_type ;
   var flag_both_reps;
   id treatment units;
run;

data flag_reps_sbys2;
  set flag_reps_sbys;
  rename _01Gy0U=flag_01Gy_0U_both_reps
         _1Gy0U=flag_1Gy_0U_both_reps
         _0Gy0U=flag_0Gy_0U_both_reps
         _01Gy100U=flag_01Gy_100U_both_reps
         _1Gy100U=flag_1Gy_100U_both_reps
         _0Gy100U=flag_0Gy_100U_both_reps;
  drop _NAME_;
run;

/* Make permanent */

data wgbslocA.flag_coverage_10x_gc;
  set flag_coverage_sbys2;
run;

data wgbslocA.flag_methylated_gc;
  set flag_meth_sbys2;
run;

data wgbslocA.flag_coverage_bothreps_gc;
  set flag_reps_sbys2;
run;


