libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

*libname wgbslocA '/blue/concannon/share/jnewman/mingqi_arab/sas_data';
ods listing;
ods html close;
options nosyntaxcheck;



proc datasets lib=work kill noprint;
run;
quit;


data all_all_dmc_tests_01gy;
  set wgbsA.all_dmc_tests_01gy_compare;
   
if flag_FET_FDR01=. then flag_FET_FDR01=0;
if flag_FET_FDR05=. then flag_FET_FDR05=0;
if flag_binomial_FDR01=. then flag_binomial_FDR01=0;
if flag_binomial_FDR05=. then flag_binomial_FDR05=0;
if flag_dss_fdr05=. then flag_dss_fdr05=0;
if flag_logit_FDR01=. then flag_logit_FDR01=0;
if flag_logit_FDR05=. then flag_logit_FDR05=0;
if flag_methylkit_q05=. then flag_methylkit_q05=0;
if flag_methylsig_bb_fdr05=. then flag_methylsig_bb_fdr05=0;
if flag_methylsig_bin_fdr05=. then flag_methylsig_bin_fdr05=0;
if flag_metilene_q05=. then flag_metilene_q05=0;


run;

data all_all_dac_tests_01gy;
  set wgbsA.all_dac_tests_01gy_compare;
   if CMH_FDR_P_01Gy_0gy_v2 ne . and CMH_FDR_P_01Gy_0gy_v2 < 0.05 then flag_CMH2_FDR_P05=1;
   else  flag_CMH2_FDR_P05=0;
   if CMH_FDR_P_01Gy_0gy_v2 ne . and CMH_FDR_P_01Gy_0gy_v2 < 0.01 then flag_CMH2_FDR_P01=1;
   else  flag_CMH2_FDR_P01=0;

   if BD_FDR_P_01Gy_0gy ne . and BD_FDR_P_01Gy_0gy < 0.05 then flag_BD_FDR_P05=1;
   else  flag_BD_FDR_P05=0;
   if BD_FDR_P_01Gy_0gy ne . and BD_FDR_P_01Gy_0gy < 0.01 then flag_BD_FDR_P01=1;
   else  flag_BD_FDR_P01=0;

   if BD_FDR_P_01Gy_0gy_v2 ne . and BD_FDR_P_01Gy_0gy_v2 < 0.05 then flag_BD2_FDR_P05=1;
   else  flag_BD_FDR_P05=0;
   if BD_FDR_P_01Gy_0gy_v2 ne . and BD_FDR_P_01Gy_0gy_v2 < 0.01 then flag_BD2_FDR_P01=1;
   else  flag_BD_FDR_P01=0;
   
if flag_FET_FDR01=. then flag_FET_FDR01=0;
if flag_FET_FDR05=. then flag_FET_FDR05=0;
if flag_binomial_FDR01=. then flag_binomial_FDR01=0;
if flag_binomial_FDR05=. then flag_binomial_FDR05=0;
if flag_logit_FDR01=. then flag_logit_FDR01=0;
if flag_logit_FDR05=. then flag_logit_FDR05=0;
if flag_methylkit_q05=. then flag_methylkit_q05=0;
if flag_metilene_q05=. then flag_metilene_q05=0;
if flag_BD2_FDR_P01=. then flag_BD2_FDR_P01=0;
if flag_BD2_FDR_P05=. then flag_BD2_FDR_P05=0;
if flag_BD_FDR_P01=. then flag_BD_FDR_P01=0;
if flag_BD_FDR_P05=. then flag_BD_FDR_P05=0;
if flag_CMH2_FDR_P01=. then flag_CMH2_FDR_P01=0;
if flag_CMH2_FDR_P05=. then flag_CMH2_FDR_P05=0;
if flag_CMH_FDR_P01=. then flag_CMH_FDR_P01=0;
if flag_CMH_FDR_P05 =. then flag_CMH_FDR_P05 =0;
if flag_FET_either_fdr01=. then flag_FET_either_fdr01=0;
if flag_FET_either_fdr05=. then flag_FET_either_fdr05=0;
if flag_dss_int_fdr05=. then flag_dss_int_fdr05=0;
if flag_dss_trt_fdr05=. then flag_dss_trt_fdr05=0;
if flag_methylkit_q05=. then flag_methylkit_q05=0;
if flag_metilene_noneg_q05=. then flag_metilene_noneg_q05=0;
if flag_metilene_q05=. then flag_metilene_q05=0;

run;




/* DMC counts */


proc contents data=all_all_dmc_tests_01gy;
run;

proc freq data=all_all_dmc_tests_01gy;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc
flag_methdiff_gt_50perc flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
tables site_type * flag_FET_FDR05 * (flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc
flag_methdiff_gt_50perc flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc freq data=all_all_dmc_tests_01gy;
 where flag_methdiff_gt_10perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc freq data=all_all_dmc_tests_01gy;
 where flag_methdiff_gt_20perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);

run;

proc freq data=all_all_dmc_tests_01gy;
 where flag_methdiff_gt_30perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc freq data=all_all_dmc_tests_01gy;
 where flag_methdiff_gt_50perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc contents data=all_all_dac_tests_01gy;
run;



proc freq data=all_all_dac_tests_01gy;
 tables 
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc
flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05
flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05 
flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05 flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05;
run;

proc freq data=all_all_dac_tests_01gy;
 where flag_FET_either_fdr05 = 1;
 tables 
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc
flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05
flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05 
flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05  flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05;
run;


proc freq data=all_all_dac_tests_01gy;
 tables (flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc) * (flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05 flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05  flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05 flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05);
run;


proc freq data=all_all_dac_tests_01gy;
	 where flag_FET_either_fdr05=1;
	 tables (flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc) * (flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05 flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05  flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
	 flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05 flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05);
 run;



data export_dac_tests;
  set all_all_dac_tests_01gy;
if flag_BD2_FDR_P01=. then delete;
if flag_BD2_FDR_P05=. then delete;
if flag_BD_FDR_P01=. then delete;
if flag_BD_FDR_P05=. then delete;
if flag_CMH2_FDR_P01=. then delete;
if flag_CMH2_FDR_P05=. then delete;
if flag_CMH_FDR_P01=. then delete;
if flag_CMH_FDR_P05 =. then delete;
if flag_FET_either_fdr01=. then delete;
if flag_FET_either_fdr05=. then delete;
if flag_dss_int_fdr05=. then delete;
if flag_dss_trt_fdr05=. then delete;
if flag_methylkit_q05=. then delete;
if flag_metilene_noneg_q05=. then delete;
if flag_metilene_q05=. then delete;
run;

*proc export data=export_dac_tests
*outfile="/TB14/TB14/sandbox/dtra_sandbox/DAC_tests_for_comparsion.csv"
*dbms=csv replace;
*run;








data all_all_dmc_tests_1gy;
  set wgbsA.all_dmc_tests_1gy_compare;
   
if flag_FET_FDR01=. then flag_FET_FDR01=0;
if flag_FET_FDR05=. then flag_FET_FDR05=0;
if flag_binomial_FDR01=. then flag_binomial_FDR01=0;
if flag_binomial_FDR05=. then flag_binomial_FDR05=0;
if flag_dss_fdr05=. then flag_dss_fdr05=0;
if flag_logit_FDR01=. then flag_logit_FDR01=0;
if flag_logit_FDR05=. then flag_logit_FDR05=0;
if flag_methylkit_q05=. then flag_methylkit_q05=0;
if flag_methylsig_bb_fdr05=. then flag_methylsig_bb_fdr05=0;
if flag_methylsig_bin_fdr05=. then flag_methylsig_bin_fdr05=0;
if flag_metilene_q05=. then flag_metilene_q05=0;

run;

data all_all_dac_tests_1gy;
  set wgbsA.all_dac_tests_1gy_compare;
   if CMH_FDR_P_1Gy_0gy_v2 ne . and CMH_FDR_P_1Gy_0gy_v2 < 0.05 then flag_CMH2_FDR_P05=1;
   else  flag_CMH2_FDR_P05=0;
   if CMH_FDR_P_1Gy_0gy_v2 ne . and CMH_FDR_P_1Gy_0gy_v2 < 0.01 then flag_CMH2_FDR_P01=1;
   else  flag_CMH2_FDR_P01=0;

   if BD_FDR_P_1Gy_0gy ne . and BD_FDR_P_1Gy_0gy < 0.05 then flag_BD_FDR_P05=1;
   else  flag_BD_FDR_P05=0;
   if BD_FDR_P_1Gy_0gy ne . and BD_FDR_P_1Gy_0gy < 0.01 then flag_BD_FDR_P01=1;
   else  flag_BD_FDR_P01=0;

   if BD_FDR_P_1Gy_0gy_v2 ne . and BD_FDR_P_1Gy_0gy_v2 < 0.05 then flag_BD2_FDR_P05=1;
   else  flag_BD_FDR_P05=0;
   if BD_FDR_P_1Gy_0gy_v2 ne . and BD_FDR_P_1Gy_0gy_v2 < 0.01 then flag_BD2_FDR_P01=1;
   else  flag_BD_FDR_P01=0;
   
  
if flag_FET_FDR01=. then flag_FET_FDR01=0;
if flag_FET_FDR05=. then flag_FET_FDR05=0;
if flag_binomial_FDR01=. then flag_binomial_FDR01=0;
if flag_binomial_FDR05=. then flag_binomial_FDR05=0;
if flag_logit_FDR01=. then flag_logit_FDR01=0;
if flag_logit_FDR05=. then flag_logit_FDR05=0;
if flag_methylkit_q05=. then flag_methylkit_q05=0;
if flag_metilene_q05=. then flag_metilene_q05=0;
if flag_BD2_FDR_P01=. then flag_BD2_FDR_P01=0;
if flag_BD2_FDR_P05=. then flag_BD2_FDR_P05=0;
if flag_BD_FDR_P01=. then flag_BD_FDR_P01=0;
if flag_BD_FDR_P05=. then flag_BD_FDR_P05=0;
if flag_CMH2_FDR_P01=. then flag_CMH2_FDR_P01=0;
if flag_CMH2_FDR_P05=. then flag_CMH2_FDR_P05=0;
if flag_CMH_FDR_P01=. then flag_CMH_FDR_P01=0;
if flag_CMH_FDR_P05 =. then flag_CMH_FDR_P05 =0;
if flag_FET_either_fdr01=. then flag_FET_either_fdr01=0;
if flag_FET_either_fdr05=. then flag_FET_either_fdr05=0;
if flag_dss_int_fdr05=. then flag_dss_int_fdr05=0;
if flag_dss_trt_fdr05=. then flag_dss_trt_fdr05=0;
if flag_methylkit_q05=. then flag_methylkit_q05=0;
if flag_metilene_noneg_q05=. then flag_metilene_noneg_q05=0;
if flag_metilene_q05=. then flag_metilene_q05=0;

run;




/* DMC counts */


proc contents data=all_all_dmc_tests_1gy;
run;

proc freq data=all_all_dmc_tests_1gy;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc
flag_methdiff_gt_50perc flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
tables site_type * flag_FET_FDR05 * (flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc
flag_methdiff_gt_50perc flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc freq data=all_all_dmc_tests_1gy;
 where flag_methdiff_gt_10perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc freq data=all_all_dmc_tests_1gy;
 where flag_methdiff_gt_20perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);

run;

proc freq data=all_all_dmc_tests_1gy;
 where flag_methdiff_gt_30perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc freq data=all_all_dmc_tests_1gy;
 where flag_methdiff_gt_50perc = 1;
 tables site_type * (flag_FET_FDR01 flag_FET_FDR05 flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
 tables site_type * flag_FET_FDR05 * ( flag_binomial_FDR01 flag_binomial_FDR05
flag_dss_fdr05 flag_logit_FDR01 flag_logit_FDR05 flag_methylkit_q05 flag_methylsig_bb_fdr05
flag_methylsig_bin_fdr05 flag_metilene_q05);
run;

proc contents data=all_all_dac_tests_1gy;
run;



proc freq data=all_all_dac_tests_1gy;
 tables 
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc
flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05 flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05 
flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05
flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05;
run;

proc freq data=all_all_dac_tests_1gy;
 where flag_FET_either_fdr05 = 1;
 tables 
flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc
flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05 flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05 
flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05
flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05;
run;


proc freq data=all_all_dac_tests_1gy;
 tables (flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc) * (flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05 flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05  flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05 flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05);
run;


proc freq data=all_all_dac_tests_1gy;
	where flag_FET_either_fdr05=1;
	 tables (flag_methdiff_gt_10perc flag_methdiff_gt_20perc flag_methdiff_gt_30perc flag_methdiff_gt_50perc) * (flag_BD2_FDR_P01 flag_BD2_FDR_P05 flag_BD_FDR_P01 flag_BD_FDR_P05 flag_CMH2_FDR_P01 flag_CMH2_FDR_P05 flag_CMH_FDR_P01 flag_CMH_FDR_P05  flag_binomial_FDR01 flag_binomial_FDR05 flag_logit_FDR01 flag_logit_FDR05
	 flag_FET_either_fdr01 flag_FET_either_fdr05 flag_dss_int_fdr05 flag_dss_trt_fdr05 flag_methylkit_q05 flag_metilene_noneg_q05 flag_metilene_q05);
 run;




