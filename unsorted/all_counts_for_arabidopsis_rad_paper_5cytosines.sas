/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


/* RNAseq counts */


/* by time and dose, up and down, include DD. Also do for logFC>1 */

data de_results;
  set arabRNA.arab_results_by_gene;
  keep gene_id flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
  flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
  flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
  mean_cpm_: 
  fdr_: 
  flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_3h_fdr05
  flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_72h_fdr05
  flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_3h_fdr05
  flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_72h_fdr05 ;
run;

proc contents data=de_results;run;quit;


data de_results2;
  set de_results;
  log2fc_01gy_Mock_1h = log2(mean_cpm_01gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_01gy_Mock_3h = log2(mean_cpm_01gy_3h ) - log2(mean_cpm_Mock_3h );
  log2fc_01gy_Mock_24h = log2(mean_cpm_01gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_01gy_Mock_72h = log2(mean_cpm_01gy_72h ) - log2(mean_cpm_Mock_72h );
  log2fc_1gy_Mock_1h = log2(mean_cpm_1gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_1gy_Mock_3h = log2(mean_cpm_1gy_3h) - log2(mean_cpm_Mock_3h );
  log2fc_1gy_Mock_24h = log2(mean_cpm_1gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_1gy_Mock_72h = log2(mean_cpm_1gy_72h ) - log2(mean_cpm_Mock_72h );
run;

data up_01_1 up_01_3 up_01_24 up_01_72
     dn_01_1 dn_01_3 dn_01_24 dn_01_72
     up_1_1  up_1_3  up_1_24  up_1_72
     dn_1_1  dn_1_3  dn_1_24  dn_1_72;
     set de_Results2;
     if (flag_01gy_v_mock_1h_fdr05=1 or (flag_01gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1;
     if (flag_01gy_v_mock_1h_fdr05=1 or (flag_01gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1;
     if (flag_01gy_v_mock_3h_fdr05=1 or (flag_01gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3;
     if (flag_01gy_v_mock_3h_fdr05=1 or (flag_01gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3;
     if (flag_01gy_v_mock_24h_fdr05=1 or (flag_01gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24;
     if (flag_01gy_v_mock_24h_fdr05=1 or (flag_01gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24;
     if (flag_01gy_v_mock_72h_fdr05=1 or (flag_01gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72;
     if (flag_01gy_v_mock_72h_fdr05=1 or (flag_01gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72;

     if (flag_1gy_v_mock_1h_fdr05=1 or (flag_1gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1;
     if (flag_1gy_v_mock_1h_fdr05=1 or (flag_1gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1;
     if (flag_1gy_v_mock_3h_fdr05=1 or (flag_1gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3;
     if (flag_1gy_v_mock_3h_fdr05=1 or (flag_1gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3;
     if (flag_1gy_v_mock_24h_fdr05=1 or (flag_1gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24;
     if (flag_1gy_v_mock_24h_fdr05=1 or (flag_1gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24;
     if (flag_1gy_v_mock_72h_fdr05=1 or (flag_1gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72;
     if (flag_1gy_v_mock_72h_fdr05=1 or (flag_1gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72;
keep gene_id;
run;


/*

There were 33977 observations read from the data set WORK.DE_RESULTS2.
The data set WORK.UP_01_1 has 3092 observations and 1 variables.
The data set WORK.UP_01_3 has 643 observations and 1 variables.
The data set WORK.UP_01_24 has 1003 observations and 1 variables.
The data set WORK.UP_01_72 has 1341 observations and 1 variables.
The data set WORK.DN_01_1 has 2674 observations and 1 variables.
The data set WORK.DN_01_3 has 571 observations and 1 variables.
The data set WORK.DN_01_24 has 642 observations and 1 variables.
The data set WORK.DN_01_72 has 1309 observations and 1 variables.
The data set WORK.UP_1_1 has 2462 observations and 1 variables.
The data set WORK.UP_1_3 has 2382 observations and 1 variables.
The data set WORK.UP_1_24 has 633 observations and 1 variables.
The data set WORK.UP_1_72 has 855 observations and 1 variables.
The data set WORK.DN_1_1 has 1813 observations and 1 variables.
The data set WORK.DN_1_3 has 821 observations and 1 variables.
The data set WORK.DN_1_24 has 501 observations and 1 variables.
The data set WORK.DN_1_72 has 463 observations and 1 variables.
*/


data up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
     dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1
     up_1_1_fc1  up_1_3_fc1  up_1_24_fc1  up_1_72_fc1
     dn_1_1_fc1  dn_1_3_fc1  dn_1_24_fc1  dn_1_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;
keep gene_id;
run;


/*

The data set WORK.UP_01_1_FC1 has 389 observations and 1 variables.
The data set WORK.UP_01_3_FC1 has 22 observations and 1 variables.
The data set WORK.UP_01_24_FC1 has 24 observations and 1 variables.
The data set WORK.UP_01_72_FC1 has 51 observations and 1 variables.
The data set WORK.DN_01_1_FC1 has 109 observations and 1 variables.
The data set WORK.DN_01_3_FC1 has 11 observations and 1 variables.
The data set WORK.DN_01_24_FC1 has 8 observations and 1 variables.
The data set WORK.DN_01_72_FC1 has 15 observations and 1 variables.
The data set WORK.UP_1_1_FC1 has 274 observations and 1 variables.
The data set WORK.UP_1_3_FC1 has 138 observations and 1 variables.
The data set WORK.UP_1_24_FC1 has 0 observations and 1 variables.
The data set WORK.UP_1_72_FC1 has 39 observations and 1 variables.
The data set WORK.DN_1_1_FC1 has 69 observations and 1 variables.
The data set WORK.DN_1_3_FC1 has 28 observations and 1 variables.
The data set WORK.DN_1_24_FC1 has 5 observations and 1 variables.
The data set WORK.DN_1_72_FC1 has 5 observations and 1 variables.
*/

proc sort data=up_01_1; by gene_id;     proc sort data=up_01_1_fc1; by gene_id;
proc sort data=up_01_3; by gene_id;     proc sort data=up_01_3_fc1; by gene_id;
proc sort data=up_01_24; by gene_id;     proc sort data=up_01_24_fc1; by gene_id;
proc sort data=up_01_72; by gene_id;     proc sort data=up_01_72_fc1; by gene_id;

proc sort data=dn_01_1; by gene_id;     proc sort data=dn_01_1_fc1; by gene_id;
proc sort data=dn_01_3; by gene_id;     proc sort data=dn_01_3_fc1; by gene_id;
proc sort data=dn_01_24; by gene_id;     proc sort data=dn_01_24_fc1; by gene_id;
proc sort data=dn_01_72; by gene_id;     proc sort data=dn_01_72_fc1; by gene_id;

proc sort data=up_1_1; by gene_id;     proc sort data=up_1_1_fc1; by gene_id;
proc sort data=up_1_3; by gene_id;     proc sort data=up_1_3_fc1; by gene_id;
proc sort data=up_1_24; by gene_id;     proc sort data=up_1_24_fc1; by gene_id;
proc sort data=up_1_72; by gene_id;     proc sort data=up_1_72_fc1; by gene_id;

proc sort data=dn_1_1; by gene_id;     proc sort data=dn_1_1_fc1; by gene_id;
proc sort data=dn_1_3; by gene_id;     proc sort data=dn_1_3_fc1; by gene_id;
proc sort data=dn_1_24; by gene_id;     proc sort data=dn_1_24_fc1; by gene_id;
proc sort data=dn_1_72; by gene_id;     proc sort data=dn_1_72_fc1; by gene_id;
run;





data deg_compare_any_updn;
  merge up_01_1_fc1 (in=in1) dn_01_1_fc1 (in=in2) up_01_3_fc1 (in=in3) dn_01_3_fc1 (in=in4)
        up_01_24_fc1 (in=in5) dn_01_24_fc1 (in=in6) up_01_72_fc1 (in=in7) dn_01_72_fc1 (in=in8)
        up_1_1_fc1 (in=in9) dn_1_1_fc1 (in=in10) up_1_3_fc1 (in=in11) dn_1_3_fc1 (in=in12)
        up_1_24_fc1 (in=in13) dn_1_24_fc1 (in=in14) up_1_72_fc1 (in=in15) dn_1_72_fc1 (in=in16);
  by gene_id;
  if in1 or in2 then flag_DEG_01_1=1; else  flag_DEG_01_1=0;
  if in3 or in4 then flag_DEG_01_3=1; else  flag_DEG_01_3=0;
  if in5 or in6 then flag_DEG_01_24=1; else  flag_DEG_01_24=0;
  if in7 or in8 then flag_DEG_01_72=1; else  flag_DEG_01_72=0;

  if in9 or in10 then flag_DEG_1_1=1; else  flag_DEG_1_1=0;
  if in11 or in12 then flag_DEG_1_3=1; else  flag_DEG_1_3=0;
  if in13 or in14 then flag_DEG_1_24=1; else  flag_DEG_1_24=0;
  if in15 or in16 then flag_DEG_1_72=1; else  flag_DEG_1_72=0;
run;
proc freq data=deg_compare_any_updn noprint;
 tables flag_DEG_01_1*flag_DEG_1_1*flag_DEG_01_72*flag_DEG_1_72 / out=ctabs_early_late;
run;
proc print data=ctabs_early_late;
run;


proc freq data=deg_compare_any_updn ;
 tables flag_DEG_01_1*flag_DEG_1_1
        flag_DEG_01_3*flag_DEG_1_3
        flag_DEG_01_24*flag_DEG_1_24
       flag_DEG_01_72*flag_DEG_1_72;
run;



data deg_compare_any;
  merge up_01_1 (in=in1) dn_01_1 (in=in2) up_01_3 (in=in3) dn_01_3 (in=in4)
        up_01_24 (in=in5) dn_01_24 (in=in6) up_01_72 (in=in7) dn_01_72 (in=in8)
        up_1_1 (in=in9) dn_1_1 (in=in10) up_1_3 (in=in11) dn_1_3 (in=in12)
        up_1_24 (in=in13) dn_1_24 (in=in14) up_1_72 (in=in15) dn_1_72 (in=in16);
  by gene_id;
  if in1 then flag_DEG_up_01_1=1; else  flag_DEG_up_01_1=0;
  if in2 then flag_DEG_dn_01_1=1; else  flag_DEG_dn_01_1=0;
  if in3 then flag_DEG_up_01_3=1; else  flag_DEG_up_01_3=0;
  if in4 then flag_DEG_dn_01_3=1; else  flag_DEG_dn_01_3=0;
  if in5 then flag_DEG_up_01_24=1; else  flag_DEG_up_01_24=0;
  if in6 then flag_DEG_dn_01_24=1; else  flag_DEG_dn_01_24=0;
  if in7 then flag_DEG_up_01_72=1; else  flag_DEG_up_01_72=0;
  if in8 then flag_DEG_dn_01_72=1; else  flag_DEG_dn_01_72=0;

  if in9 then flag_DEG_up_1_1=1; else  flag_DEG_up_1_1=0;
  if in10 then flag_DEG_dn_1_1=1; else  flag_DEG_dn_1_1=0;
  if in11 then flag_DEG_up_1_3=1; else  flag_DEG_up_1_3=0;
  if in12 then flag_DEG_dn_1_3=1; else  flag_DEG_dn_1_3=0;
  if in13 then flag_DEG_up_1_24=1; else  flag_DEG_up_1_24=0;
  if in14 then flag_DEG_dn_1_24=1; else  flag_DEG_dn_1_24=0;
  if in15 then flag_DEG_up_1_72=1; else  flag_DEG_up_1_72=0;
  if in16 then flag_DEG_dn_1_72=1; else  flag_DEG_dn_1_72=0;
run;


data deg_compare_fc1;
  merge up_01_1_fc1 (in=in1) dn_01_1_fc1 (in=in2) up_01_3_fc1 (in=in3) dn_01_3_fc1 (in=in4)
        up_01_24_fc1 (in=in5) dn_01_24_fc1 (in=in6) up_01_72_fc1 (in=in7) dn_01_72_fc1 (in=in8)
        up_1_1_fc1 (in=in9) dn_1_1_fc1 (in=in10) up_1_3_fc1 (in=in11) dn_1_3_fc1 (in=in12)
        up_1_24_fc1 (in=in13) dn_1_24_fc1 (in=in14) up_1_72_fc1 (in=in15) dn_1_72_fc1 (in=in16);
  by gene_id;
  if in1 then flag_DEG_up_01_1=1; else  flag_DEG_up_01_1=0;
  if in2 then flag_DEG_dn_01_1=1; else  flag_DEG_dn_01_1=0;
  if in3 then flag_DEG_up_01_3=1; else  flag_DEG_up_01_3=0;
  if in4 then flag_DEG_dn_01_3=1; else  flag_DEG_dn_01_3=0;
  if in5 then flag_DEG_up_01_24=1; else  flag_DEG_up_01_24=0;
  if in6 then flag_DEG_dn_01_24=1; else  flag_DEG_dn_01_24=0;
  if in7 then flag_DEG_up_01_72=1; else  flag_DEG_up_01_72=0;
  if in8 then flag_DEG_dn_01_72=1; else  flag_DEG_dn_01_72=0;

  if in9 then flag_DEG_up_1_1=1; else  flag_DEG_up_1_1=0;
  if in10 then flag_DEG_dn_1_1=1; else  flag_DEG_dn_1_1=0;
  if in11 then flag_DEG_up_1_3=1; else  flag_DEG_up_1_3=0;
  if in12 then flag_DEG_dn_1_3=1; else  flag_DEG_dn_1_3=0;
  if in13 then flag_DEG_up_1_24=1; else  flag_DEG_up_1_24=0;
  if in14 then flag_DEG_dn_1_24=1; else  flag_DEG_dn_1_24=0;
  if in15 then flag_DEG_up_1_72=1; else  flag_DEG_up_1_72=0;
  if in16 then flag_DEG_dn_1_72=1; else  flag_DEG_dn_1_72=0;
run;

proc freq data=deg_compare_any noprint;
   tables flag_DEG_up_01_1*flag_DEG_up_01_3*flag_DEG_up_01_24*flag_DEG_up_01_72 / out=ctabs_01_up;
   tables flag_DEG_up_1_1*flag_DEG_up_1_3*flag_DEG_up_1_24*flag_DEG_up_1_72 / out=ctabs_1_up;
   tables flag_DEG_dn_01_1*flag_DEG_dn_01_3*flag_DEG_dn_01_24*flag_DEG_dn_01_72 / out=ctabs_01_dn;
   tables flag_DEG_dn_1_1*flag_DEG_dn_1_3*flag_DEG_dn_1_24*flag_DEG_dn_1_72 / out=ctabs_1_dn;
   tables flag_DEG_up_01_1*flag_DEG_up_1_1 / out=ctabs_1h_up;
   tables flag_DEG_dn_01_1*flag_DEG_dn_1_1 / out=ctabs_1h_dn;
   tables flag_DEG_up_01_3*flag_DEG_up_1_3 / out=ctabs_3h_up;
   tables flag_DEG_dn_01_3*flag_DEG_dn_1_3 / out=ctabs_3h_dn;
   tables flag_DEG_up_01_24*flag_DEG_up_1_24 / out=ctabs_24h_up;
   tables flag_DEG_dn_01_24*flag_DEG_dn_1_24 / out=ctabs_24h_dn;
   tables flag_DEG_up_01_72*flag_DEG_up_1_72 / out=ctabs_72h_up;
   tables flag_DEG_dn_01_72*flag_DEG_dn_1_72 / out=ctabs_72h_dn;
   tables flag_DEG_up_01_1*flag_DEG_up_1_1*flag_DEG_up_01_72*flag_DEG_up_1_72 / out=ctabs_early_v_late_up;
   tables flag_DEG_dn_01_1*flag_DEG_dn_1_1*flag_DEG_dn_01_72*flag_DEG_dn_1_72 / out=ctabs_early_v_late_dn;
run;

proc print data=ctabs_01_up; run;
/*
  flag_      flag_      flag_      flag_
 DEG_up_    DEG_up_    DEG_up_    DEG_up_
   01_1       01_3      01_24      01_72     COUNT    PERCENT

    0          0          0          0        5915    53.7923
    0          0          0          1         836     7.6028
    0          0          1          0         558     5.0746
    0          0          1          1          91     0.8276
    0          1          0          0         350     3.1830
    0          1          0          1          71     0.6457
    0          1          1          0          69     0.6275
    0          1          1          1          14     0.1273
    1          0          0          0        2462    22.3900
    1          0          0          1         243     2.2099
    1          0          1          0         198     1.8007
    1          0          1          1          50     0.4547
    1          1          0          0          88     0.8003
    1          1          0          1          28     0.2546
    1          1          1          0          15     0.1364
    1          1          1          1           8     0.0728

*/
proc print data=ctabs_1_up; run;
/*
  flag_      flag_      flag_      flag_
 DEG_up_    DEG_up_    DEG_up_    DEG_up_
   01_1       01_3      01_24      01_72     COUNT    PERCENT

    0          0          0          0        5915    53.7923
    0          0          0          1         836     7.6028
    0          0          1          0         558     5.0746
    0          0          1          1          91     0.8276
    0          1          0          0         350     3.1830
    0          1          0          1          71     0.6457
    0          1          1          0          69     0.6275
    0          1          1          1          14     0.1273
    1          0          0          0        2462    22.3900
    1          0          0          1         243     2.2099
    1          0          1          0         198     1.8007
    1          0          1          1          50     0.4547
    1          1          0          0          88     0.8003
    1          1          0          1          28     0.2546
    1          1          1          0          15     0.1364
    1          1          1          1           8     0.0728

*/
proc print data=ctabs_01_dn; run;
/*
  flag_      flag_      flag_      flag_
 DEG_dn_    DEG_dn_    DEG_dn_    DEG_dn_
   01_1       01_3      01_24      01_72     COUNT    PERCENT

    0          0          0          0        6429    58.4667
    0          0          0          1         992     9.0215
    0          0          1          0         366     3.3285
    0          0          1          1          76     0.6912
    0          1          0          0         333     3.0284
    0          1          0          1          57     0.5184
    0          1          1          0          58     0.5275
    0          1          1          1          11     0.1000
    1          0          0          0        2322    21.1168
    1          0          0          1         126     1.1459
    1          0          1          0          88     0.8003
    1          0          1          1          26     0.2364
    1          1          0          0          79     0.7184
    1          1          0          1          16     0.1455
    1          1          1          0          12     0.1091
    1          1          1          1           5     0.0455

*/
proc print data=ctabs_1_dn; run;
/*
  flag_      flag_      flag_      flag_
 DEG_dn_    DEG_dn_    DEG_dn_    DEG_dn_
   1_1        1_3        1_24       1_72     COUNT    PERCENT

    0          0          0          0        7836    71.2623
    0          0          0          1         296     2.6919
    0          0          1          0         317     2.8829
    0          0          1          1          41     0.3729
    0          1          0          0         581     5.2837
    0          1          0          1          47     0.4274
    0          1          1          0          51     0.4638
    0          1          1          1          14     0.1273
    1          0          0          0        1578    14.3507
    1          0          0          1          44     0.4001
    1          0          1          0          52     0.4729
    1          0          1          1          11     0.1000
    1          1          0          0         104     0.9458
    1          1          0          1           9     0.0818
    1          1          1          0          14     0.1273
    1          1          1          1           1     0.0091

*/
proc print data=ctabs_1h_up; run;
/*
 flag_      flag_
DEG_up_    DEG_up_
  01_1       1_1      COUNT    PERCENT

   0          0        7229    65.7421
   0          1         675     6.1386
   1          0        1305    11.8680
   1          1        1787    16.2514

*/
proc print data=ctabs_1h_dn; run;
/*
  flag_      flag_
 DEG_dn_    DEG_dn_
   01_1       1_1      COUNT    PERCENT

    0          0        7848    71.3714
    0          1         474     4.3107
    1          0        1335    12.1408
    1          1        1339    12.1772

*/
proc print data=ctabs_3h_up; run;
/*
  flag_      flag_
 DEG_up_    DEG_up_
   01_3       1_3      COUNT    PERCENT

    0          0        8275    75.2546
    0          1        2078    18.8978
    1          0         339     3.0829
    1          1         304     2.7646

*/
proc print data=ctabs_3h_dn; run;
/*
   flag_      flag_
  DEG_dn_    DEG_dn_
    01_3       1_3      COUNT    PERCENT

     0          0        9880    89.8509
     0          1         545     4.9563
     1          0         295     2.6828
     1          1         276     2.5100

           The SAS System
*/
proc print data=ctabs_24h_up; run;
/*

    flag_      flag_
   DEG_up_    DEG_up_
    01_24       1_24     COUNT    PERCENT

      0          0        9636    87.6319
      0          1         357     3.2466
      1          0         727     6.6115
      1          1         276     2.5100
*/
proc print data=ctabs_24h_dn; run;
/*
  flag_      flag_
 DEG_dn_    DEG_dn_
  01_24       1_24     COUNT    PERCENT

    0          0       10098    91.8334
    0          1         256     2.3281
    1          0         397     3.6104
    1          1         245     2.2281
*/
proc print data=ctabs_72h_up; run;
/*
  flag_      flag_
 DEG_up_    DEG_up_
  01_72       1_72     COUNT    PERCENT

    0          0        9117    82.9120
    0          1         538     4.8927
    1          0        1024     9.3125
    1          1         317     2.8829

*/
proc print data=ctabs_72h_dn; run;
/*
    flag_      flag_
   DEG_dn_    DEG_dn_
    01_72       1_72     COUNT    PERCENT

      0          0        9469    86.1131
      0          1         218     1.9825
      1          0        1064     9.6762
      1          1         245     2.2281
*/
proc print data=ctabs_early_v_late_up; run;
/*
  flag_      flag_      flag_      flag_
 DEG_up_    DEG_up_    DEG_up_    DEG_up_
   01_1       1_1       01_72       1_72     COUNT    PERCENT

    0          0          0          0        5887    53.5377
    0          0          0          1         395     3.5922
    0          0          1          0         731     6.6479
    0          0          1          1         216     1.9644
    0          1          0          0         558     5.0746
    0          1          0          1          52     0.4729
    0          1          1          0          41     0.3729
    0          1          1          1          24     0.2183
    1          0          0          0        1067     9.7035
    1          0          0          1          62     0.5638
    1          0          1          0         131     1.1913
    1          0          1          1          45     0.4092
    1          1          0          0        1605    14.5962
    1          1          0          1          29     0.2637
    1          1          1          0         121     1.1004
    1          1          1          1          32     0.2910

*/
proc print data=ctabs_early_v_late_dn; run;
/*
  flag_      flag_      flag_      flag_
 DEG_dn_    DEG_dn_    DEG_dn_    DEG_dn_
   01_1       1_1       01_72       1_72     COUNT    PERCENT

    0          0          0          0        6603    60.0491
    0          0          0          1         154     1.4005
    0          0          1          0         888     8.0757
    0          0          1          1         203     1.8461
    0          1          0          0         413     3.7559
    0          1          0          1          16     0.1455
    0          1          1          0          33     0.3001
    0          1          1          1          12     0.1091
    1          0          0          0        1223    11.1222
    1          0          0          1          23     0.2092
    1          0          1          0          71     0.6457
    1          0          1          1          18     0.1637
    1          1          0          0        1230    11.1859
    1          1          0          1          25     0.2274
    1          1          1          0          72     0.6548
    1          1          1          1          12     0.1091


*/


proc freq data=deg_compare_fc1 noprint;
   tables flag_DEG_up_01_1*flag_DEG_up_01_3*flag_DEG_up_01_24*flag_DEG_up_01_72 / out=ctabs_01_up;
   tables flag_DEG_up_1_1*flag_DEG_up_1_3*flag_DEG_up_1_24*flag_DEG_up_1_72 / out=ctabs_1_up;
   tables flag_DEG_dn_01_1*flag_DEG_dn_01_3*flag_DEG_dn_01_24*flag_DEG_dn_01_72 / out=ctabs_01_dn;
   tables flag_DEG_dn_1_1*flag_DEG_dn_1_3*flag_DEG_dn_1_24*flag_DEG_dn_1_72 / out=ctabs_1_dn;
   tables flag_DEG_up_01_1*flag_DEG_up_1_1 / out=ctabs_1h_up;
   tables flag_DEG_dn_01_1*flag_DEG_dn_1_1 / out=ctabs_1h_dn;
   tables flag_DEG_up_01_3*flag_DEG_up_1_3 / out=ctabs_3h_up;
   tables flag_DEG_dn_01_3*flag_DEG_dn_1_3 / out=ctabs_3h_dn;
   tables flag_DEG_up_01_24*flag_DEG_up_1_24 / out=ctabs_24h_up;
   tables flag_DEG_dn_01_24*flag_DEG_dn_1_24 / out=ctabs_24h_dn;
   tables flag_DEG_up_01_72*flag_DEG_up_1_72 / out=ctabs_72h_up;
   tables flag_DEG_dn_01_72*flag_DEG_dn_1_72 / out=ctabs_72h_dn;
   tables flag_DEG_up_01_1*flag_DEG_up_1_1*flag_DEG_up_01_72*flag_DEG_up_1_72 / out=ctabs_early_v_late_up;
   tables flag_DEG_dn_01_1*flag_DEG_dn_1_1*flag_DEG_dn_01_72*flag_DEG_dn_1_72 / out=ctabs_early_v_late_dn;
run;

proc print data=ctabs_01_up; run;
/*

   flag_      flag_      flag_      flag_
  DEG_up_    DEG_up_    DEG_up_    DEG_up_
    01_1       01_3      01_24      01_72     COUNT    PERCENT

     0          0          0          0        388     45.5399
     0          0          0          1         43      5.0469
     0          0          1          0         14      1.6432
     0          0          1          1          1      0.1174
     0          1          0          0         15      1.7606
     0          1          0          1          1      0.1174
     0          1          1          1          1      0.1174
     1          0          0          0        375     44.0141
     1          0          0          1          2      0.2347
     1          0          1          0          6      0.7042
     1          0          1          1          1      0.1174
     1          1          0          0          3      0.3521
     1          1          0          1          1      0.1174
     1          1          1          1          1      0.1174

*/
proc print data=ctabs_1_up; run;
/*
 flag_      flag_      flag_      flag_
DEG_up_    DEG_up_    DEG_up_    DEG_up_
  1_1        1_3        1_24       1_72     COUNT    PERCENT

   0          0          0          0        420     49.2958
   0          0          0          1         34      3.9906
   0          1          0          0        121     14.2019
   0          1          0          1          3      0.3521
   1          0          0          0        259     30.3991
   1          0          0          1          1      0.1174
   1          1          0          0         13      1.5258
   1          1          0          1          1      0.1174

*/
proc print data=ctabs_01_dn; run;
/*
 flag_      flag_      flag_      flag_
DEG_dn_    DEG_dn_    DEG_dn_    DEG_dn_
  01_1       01_3      01_24      01_72     COUNT    PERCENT

   0          0          0          0        713     83.6854
   0          0          0          1         11      1.2911
   0          0          1          0          6      0.7042
   0          0          1          1          2      0.2347
   0          1          0          0          9      1.0563
   0          1          0          1          2      0.2347
   1          0          0          0        109     12.7934
*/
proc print data=ctabs_1_dn; run;
/*
  flag_      flag_      flag_      flag_
 DEG_dn_    DEG_dn_    DEG_dn_    DEG_dn_
   1_1        1_3        1_24       1_72     COUNT    PERCENT

    0          0          0          0        748     87.7934
    0          0          0          1          2      0.2347
    0          0          1          0          5      0.5869
    0          1          0          0         25      2.9343
    0          1          0          1          3      0.3521
    1          0          0          0         69      8.0986

*/
proc print data=ctabs_1h_up; run;
/*
  flag_      flag_
 DEG_up_    DEG_up_
   01_1       1_1      COUNT    PERCENT

    0          0        403     47.3005
    0          1         60      7.0423
    1          0        175     20.5399
    1          1        214     25.1174

*/
proc print data=ctabs_1h_dn; run;
/*
   flag_      flag_
  DEG_dn_    DEG_dn_
    01_1       1_1      COUNT    PERCENT

     0          0        710     83.3333
     0          1         33      3.8732
     1          0         73      8.5681
     1          1         36      4.2254
*/
proc print data=ctabs_3h_up; run;
/*
   flag_      flag_
  DEG_up_    DEG_up_
    01_3       1_3      COUNT    PERCENT

     0          0        703     82.5117
     0          1        127     14.9061
     1          0         11      1.2911
     1          1         11      1.2911

*/
proc print data=ctabs_3h_dn; run;
/*
   flag_      flag_
  DEG_dn_    DEG_dn_
    01_3       1_3      COUNT    PERCENT

     0          0        817     95.8920
     0          1         24      2.8169
     1          0          7      0.8216
     1          1          4      0.4695
*/
proc print data=ctabs_24h_up; run;
/*
   flag_      flag_
  DEG_up_    DEG_up_
   01_24       1_24     COUNT    PERCENT

     0          0        828     97.1831
     1          0         24      2.8169

*/
proc print data=ctabs_24h_dn; run;
/*
  flag_      flag_
 DEG_dn_    DEG_dn_
  01_24       1_24     COUNT    PERCENT

    0          0        839     98.4742
    0          1          5      0.5869
    1          0          8      0.9390

*/
proc print data=ctabs_72h_up; run;
/*
  flag_      flag_
 DEG_up_    DEG_up_
  01_72       1_72     COUNT    PERCENT

    0          0        767     90.0235
    0          1         34      3.9906
    1          0         46      5.3991
    1          1          5      0.5869
*/
proc print data=ctabs_72h_dn; run;
/*
   flag_      flag_
  DEG_dn_    DEG_dn_
   01_72       1_72     COUNT    PERCENT

     0          0        836     98.1221
     0          1          1      0.1174
     1          0         11      1.2911
     1          1          4      0.4695

*/
proc print data=ctabs_early_v_late_up; run;
/*
   flag_      flag_      flag_      flag_
  DEG_up_    DEG_up_    DEG_up_    DEG_up_
    01_1       1_1       01_72       1_72     COUNT    PERCENT

     0          0          0          0        325     38.1455
     0          0          0          1         32      3.7559
     0          0          1          0         42      4.9296
     0          0          1          1          4      0.4695
     0          1          0          0         60      7.0423
     1          0          0          0        171     20.0704
     1          0          1          0          3      0.3521
     1          0          1          1          1      0.1174
     1          1          0          0        211     24.7653
     1          1          0          1          2      0.2347
     1          1          1          0          1      0.1174
*/
proc print data=ctabs_early_v_late_dn; run;
/*
   flag_      flag_      flag_      flag_
  DEG_dn_    DEG_dn_    DEG_dn_    DEG_dn_
    01_1       1_1       01_72       1_72     COUNT    PERCENT

     0          0          0          0        695     81.5728
     0          0          0          1          1      0.1174
     0          0          1          0         10      1.1737
     0          0          1          1          4      0.4695
     0          1          0          0         32      3.7559
     0          1          1          0          1      0.1174
     1          0          0          0         73      8.5681
     1          1          0          0         36      4.2254


*/

/* Output data for heatmaps:
    (1) logCPM
    (2) logFC */

data de_genes;
  set deg_compare_any;
  keep gene_id;
run;

data de_genes_fc1;
  set deg_compare_fc1;
  keep gene_id;
run;


data means;
  set arabRNA.arab_gene_mean_cpm_by_trt_time;
  log2_cpm_mock_1h= log2( mean_cpm_mock_1h + 1);
  log2_cpm_mock_3h= log2( mean_cpm_mock_3h + 1);
  log2_cpm_mock_24h= log2( mean_cpm_mock_24h + 1);
  log2_cpm_mock_72h= log2( mean_cpm_mock_72h + 1);

  log2_cpm_01gy_1h= log2( mean_cpm_01gy_1h + 1);
  log2_cpm_01gy_3h= log2( mean_cpm_01gy_3h + 1);
  log2_cpm_01gy_24h= log2( mean_cpm_01gy_24h + 1);
  log2_cpm_01gy_72h= log2( mean_cpm_01gy_72h + 1);

  log2_cpm_1gy_1h= log2( mean_cpm_1gy_1h + 1);
  log2_cpm_1gy_3h= log2( mean_cpm_1gy_3h + 1);
  log2_cpm_1gy_24h= log2( mean_cpm_1gy_24h + 1);
  log2_cpm_1gy_72h= log2( mean_cpm_1gy_72h + 1);

  log2fc_01gy_Mock_1h = log2_cpm_01gy_1h - log2_cpm_Mock_1h ;
  log2fc_01gy_Mock_3h = log2_cpm_01gy_3h - log2_cpm_Mock_3h ;
  log2fc_01gy_Mock_24h = log2_cpm_01gy_24h - log2_cpm_Mock_24h ;
  log2fc_01gy_Mock_72h = log2_cpm_01gy_72h - log2_cpm_Mock_72h ;

  log2fc_1gy_Mock_1h = log2_cpm_1gy_1h - log2_cpm_Mock_1h ;
  log2fc_1gy_Mock_3h = log2_cpm_1gy_3h - log2_cpm_Mock_3h ;
  log2fc_1gy_Mock_24h = log2_cpm_1gy_24h - log2_cpm_Mock_24h ;
  log2fc_1gy_Mock_72h = log2_cpm_1gy_72h - log2_cpm_Mock_72h ;

  keep gene_id log2_cpm_: log2fc_: ;
run;


proc corr data=means pearson;
  var log2fc_01gy_Mock_1h log2fc_1gy_Mock_1h
      log2fc_01gy_Mock_3h log2fc_1gy_Mock_3h
      log2fc_01gy_Mock_24h log2fc_1gy_Mock_24h
      log2fc_01gy_Mock_72h log2fc_1gy_Mock_72h;
run;


/*


                               log2fc_        log2fc_       log2fc_        log2fc_         log2fc_        log2fc_         log2fc_        log2fc_
                                 01gy_      1gy_Mock_         01gy_      1gy_Mock_      01gy_Mock_      1gy_Mock_      01gy_Mock_      1gy_Mock_
                               Mock_1h             1h       Mock_3h             3h             24h            24h             72h            72h

    log2fc_01gy_Mock_1h        1.00000        0.84241       0.18017        0.24279         0.17122        0.06078         0.15641        0.24168
                                               <.0001        <.0001         <.0001          <.0001         <.0001          <.0001         <.0001

    log2fc_1gy_Mock_1h         0.84241        1.00000       0.22528        0.25921         0.19751        0.04416         0.14733        0.25025
                                <.0001                       <.0001         <.0001          <.0001         <.0001          <.0001         <.0001

    log2fc_01gy_Mock_3h        0.18017        0.22528       1.00000        0.63554         0.20744       -0.03424         0.16917        0.42631
                                <.0001         <.0001                       <.0001          <.0001         <.0001          <.0001         <.0001

    log2fc_1gy_Mock_3h         0.24279        0.25921       0.63554        1.00000         0.39425        0.02801         0.30870        0.45141
                                <.0001         <.0001        <.0001                         <.0001         <.0001          <.0001         <.0001

                            log2fc_        log2fc_       log2fc_        log2fc_         log2fc_        log2fc_         log2fc_        log2fc_
                              01gy_      1gy_Mock_         01gy_      1gy_Mock_      01gy_Mock_      1gy_Mock_      01gy_Mock_      1gy_Mock_
                            Mock_1h             1h       Mock_3h             3h             24h            24h             72h            72h

 log2fc_01gy_Mock_24h       0.17122        0.19751       0.20744        0.39425         1.00000        0.41631         0.38411        0.45486
                             <.0001         <.0001        <.0001         <.0001                         <.0001          <.0001         <.0001

 log2fc_1gy_Mock_24h        0.06078        0.04416      -0.03424        0.02801         0.41631        1.00000         0.35285        0.23994
                             <.0001         <.0001        <.0001         <.0001          <.0001                         <.0001         <.0001

 log2fc_01gy_Mock_72h       0.15641        0.14733       0.16917        0.30870         0.38411        0.35285         1.00000        0.55631
                             <.0001         <.0001        <.0001         <.0001          <.0001         <.0001                         <.0001

 log2fc_1gy_Mock_72h        0.24168        0.25025       0.42631        0.45141         0.45486        0.23994         0.55631        1.00000
                             <.0001         <.0001        <.0001         <.0001          <.0001         <.0001          <.0001


*/


proc sgplot data=means;
   scatter x=log2fc_01gy_Mock_1h y=log2fc_1gy_Mock_1h;
run;
proc sgplot data=means;
   scatter x=log2fc_01gy_Mock_3h y=log2fc_1gy_Mock_3h;
run;
proc sgplot data=means;
   scatter x=log2fc_01gy_Mock_24h y=log2fc_1gy_Mock_24h;
run;
proc sgplot data=means;
   scatter x=log2fc_01gy_Mock_72h y=log2fc_1gy_Mock_72h;
run;




proc sort data=means;
  by gene_id;
proc sort data=deg_compare_fc1;
  by gene_id;
run;

data means_de;
  merge means (in=in1) deg_compare_fc1 (in=in2);
  by gene_id;
  if in1;
run;


data compare_1h;
  set means_de;
  length DE $12.;
  if (flag_DEG_up_01_1=1 or flag_DEG_dn_01_1=1) and (flag_DEG_up_1_1=1 or flag_DEG_dn_1_1=1) then DE="3-both";
  else if (flag_DEG_up_01_1=1 or flag_DEG_dn_01_1=1) and (flag_DEG_up_1_1=0 and flag_DEG_dn_1_1=0) then DE="2-10cGy";
  else if (flag_DEG_up_01_1=0 and flag_DEG_dn_01_1=0) and (flag_DEG_up_1_1=1 or flag_DEG_dn_1_1=1) then DE="1-100cGy";
  else DE="0-none";
  keep gene_id log2fc_01gy_mock_1h log2fc_1gy_mock_1h DE;
run;


data compare_3h;
  set means_de;
  length DE $12.;
  if (flag_DEG_up_01_3=1 or flag_DEG_dn_01_3=1) and (flag_DEG_up_1_3=1 or flag_DEG_dn_1_3=1) then DE="3-both";
  else if (flag_DEG_up_01_3=1 or flag_DEG_dn_01_3=1) and (flag_DEG_up_1_3=0 and flag_DEG_dn_1_3=0) then DE="2-10cGy";
  else if (flag_DEG_up_01_3=0 and flag_DEG_dn_01_3=0) and (flag_DEG_up_1_3=1 or flag_DEG_dn_1_3=1) then DE="1-100cGy";
  else DE="0-none";
  keep gene_id log2fc_01gy_mock_3h log2fc_1gy_mock_3h DE;
run;


data compare_24h;
  set means_de;
  length DE $12.;
  if (flag_DEG_up_01_24=1 or flag_DEG_dn_01_24=1) and (flag_DEG_up_1_24=1 or flag_DEG_dn_1_24=1) then DE="3-both";
  else if (flag_DEG_up_01_24=1 or flag_DEG_dn_01_24=1) and (flag_DEG_up_1_24=0 and flag_DEG_dn_1_24=0) then DE="2-10cGy";
  else if (flag_DEG_up_01_24=0 and flag_DEG_dn_01_24=0) and (flag_DEG_up_1_24=1 or flag_DEG_dn_1_24=1) then DE="1-100cGy";
  else DE="0-none";
  keep gene_id log2fc_01gy_mock_24h log2fc_1gy_mock_24h DE;
run;


data compare_72h;
  set means_de;
  length DE $12.;
  if (flag_DEG_up_01_72=1 or flag_DEG_dn_01_72=1) and (flag_DEG_up_1_72=1 or flag_DEG_dn_1_72=1) then DE="3-both";
  else if (flag_DEG_up_01_72=1 or flag_DEG_dn_01_72=1) and (flag_DEG_up_1_72=0 and flag_DEG_dn_1_72=0) then DE="2-10cGy";
  else if (flag_DEG_up_01_72=0 and flag_DEG_dn_01_72=0) and (flag_DEG_up_1_72=1 or flag_DEG_dn_1_72=1) then DE="1-100cGy";
  else DE="0-none";
  keep gene_id log2fc_01gy_mock_72h log2fc_1gy_mock_72h DE;
run;


proc sort data=compare_1h; by DE gene_id;
proc sort data=compare_3h; by DE gene_id;
proc sort data=compare_24h; by DE gene_id;
proc sort data=compare_72h; by DE gene_id;
run;

proc export data=compare_1h outfile="!HOME/concannon/DTRA/arabidopsis_logFC_compare_1h.csv" dbms=csv replace; run;
proc export data=compare_3h outfile="!HOME/concannon/DTRA/arabidopsis_logFC_compare_3h.csv" dbms=csv replace; run;
proc export data=compare_24h outfile="!HOME/concannon/DTRA/arabidopsis_logFC_compare_24h.csv" dbms=csv replace; run;
proc export data=compare_72h outfile="!HOME/concannon/DTRA/arabidopsis_logFC_compare_72h.csv" dbms=csv replace; run;




data logfc_diff;
  set means;
  diff_1h=log2fc_01gy_Mock_1h - log2fc_1gy_Mock_1h;
  diff_3h=log2fc_01gy_Mock_3h - log2fc_1gy_Mock_3h;
  diff_24h=log2fc_01gy_Mock_24h - log2fc_1gy_Mock_24h;
  diff_72h=log2fc_01gy_Mock_72h - log2fc_1gy_Mock_72h;
run;

proc means data=logfc_diff noprint;
  var diff_1h;
  output out=diff_1h_distrib mean=mean stddev=sd q1=q1 median=median q3=q3;
run;

proc means data=logfc_diff noprint;
  var diff_3h;
  output out=diff_3h_distrib mean=mean stddev=sd q1=q1 median=median q3=q3;
run;

proc means data=logfc_diff noprint;
  var diff_24h;
  output out=diff_24h_distrib mean=mean stddev=sd q1=q1 median=median q3=q3;
run;

proc means data=logfc_diff noprint;
  var diff_72h;
  output out=diff_72h_distrib mean=mean stddev=sd q1=q1 median=median q3=q3;
run;


proc sort data=de_genes;
  by gene_id;
proc sort data=de_genes_fc1;
  by gene_id;
proc sort data=means;
  by gene_id;
run;


data de_genes_means_fc1;
  merge de_genes_fc1 (in=in1) means (in=in2);
  by gene_id;
  if in1 and in2;
run;




data de_genes_means;
  merge de_genes (in=in1) means (in=in2);
  by gene_id;
  if in1 and in2;
run;

data de_genes_means_fc1;
  merge de_genes_fc1 (in=in1) means (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=de_genes_means;
  by log2fc_01gy_mock_1h log2fc_01gy_mock_3h log2fc_01gy_mock_24h log2fc_01gy_mock_72h;
run;

proc sort data=de_genes_means_fc1;
  by log2fc_01gy_Mock_1h log2fc_01gy_Mock_3h log2fc_01gy_mock_24h log2fc_01gy_mock_72h;
run;


data de_genes_any_logCPM_heatmap;
  retain gene_id log2_cpm_mock_1h log2_cpm_mock_3h log2_cpm_mock_24h log2_cpm_mock_72h
                 log2_cpm_01gy_1h log2_cpm_01gy_3h log2_cpm_01gy_24h log2_cpm_01gy_72h
                 log2_cpm_1gy_1h log2_cpm_1gy_3h log2_cpm_1gy_24h log2_cpm_1gy_72h;
  set de_genes_means;
  keep gene_id log2_cpm_mock_1h log2_cpm_mock_3h log2_cpm_mock_24h log2_cpm_mock_72h
                 log2_cpm_01gy_1h log2_cpm_01gy_3h log2_cpm_01gy_24h log2_cpm_01gy_72h
                 log2_cpm_1gy_1h log2_cpm_1gy_3h log2_cpm_1gy_24h log2_cpm_1gy_72h;
run;

data de_genes_any_logFC_heatmap;
  retain gene_id log2fc_01gy_Mock_1h log2fc_1gy_Mock_1h log2fc_01gy_Mock_3h log2fc_1gy_Mock_3h
                 log2fc_01gy_Mock_24h log2fc_1gy_Mock_24h log2fc_01gy_Mock_72h log2fc_1gy_Mock_72h;
  set de_genes_means;
  keep gene_id log2fc_01gy_Mock_1h log2fc_1gy_Mock_1h log2fc_01gy_Mock_3h log2fc_1gy_Mock_3h
                 log2fc_01gy_Mock_24h log2fc_1gy_Mock_24h log2fc_01gy_Mock_72h log2fc_1gy_Mock_72h;
run;



data de_genes_fc1_logCPM_heatmap;
  retain gene_id log2_cpm_mock_1h log2_cpm_mock_3h log2_cpm_mock_24h log2_cpm_mock_72h
                 log2_cpm_01gy_1h log2_cpm_01gy_3h log2_cpm_01gy_24h log2_cpm_01gy_72h
                 log2_cpm_1gy_1h log2_cpm_1gy_3h log2_cpm_1gy_24h log2_cpm_1gy_72h;
  set de_genes_means_fc1;
  keep gene_id log2_cpm_mock_1h log2_cpm_mock_3h log2_cpm_mock_24h log2_cpm_mock_72h
                 log2_cpm_01gy_1h log2_cpm_01gy_3h log2_cpm_01gy_24h log2_cpm_01gy_72h
                 log2_cpm_1gy_1h log2_cpm_1gy_3h log2_cpm_1gy_24h log2_cpm_1gy_72h;
run;

data de_genes_fc1_logFC_heatmap;
  retain gene_id log2fc_01gy_Mock_1h log2fc_1gy_Mock_1h log2fc_01gy_Mock_3h log2fc_1gy_Mock_3h
                 log2fc_01gy_Mock_24h log2fc_1gy_Mock_24h log2fc_01gy_Mock_72h log2fc_1gy_Mock_72h;
  set de_genes_means_fc1;
  keep gene_id log2fc_01gy_Mock_1h log2fc_1gy_Mock_1h log2fc_01gy_Mock_3h log2fc_1gy_Mock_3h
                 log2fc_01gy_Mock_24h log2fc_1gy_Mock_24h log2fc_01gy_Mock_72h log2fc_1gy_Mock_72h;
run;

proc export data=de_genes_any_logCPM_heatmap
    outfile="!HOME/concannon/DTRA/arabidopsis_heatmap_data_logCPM_anyDEG.csv"
    dbms=csv replace;
run;

proc export data=de_genes_any_logFC_heatmap
    outfile="!HOME/concannon/DTRA/arabidopsis_heatmap_data_logFC_anyDEG.csv"
    dbms=csv replace;
run;

proc export data=de_genes_fc1_logCPM_heatmap
    outfile="!HOME/concannon/DTRA/arabidopsis_heatmap_data_logCPM_logFC1.csv"
    dbms=csv replace;
run;

proc export data=de_genes_fc1_logFC_heatmap
    outfile="!HOME/concannon/DTRA/arabidopsis_heatmap_data_logFC_logFC1.csv"
    dbms=csv replace;
run;





data de_genes_any_logCPM_hmap_cent;
  retain gene_id log2_cpm_mock_1h_c log2_cpm_mock_3h_c log2_cpm_mock_24h_c log2_cpm_mock_72h_c
                 log2_cpm_01gy_1h_c log2_cpm_01gy_3h_c log2_cpm_01gy_24h_c log2_cpm_01gy_72h_c
                 log2_cpm_1gy_1h_c log2_cpm_1gy_3h_c log2_cpm_1gy_24h_c log2_cpm_1gy_72h_c;
  set de_genes_means;
  mean_row=mean(OF log2_cpm_:);
  log2_cpm_mock_1h_c= log2_cpm_mock_1h - mean_row;
  log2_cpm_mock_3h_c= log2_cpm_mock_3h - mean_row;
  log2_cpm_mock_24h_c= log2_cpm_mock_24h - mean_row;
  log2_cpm_mock_72h_c= log2_cpm_mock_72h - mean_row;
  log2_cpm_01gy_1h_c= log2_cpm_01gy_1h - mean_row;
  log2_cpm_01gy_3h_c= log2_cpm_01gy_3h - mean_row;
  log2_cpm_01gy_24h_c= log2_cpm_01gy_24h - mean_row;
  log2_cpm_01gy_72h_c= log2_cpm_01gy_72h - mean_row;
  log2_cpm_1gy_1h_c= log2_cpm_1gy_1h - mean_row;
  log2_cpm_1gy_3h_c= log2_cpm_1gy_3h - mean_row;
  log2_cpm_1gy_24h_c= log2_cpm_1gy_24h - mean_row;
  log2_cpm_1gy_72h_c= log2_cpm_1gy_72h - mean_row;

  keep gene_id log2_cpm_mock_1h_c log2_cpm_mock_3h_c log2_cpm_mock_24h_c log2_cpm_mock_72h_c
                 log2_cpm_01gy_1h_c log2_cpm_01gy_3h_c log2_cpm_01gy_24h_c log2_cpm_01gy_72h_c
                 log2_cpm_1gy_1h_c log2_cpm_1gy_3h_c log2_cpm_1gy_24h_c log2_cpm_1gy_72h_c;
  rename log2_cpm_mock_1h_c=log2_cpm_mock_1h
         log2_cpm_mock_3h_c=log2_cpm_mock_3h
         log2_cpm_mock_24h_c=log2_cpm_mock_24h
         log2_cpm_mock_72h_c=log2_cpm_mock_72h
         log2_cpm_01gy_1h_c=log2_cpm_01gy_1h
         log2_cpm_01gy_3h_c=log2_cpm_01gy_3h
         log2_cpm_01gy_24h_c=log2_cpm_01gy_24h
         log2_cpm_01gy_72h_c=log2_cpm_01gy_72h
         log2_cpm_1gy_1h_c=log2_cpm_1gy_1h
         log2_cpm_1gy_3h_c=log2_cpm_1gy_3h
         log2_cpm_1gy_24h_c=log2_cpm_1gy_24h
         log2_cpm_1gy_72h_c=log2_cpm_1gy_72h;
run;


data de_genes_fc1_logCPM_hmap_cent;
  retain gene_id log2_cpm_mock_1h_c log2_cpm_mock_3h_c log2_cpm_mock_24h_c log2_cpm_mock_72h_c
                 log2_cpm_01gy_1h_c log2_cpm_01gy_3h_c log2_cpm_01gy_24h_c log2_cpm_01gy_72h_c
                 log2_cpm_1gy_1h_c log2_cpm_1gy_3h_c log2_cpm_1gy_24h_c log2_cpm_1gy_72h_c;
  set de_genes_means_fc1;
  mean_row=mean(OF log2_cpm_:);
  log2_cpm_mock_1h_c= log2_cpm_mock_1h - mean_row;
  log2_cpm_mock_3h_c= log2_cpm_mock_3h - mean_row;
  log2_cpm_mock_24h_c= log2_cpm_mock_24h - mean_row;
  log2_cpm_mock_72h_c= log2_cpm_mock_72h - mean_row;
  log2_cpm_01gy_1h_c= log2_cpm_01gy_1h - mean_row;
  log2_cpm_01gy_3h_c= log2_cpm_01gy_3h - mean_row;
  log2_cpm_01gy_24h_c= log2_cpm_01gy_24h - mean_row;
  log2_cpm_01gy_72h_c= log2_cpm_01gy_72h - mean_row;
  log2_cpm_1gy_1h_c= log2_cpm_1gy_1h - mean_row;
  log2_cpm_1gy_3h_c= log2_cpm_1gy_3h - mean_row;
  log2_cpm_1gy_24h_c= log2_cpm_1gy_24h - mean_row;
  log2_cpm_1gy_72h_c= log2_cpm_1gy_72h - mean_row;

  keep gene_id log2_cpm_mock_1h_c log2_cpm_mock_3h_c log2_cpm_mock_24h_c log2_cpm_mock_72h_c
                 log2_cpm_01gy_1h_c log2_cpm_01gy_3h_c log2_cpm_01gy_24h_c log2_cpm_01gy_72h_c
                 log2_cpm_1gy_1h_c log2_cpm_1gy_3h_c log2_cpm_1gy_24h_c log2_cpm_1gy_72h_c;
  rename log2_cpm_mock_1h_c=log2_cpm_mock_1h
         log2_cpm_mock_3h_c=log2_cpm_mock_3h
         log2_cpm_mock_24h_c=log2_cpm_mock_24h
         log2_cpm_mock_72h_c=log2_cpm_mock_72h
         log2_cpm_01gy_1h_c=log2_cpm_01gy_1h
         log2_cpm_01gy_3h_c=log2_cpm_01gy_3h
         log2_cpm_01gy_24h_c=log2_cpm_01gy_24h
         log2_cpm_01gy_72h_c=log2_cpm_01gy_72h
         log2_cpm_1gy_1h_c=log2_cpm_1gy_1h
         log2_cpm_1gy_3h_c=log2_cpm_1gy_3h
         log2_cpm_1gy_24h_c=log2_cpm_1gy_24h
         log2_cpm_1gy_72h_c=log2_cpm_1gy_72h;
run;




proc export data=de_genes_any_logCPM_hmap_cent
    outfile="!HOME/concannon/DTRA/arabidopsis_heatmap_data_logFC_anyDEG_centered.csv"
    dbms=csv replace;
run;

proc export data=de_genes_fc1_logCPM_hmap_cent
    outfile="!HOME/concannon/DTRA/arabidopsis_heatmap_data_logFC_logFC1_centered.csv"
    dbms=csv replace;
run;

/*   output data for scatterplots   */

proc sort data=deg_compare_any;
   by gene_id;
proc sort data=deg_compare_fc1;
   by gene_id;
proc sort data=means;
   by gene_id;
run;

data deg_any_w_means;
  merge deg_compare_any (in=in1) means (in=in2);
  by gene_id;
  if in2;
run;

data deg_fc1_w_means;
  merge deg_compare_fc1 (in=in1) means (in=in2);
  by gene_id;
  if in2;
run;


data deg_any_w_means2;
  merge deg_compare_any (in=in1) means (in=in2);
  by gene_id;
  if in1 and in2;
run;

data deg_fc1_w_means2;
  merge deg_compare_fc1 (in=in1) means (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc contents data=deg_fc1_w_means;
run;
quit;



%macro exportScatter(inData, xData, yData, upFlag, downFlag, type, outName);

data scatter_Data;
  retain gene_id &xData. &yData. DE_color ;
  set &inData.;
  length DE_color $32.;
  %if &type.=compare %then %do;
  if &upFlag. = 1 then DE_color="2_upreg";
  else if &downFlag. = 1 then DE_color="1_downreg";
  else DE_color = "0_notDE";
  %end;
  %else %do;
  if &xData. > &yData. then DE_color="X_&xdata.";
  else if &xData. < &yData. then DE_color="Y_&ydata.";
  else DE_color="0_no_bias";
  %end;
  keep gene_id &xData. &yData. DE_color ;
run;

/*proc corr data=scatter_data pearson;
   var &xData. &yData.;
run;*/

proc sort data=scatter_data;
  by DE_color;
run;

proc export data=scatter_data
     outfile="!HOME/concannon/DTRA/arabidopsis_scatterplot_&outName..csv"
     dbms=csv replace;
run;

%mend;

*proc printto print="!HOME/concannon/DTRA/arabidopsis_scatterplot_pearson.lst.txt";
*run;

%exportScatter(deg_fc1_w_means, log2_cpm_mock_1h, log2_cpm_01gy_1h, flag_DEG_up_01_1, flag_DEG_dn_01_1, compare, DEG_01Gy_Mock_1h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_3h, log2_cpm_01gy_3h, flag_DEG_up_01_3, flag_DEG_dn_01_3, compare, DEG_01Gy_Mock_3h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_24h, log2_cpm_01gy_24h, flag_DEG_up_01_24, flag_DEG_dn_01_24, compare, DEG_01Gy_Mock_24h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_72h, log2_cpm_01gy_72h, flag_DEG_up_01_72, flag_DEG_dn_01_72, compare, DEG_01Gy_Mock_72h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_1h, log2_cpm_1gy_1h, flag_DEG_up_1_1, flag_DEG_dn_1_1, compare, DEG_1Gy_Mock_1h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_3h, log2_cpm_1gy_3h, flag_DEG_up_1_3, flag_DEG_dn_1_3, compare, DEG_1Gy_Mock_3h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_24h, log2_cpm_1gy_24h, flag_DEG_up_1_24, flag_DEG_dn_1_24, compare, DEG_1Gy_Mock_24h_logCPM_logFC1);
%exportScatter(deg_fc1_w_means, log2_cpm_mock_72h, log2_cpm_1gy_72h, flag_DEG_up_1_72, flag_DEG_dn_1_72, compare, DEG_1Gy_Mock_72h_logCPM_logFC1);

%exportScatter(deg_fc1_w_means2, log2fc_01gy_Mock_1h, log2fc_1gy_Mock_1h, _null_, _null_, logFC, 01Gy_1Gy_1h_logFC_logFC1);
%exportScatter(deg_fc1_w_means2, log2fc_01gy_Mock_3h, log2fc_1gy_Mock_3h, _null_, _null_, logFC, 01Gy_1Gy_3h_logFC_logFC1);
%exportScatter(deg_fc1_w_means2, log2fc_01gy_Mock_24h, log2fc_1gy_Mock_24h, _null_, _null_, logFC, 01Gy_1Gy_24h_logFC_logFC1);
%exportScatter(deg_fc1_w_means2, log2fc_01gy_Mock_72h, log2fc_1gy_Mock_72h, _null_, _null_, logFC, 01Gy_1Gy_72h_logFC_logFC1);

%exportScatter(deg_any_w_means, log2_cpm_mock_1h, log2_cpm_01gy_1h, flag_DEG_up_01_1, flag_DEG_dn_01_1, compare, DEG_01Gy_Mock_1h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_3h, log2_cpm_01gy_3h, flag_DEG_up_01_3, flag_DEG_dn_01_3, compare, DEG_01Gy_Mock_3h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_24h, log2_cpm_01gy_24h, flag_DEG_up_01_24, flag_DEG_dn_01_24, compare, DEG_01Gy_Mock_24h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_72h, log2_cpm_01gy_72h, flag_DEG_up_01_72, flag_DEG_dn_01_72, compare, DEG_01Gy_Mock_72h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_1h, log2_cpm_1gy_1h, flag_DEG_up_1_1, flag_DEG_dn_1_1, compare, DEG_1Gy_Mock_1h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_3h, log2_cpm_1gy_3h, flag_DEG_up_1_3, flag_DEG_dn_1_3, compare, DEG_1Gy_Mock_3h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_24h, log2_cpm_1gy_24h, flag_DEG_up_1_24, flag_DEG_dn_1_24, compare, DEG_1Gy_Mock_24h_logCPM_any);
%exportScatter(deg_any_w_means, log2_cpm_mock_72h, log2_cpm_1gy_72h, flag_DEG_up_1_72, flag_DEG_dn_1_72, compare, DEG_1Gy_Mock_72h_logCPM_any);

%exportScatter(deg_any_w_means2, log2fc_01gy_Mock_1h, log2fc_1gy_Mock_1h, _null_, _null_, logFC, 01Gy_1Gy_1h_logFC_any);
%exportScatter(deg_any_w_means2, log2fc_01gy_Mock_3h, log2fc_1gy_Mock_3h, _null_, _null_, logFC, 01Gy_1Gy_3h_logFC_any);
%exportScatter(deg_any_w_means2, log2fc_01gy_Mock_24h, log2fc_1gy_Mock_24h, _null_, _null_, logFC, 01Gy_1Gy_24h_logFC_any);
%exportScatter(deg_any_w_means2, log2fc_01gy_Mock_72h, log2fc_1gy_Mock_72h, _null_, _null_, logFC, 01Gy_1Gy_72h_logFC_any);

proc printto ;
run;


/****************************************8 Export DARs and DMRs for DREME motif analysis ***************************************************************/
 


data dmr;
  set arabMAP.results_by_dmr_annot_v2;
  where num_sites_in_region >=5 ;
run;



proc sort data=dmr;
   by site_type comparison;
proc multtest inpvalues(binomial_P)=dmr fdr out=dmr_fdr noprint;
  by site_type comparison;
run;


data up_dmr_01_72 up_dmr_1_72 dn_dmr_01_72 dn_dmr_1_72
     up_dmr_01_72_cg up_dmr_1_72_cg dn_dmr_01_72_cg dn_dmr_1_72_cg
     up_dmr_01_72_chg up_dmr_1_72_chg dn_dmr_01_72_chg dn_dmr_1_72_chg
     up_dmr_01_72_chh up_dmr_1_72_chh dn_dmr_01_72_chh dn_dmr_1_72_chh ;
     set dmr_fdr;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;

     if site_type="CG" then do;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_cg;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_cg;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_cg;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_cg;

        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_cg;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_cg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_cg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_cg;

    end;

     if site_type="CHG" then do;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chg;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chg;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chg;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chg;

        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chg;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chg;

     end;

     if site_type="CHH" then do;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chh;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chh;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chh;
        if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chh;

        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chh;
     end;
     keep comparison site_type chr region_Start region_stop ;
run;


/*
 The data set WORK.UP_DMR_01_72 has 1222 observations and 5 variables.
 The data set WORK.UP_DMR_1_72 has 10511 observations and 5 variables.
 The data set WORK.DN_DMR_01_72 has 25795 observations and 5 variables.
 The data set WORK.DN_DMR_1_72 has 3651 observations and 5 variables.
 The data set WORK.UP_DMR_01_72_CG has 29 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CG has 2 observations and 5 variables.
 The data set WORK.DN_DMR_01_72_CG has 120 observations and 5 variables.
 The data set WORK.DN_DMR_1_72_CG has 20 observations and 5 variables.
 The data set WORK.UP_DMR_01_72_CHG has 39 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CHG has 7 observations and 5 variables.
 The data set WORK.DN_DMR_01_72_CHG has 323 observations and 5 variables.
 The data set WORK.DN_DMR_1_72_CHG has 3 observations and 5 variables.
 The data set WORK.UP_DMR_01_72_CHH has 1154 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CHH has 10502 observations and 5 variables.
 The data set WORK.DN_DMR_01_72_CHH has 25352 observations and 5 variables.
 The data set WORK.DN_DMR_1_72_CHH has 3628 observations and 5 variables.

*/

data dar_results_5sites;
   set arabMAP.results_by_dar_annot_v2;
   where num_sites_in_region >= 5 ;
run;

proc sort data=dar_results_5sites;
   by comparison;
proc multtest inpvalues(binomial_P)=dar_results_5sites fdr out=dar_results_5sites_fdr noprint;
  by comparison;
run;

data up_dar_01_72 up_dar_1_72 dn_dar_01_72 dn_dar_1_72;
     set dar_results_5sites_fdr;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" and num_sites_FDR05_CTL_or_TRT >=2 then output up_dar_01_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" and num_sites_FDR05_CTL_or_TRT >=2  then output dn_dar_01_72;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" and num_sites_FDR05_CTL_or_TRT >=2 then output up_dar_1_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" and num_sites_FDR05_CTL_or_TRT >=2 then output dn_dar_1_72;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;

     keep comparison site_type chr region_Start region_stop ;
run;



/*
 There were 2346 observations read from the data set WORK.DAR_RESULTS_5SITES_FDR.
The data set WORK.UP_DAR_01_72 has 1967 observations and 5 variables.
The data set WORK.UP_DAR_1_72 has 180 observations and 5 variables.
The data set WORK.DN_DAR_01_72 has 3 observations and 5 variables.
The data set WORK.DN_DAR_1_72 has 20 observations and 5 variables.

*/

proc sort data=up_dmr_01_72 nodup;  by _all_; run;
proc sort data=up_dmr_1_72 nodup;  by _all_; run;
proc sort data=dn_dmr_01_72 nodup;  by _all_; run;
proc sort data=dn_dmr_1_72 nodup;  by _all_; run;

proc sort data=up_dmr_01_72_cg nodup;  by _all_; run;
proc sort data=up_dmr_1_72_cg nodup;  by _all_; run;
proc sort data=dn_dmr_01_72_cg nodup;  by _all_; run;
proc sort data=dn_dmr_1_72_cg nodup;  by _all_; run;

proc sort data=up_dmr_01_72_chg nodup;  by _all_; run;
proc sort data=up_dmr_1_72_chg nodup;  by _all_; run;
proc sort data=dn_dmr_01_72_chg nodup;  by _all_; run;
proc sort data=dn_dmr_1_72_chg nodup;  by _all_; run;

proc sort data=up_dmr_01_72_chh nodup;  by _all_; run;
proc sort data=up_dmr_1_72_chh nodup;  by _all_; run;
proc sort data=dn_dmr_01_72_chh nodup;  by _all_; run;
proc sort data=dn_dmr_1_72_chh nodup;  by _all_; run;


proc sort data=up_dar_01_72 nodup;  by _all_; run;
proc sort data=up_dar_1_72 nodup;  by _all_; run;
proc sort data=dn_dar_01_72 nodup;  by _all_; run;
proc sort data=dn_dar_1_72 nodup;  by _all_; run;



/*
: The data set WORK.UP_DMR_01_72 has 1215 observations and 5 variables.
: The data set WORK.UP_DMR_1_72 has 10450 observations and 5 variables.
The data set WORK.DN_DMR_01_72 has 25765 observations and 5 variables.
The data set WORK.DN_DMR_1_72 has 3651 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CG has 26 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CG has 2 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CG has 105 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CG has 20 observations and 5 variables.

 The data set WORK.UP_DMR_01_72_CHG has 37 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CHG has 7 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CHG has 320 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CHG has 3 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CHH has 1152 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CHH has 10441 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CHH has 25340 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CHH has 3628 observations and 5 variables.

The data set WORK.UP_DAR_01_72 has 1965 observations and 5 variables.
The data set WORK.UP_DAR_1_72 has 180 observations and 5 variables.
The data set WORK.DN_DAR_01_72 has 3 observations and 5 variables.
The data set WORK.DN_DAR_1_72 has 20 observations and 5 variables.


*/


/* for single site * condition outputs, we can leave these alone 
   for combined DMRs, and for hyper/hypo-DMRs/DARs shared between 0.1Gy and 1Gy, we need to merge these into superregions 

   easiest way is to STACK sites, merge yte region with previous is overlapping ( lag code!! )  and track siteType and comparison
   
   for SHARED between 0.1 and 1 Gy y, count the number of occurances of super-region by unique region. If 2 then shared superregion
   */

%macro mergeDMR(dataIN);

proc sort data=&dataIN.;
   by chr region_Start region_stop;
run;

data &dataIN._make_super;
  retain superregion_num;
  set &dataIN.;
  by chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop);
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

proc sort data=&dataIN._make_super;
   by chr superregion_num region_start region_stop;
proc means data=&dataIN._make_super noprint;
   by chr superregion_num;
   var region_start region_stop;
   output out=&dataIN._merged (drop=_TYPE_ _FREQ_) min(region_start)=region_start max(region_stop)=region_stop;
run;

data &dataIN._merged2;
  set &dataIN._merged;
  length site_type $6.;
  site_type="merged";
run;

%mend;

%mergeDMR(up_dmr_01_72);
%mergeDMR(up_dmr_1_72);
%mergeDMR(dn_dmr_01_72);
%mergeDMR(dn_dmr_1_72);


%macro commonSuper(dataA, dataB, outName);

data stack_&outName.;
  set &dataA. &dataB.;
  keep site_type chr  region_start region_stop ;
run;

proc sort data=stack_&outName. nodup;
  by site_type chr  region_start region_stop;
run;


data stack_&outName._make_super;
  retain superregion_num;
  set stack_&outName.;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop);
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;


proc sort data=stack_&outName._make_super;
   by site_type chr superregion_num region_start region_stop;
proc means data=stack_&outName._make_super noprint;
   by site_type chr superregion_num;
   var region_start region_stop;
   output out=&outName._shared min(region_start)=region_start max(region_stop)=region_stop;
run;



data &outName._shared2;
  set &outName._shared;
  where _FREQ_ > 1;
  drop _TYPE_ _FREQ_;
run;

%mend;


%commonSuper(up_dmr_01_72, up_dmr_1_72, up_dmr_72);
%commonSuper(dn_dmr_01_72, dn_dmr_1_72, dn_dmr_72);
%commonSuper(up_dmr_01_72_cg, up_dmr_1_72_cg, up_dmr_72_CG);
%commonSuper(dn_dmr_01_72_cg, dn_dmr_1_72_cg, dn_dmr_72_CG);

%commonSuper(up_dmr_01_72_chg, up_dmr_1_72_chg, up_dmr_72_CHG);
%commonSuper(dn_dmr_01_72_chg, dn_dmr_1_72_chg, dn_dmr_72_CHG);
%commonSuper(up_dmr_01_72_chh, up_dmr_1_72_chh, up_dmr_72_CHH);
%commonSuper(dn_dmr_01_72_chh, dn_dmr_1_72_chh, dn_dmr_72_CHH);

%commonSuper(up_dar_01_72, up_dar_1_72, up_dar_72);
%commonSuper(dn_dar_01_72, dn_dar_1_72, dn_dar_72);




/* Export BED files for DREME analysis (which I am going to do via the
   online submission portal, because lazy and I don't want to have to
   reinstall MEME suite
*/

/* Merge hyper/hypo DMR by site type and dose */


data dmr_01_72_cg; set dn_dmr_01_72_cg up_dmr_01_72_cg; run;
data dmr_01_72_chg; set dn_dmr_01_72_chg up_dmr_01_72_chg; run;
data dmr_01_72_chh; set dn_dmr_01_72_chh up_dmr_01_72_chh; run;

data dmr_1_72_cg; set dn_dmr_1_72_cg up_dmr_1_72_cg; run;
data dmr_1_72_chg; set dn_dmr_1_72_chg up_dmr_1_72_chg; run;
data dmr_1_72_chh; set dn_dmr_1_72_chh up_dmr_1_72_chh; run;


proc sort data=dmr_01_72_cg nodup; by _all_; run;
proc sort data=dmr_01_72_chg nodup; by _all_; run;
proc sort data=dmr_01_72_chh nodup; by _all_; run;
proc sort data=dmr_1_72_cg nodup; by _all_; run;
proc sort data=dmr_1_72_chg nodup; by _all_; run;
proc sort data=dmr_1_72_chh nodup; by _all_; run;


%macro exportBED(inputData);

data format_for_bedfile;
   retain chr region_start2 region_stop2 region_id;
   set &inputData.;
   length region_id $100.;
   region_start2 = region_start - 6;
   region_stop2 = region_stop + 6;
   region_id=compress(catx("_",chr,region_start,region_stop));
   keep chr region_start2 region_stop2 region_id;
run;

proc sort data=format_for_bedfile nodup;
  by chr region_start2 region_stop2 region_id;
run;


proc export data=format_for_bedfile
     outfile="!HOME/concannon/DTRA/at_rad_motif_analysis/input_data/min5sites_2sig_&inputData..bed"
     dbms=tab replace;
     putnames=no;
run;

%mend;


%exportBED(up_dmr_01_72_cg);
%exportBED(up_dmr_1_72_cg);
%exportBED(dn_dmr_01_72_cg);
%exportBED(dn_dmr_1_72_cg); 

%exportBED(up_dmr_01_72_chg);
%exportBED(up_dmr_1_72_chg); 
%exportBED(dn_dmr_01_72_chg);
%exportBED(dn_dmr_1_72_chg); 

%exportBED(up_dmr_01_72_chh);
%exportBED(up_dmr_1_72_chh); 
%exportBED(dn_dmr_01_72_chh);
%exportBED(dn_dmr_1_72_chh); 

%exportBED(up_dar_01_72); 
%exportBED(up_dar_1_72); 
%exportBED(dn_dar_01_72);
%exportBED(dn_dar_1_72); 

%exportBED(up_dmr_01_72_merged2);
%exportBED(up_dmr_1_72_merged2);
%exportBED(dn_dmr_01_72_merged2);
%exportBED(dn_dmr_1_72_merged2);

%exportBED(up_dmr_72_shared2);
%exportBED(dn_dmr_72_shared2);
%exportBED(up_dmr_72_CG_shared2);
%exportBED(dn_dmr_72_CG_shared2);

%exportBED(up_dmr_72_CHG_shared2);
%exportBED(dn_dmr_72_CHG_shared2);
%exportBED(up_dmr_72_CHH_shared2);
%exportBED(dn_dmr_72_CHH_shared2);

%exportBED(up_dar_72_shared2);
%exportBED(dn_dar_72_shared2);



%exportBED(dmr_01_72_cg);
%exportBED(dmr_1_72_cg);
%exportBED(dmr_01_72_chg);
%exportBED(dmr_1_72_chg);
%exportBED(dmr_01_72_chh);
%exportBED(dmr_1_72_chh);

/* DAR counts:
   Hyper all, exon, intron, promoter, downstream; 0.1G, 1G
   Hypo all, exon, intron, promoter, downstream; 0.1G, 1G

 */


/****************************************8 Export DARs and DMRs for DREME motif analysis ***************************************************************/
 


data dmr;
  set arabMAP.results_by_dmr_annot_v2;
  where num_sites_in_region >=5 ;
run;



proc sort data=dmr;
   by site_type comparison;
proc multtest inpvalues(binomial_P)=dmr fdr out=dmr_fdr noprint;
  by site_type comparison;
run;


data up_dmr_01_72 up_dmr_1_72 dn_dmr_01_72 dn_dmr_1_72
     up_dmr_01_72_cg up_dmr_1_72_cg dn_dmr_01_72_cg dn_dmr_1_72_cg
     up_dmr_01_72_chg up_dmr_1_72_chg dn_dmr_01_72_chg dn_dmr_1_72_chg
     up_dmr_01_72_chh up_dmr_1_72_chh dn_dmr_01_72_chh dn_dmr_1_72_chh ;
     set dmr_fdr;
     /* FET */
     if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;

     if site_type="CG" then do;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_cg;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_cg;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_cg;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_cg;

        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_cg;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_cg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_cg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_cg;

    end;

     if site_type="CHG" then do;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chg;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chg;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chg;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chg;

        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chg;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chg;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chg;

     end;

     if site_type="CHH" then do;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chh;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chh;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chh;
        if num_sites_FDR05_diff_10perc>=5 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_chh;
        if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_chh;
     end;
     keep comparison site_type chr region_Start region_stop ;
run;


/*
 The data set WORK.UP_DMR_01_72 has 1222 observations and 5 variables.
 The data set WORK.UP_DMR_1_72 has 10511 observations and 5 variables.
 The data set WORK.DN_DMR_01_72 has 25795 observations and 5 variables.
 The data set WORK.DN_DMR_1_72 has 3651 observations and 5 variables.
 The data set WORK.UP_DMR_01_72_CG has 29 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CG has 2 observations and 5 variables.
 The data set WORK.DN_DMR_01_72_CG has 120 observations and 5 variables.
 The data set WORK.DN_DMR_1_72_CG has 20 observations and 5 variables.
 The data set WORK.UP_DMR_01_72_CHG has 39 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CHG has 7 observations and 5 variables.
 The data set WORK.DN_DMR_01_72_CHG has 323 observations and 5 variables.
 The data set WORK.DN_DMR_1_72_CHG has 3 observations and 5 variables.
 The data set WORK.UP_DMR_01_72_CHH has 1154 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CHH has 10502 observations and 5 variables.
 The data set WORK.DN_DMR_01_72_CHH has 25352 observations and 5 variables.
 The data set WORK.DN_DMR_1_72_CHH has 3628 observations and 5 variables.

*/

data dar_results_5sites;
   set arabMAP.results_by_dar_annot_v2;
   where num_sites_in_region >= 5 ;
run;

proc sort data=dar_results_5sites;
   by comparison;
proc multtest inpvalues(binomial_P)=dar_results_5sites fdr out=dar_results_5sites_fdr noprint;
  by comparison;
run;

data up_dar_01_72 up_dar_1_72 dn_dar_01_72 dn_dar_1_72;
     set dar_results_5sites_fdr;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" and num_sites_FDR05_CTL_or_TRT >=5 then output up_dar_01_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" and num_sites_FDR05_CTL_or_TRT >=5  then output dn_dar_01_72;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" and num_sites_FDR05_CTL_or_TRT >=5 then output up_dar_1_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" and num_sites_FDR05_CTL_or_TRT >=5 then output dn_dar_1_72;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;

     keep comparison site_type chr region_Start region_stop ;
run;



/*
 There were 2346 observations read from the data set WORK.DAR_RESULTS_5SITES_FDR.
The data set WORK.UP_DAR_01_72 has 1967 observations and 5 variables.
The data set WORK.UP_DAR_1_72 has 180 observations and 5 variables.
The data set WORK.DN_DAR_01_72 has 3 observations and 5 variables.
The data set WORK.DN_DAR_1_72 has 20 observations and 5 variables.

*/

proc sort data=up_dmr_01_72 nodup;  by _all_; run;
proc sort data=up_dmr_1_72 nodup;  by _all_; run;
proc sort data=dn_dmr_01_72 nodup;  by _all_; run;
proc sort data=dn_dmr_1_72 nodup;  by _all_; run;

proc sort data=up_dmr_01_72_cg nodup;  by _all_; run;
proc sort data=up_dmr_1_72_cg nodup;  by _all_; run;
proc sort data=dn_dmr_01_72_cg nodup;  by _all_; run;
proc sort data=dn_dmr_1_72_cg nodup;  by _all_; run;

proc sort data=up_dmr_01_72_chg nodup;  by _all_; run;
proc sort data=up_dmr_1_72_chg nodup;  by _all_; run;
proc sort data=dn_dmr_01_72_chg nodup;  by _all_; run;
proc sort data=dn_dmr_1_72_chg nodup;  by _all_; run;

proc sort data=up_dmr_01_72_chh nodup;  by _all_; run;
proc sort data=up_dmr_1_72_chh nodup;  by _all_; run;
proc sort data=dn_dmr_01_72_chh nodup;  by _all_; run;
proc sort data=dn_dmr_1_72_chh nodup;  by _all_; run;


proc sort data=up_dar_01_72 nodup;  by _all_; run;
proc sort data=up_dar_1_72 nodup;  by _all_; run;
proc sort data=dn_dar_01_72 nodup;  by _all_; run;
proc sort data=dn_dar_1_72 nodup;  by _all_; run;



/*
: The data set WORK.UP_DMR_01_72 has 1215 observations and 5 variables.
: The data set WORK.UP_DMR_1_72 has 10450 observations and 5 variables.
The data set WORK.DN_DMR_01_72 has 25765 observations and 5 variables.
The data set WORK.DN_DMR_1_72 has 3651 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CG has 26 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CG has 2 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CG has 105 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CG has 20 observations and 5 variables.

 The data set WORK.UP_DMR_01_72_CHG has 37 observations and 5 variables.
 The data set WORK.UP_DMR_1_72_CHG has 7 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CHG has 320 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CHG has 3 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CHH has 1152 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CHH has 10441 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CHH has 25340 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CHH has 3628 observations and 5 variables.

The data set WORK.UP_DAR_01_72 has 1965 observations and 5 variables.
The data set WORK.UP_DAR_1_72 has 180 observations and 5 variables.
The data set WORK.DN_DAR_01_72 has 3 observations and 5 variables.
The data set WORK.DN_DAR_1_72 has 20 observations and 5 variables.



The data set WORK.UP_DMR_01_72 has 1214 observations and 5 variables.
The data set WORK.UP_DMR_1_72 has 10447 observations and 5 variables.
The data set WORK.DN_DMR_01_72 has 25764 observations and 5 variables.
The data set WORK.DN_DMR_1_72 has 3621 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CG has 25 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CG has 2 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CG has 105 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CG has 20 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CHG has 37 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CHG has 7 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CHG has 319 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CHG has 3 observations and 5 variables.

The data set WORK.UP_DMR_01_72_CHH has 1152 observations and 5 variables.
The data set WORK.UP_DMR_1_72_CHH has 10438 observations and 5 variables.
The data set WORK.DN_DMR_01_72_CHH has 25340 observations and 5 variables.
The data set WORK.DN_DMR_1_72_CHH has 3598 observations and 5 variables.

The data set WORK.UP_DAR_01_72 has 1033 observations and 5 variables.
The data set WORK.UP_DAR_1_72 has 18 observations and 5 variables.
The data set WORK.DN_DAR_01_72 has 1 observations and 5 variables.
The data set WORK.DN_DAR_1_72 has 2 observations and 5 variables.

*/


/* for single site * condition outputs, we can leave these alone 
   for combined DMRs, and for hyper/hypo-DMRs/DARs shared between 0.1Gy and 1Gy, we need to merge these into superregions 

   easiest way is to STACK sites, merge yte region with previous is overlapping ( lag code!! )  and track siteType and comparison
   
   for SHARED between 0.1 and 1 Gy y, count the number of occurances of super-region by unique region. If 2 then shared superregion
   */

%macro mergeDMR(dataIN);

proc sort data=&dataIN.;
   by chr region_Start region_stop;
run;

data &dataIN._make_super;
  retain superregion_num;
  set &dataIN.;
  by chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop);
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

proc sort data=&dataIN._make_super;
   by chr superregion_num region_start region_stop;
proc means data=&dataIN._make_super noprint;
   by chr superregion_num;
   var region_start region_stop;
   output out=&dataIN._merged (drop=_TYPE_ _FREQ_) min(region_start)=region_start max(region_stop)=region_stop;
run;

data &dataIN._merged2;
  set &dataIN._merged;
  length site_type $6.;
  site_type="merged";
run;

%mend;

%mergeDMR(up_dmr_01_72);
%mergeDMR(up_dmr_1_72);
%mergeDMR(dn_dmr_01_72);
%mergeDMR(dn_dmr_1_72);


%macro commonSuper(dataA, dataB, outName);

data stack_&outName.;
  set &dataA. &dataB.;
  keep site_type chr  region_start region_stop ;
run;

proc sort data=stack_&outName. nodup;
  by site_type chr  region_start region_stop;
run;


data stack_&outName._make_super;
  retain superregion_num;
  set stack_&outName.;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop);
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;


proc sort data=stack_&outName._make_super;
   by site_type chr superregion_num region_start region_stop;
proc means data=stack_&outName._make_super noprint;
   by site_type chr superregion_num;
   var region_start region_stop;
   output out=&outName._shared min(region_start)=region_start max(region_stop)=region_stop;
run;



data &outName._shared2;
  set &outName._shared;
  where _FREQ_ > 1;
  drop _TYPE_ _FREQ_;
run;

%mend;


%commonSuper(up_dmr_01_72, up_dmr_1_72, up_dmr_72);
%commonSuper(dn_dmr_01_72, dn_dmr_1_72, dn_dmr_72);
%commonSuper(up_dmr_01_72_cg, up_dmr_1_72_cg, up_dmr_72_CG);
%commonSuper(dn_dmr_01_72_cg, dn_dmr_1_72_cg, dn_dmr_72_CG);

%commonSuper(up_dmr_01_72_chg, up_dmr_1_72_chg, up_dmr_72_CHG);
%commonSuper(dn_dmr_01_72_chg, dn_dmr_1_72_chg, dn_dmr_72_CHG);
%commonSuper(up_dmr_01_72_chh, up_dmr_1_72_chh, up_dmr_72_CHH);
%commonSuper(dn_dmr_01_72_chh, dn_dmr_1_72_chh, dn_dmr_72_CHH);

%commonSuper(up_dar_01_72, up_dar_1_72, up_dar_72);
%commonSuper(dn_dar_01_72, dn_dar_1_72, dn_dar_72);




/* Export BED files for DREME analysis (which I am going to do via the
   online submission portal, because lazy and I don't want to have to
   reinstall MEME suite
*/

/* Merge hyper/hypo DMR by site type and dose */


data dmr_01_72_cg; set dn_dmr_01_72_cg up_dmr_01_72_cg; run;
data dmr_01_72_chg; set dn_dmr_01_72_chg up_dmr_01_72_chg; run;
data dmr_01_72_chh; set dn_dmr_01_72_chh up_dmr_01_72_chh; run;

data dmr_1_72_cg; set dn_dmr_1_72_cg up_dmr_1_72_cg; run;
data dmr_1_72_chg; set dn_dmr_1_72_chg up_dmr_1_72_chg; run;
data dmr_1_72_chh; set dn_dmr_1_72_chh up_dmr_1_72_chh; run;


proc sort data=dmr_01_72_cg nodup; by _all_; run;
proc sort data=dmr_01_72_chg nodup; by _all_; run;
proc sort data=dmr_01_72_chh nodup; by _all_; run;
proc sort data=dmr_1_72_cg nodup; by _all_; run;
proc sort data=dmr_1_72_chg nodup; by _all_; run;
proc sort data=dmr_1_72_chh nodup; by _all_; run;


%macro exportBED(inputData);

data format_for_bedfile;
   retain chr region_start2 region_stop2 region_id;
   set &inputData.;
   length region_id $100.;
   region_start2 = region_start - 6;
   region_stop2 = region_stop + 6;
   region_id=compress(catx("_",chr,region_start,region_stop));
   keep chr region_start2 region_stop2 region_id;
run;

proc sort data=format_for_bedfile nodup;
  by chr region_start2 region_stop2 region_id;
run;


proc export data=format_for_bedfile
     outfile="!HOME/concannon/DTRA/at_rad_motif_analysis/input_data/min5sites_5sig_&inputData..bed"
     dbms=tab replace;
     putnames=no;
run;

%mend;


%exportBED(up_dmr_01_72_cg);
%exportBED(up_dmr_1_72_cg);
%exportBED(dn_dmr_01_72_cg);
%exportBED(dn_dmr_1_72_cg); 

%exportBED(up_dmr_01_72_chg);
%exportBED(up_dmr_1_72_chg); 
%exportBED(dn_dmr_01_72_chg);
%exportBED(dn_dmr_1_72_chg); 

%exportBED(up_dmr_01_72_chh);
%exportBED(up_dmr_1_72_chh); 
%exportBED(dn_dmr_01_72_chh);
%exportBED(dn_dmr_1_72_chh); 

%exportBED(up_dar_01_72); 
%exportBED(up_dar_1_72); 
%exportBED(dn_dar_01_72);
%exportBED(dn_dar_1_72); 

%exportBED(up_dmr_01_72_merged2);
%exportBED(up_dmr_1_72_merged2);
%exportBED(dn_dmr_01_72_merged2);
%exportBED(dn_dmr_1_72_merged2);

%exportBED(up_dmr_72_shared2);
%exportBED(dn_dmr_72_shared2);
%exportBED(up_dmr_72_CG_shared2);
%exportBED(dn_dmr_72_CG_shared2);

%exportBED(up_dmr_72_CHG_shared2);
%exportBED(dn_dmr_72_CHG_shared2);
%exportBED(up_dmr_72_CHH_shared2);
%exportBED(dn_dmr_72_CHH_shared2);

%exportBED(up_dar_72_shared2);
%exportBED(dn_dar_72_shared2);



%exportBED(dmr_01_72_cg);
%exportBED(dmr_1_72_cg);
%exportBED(dmr_01_72_chg);
%exportBED(dmr_1_72_chg);
%exportBED(dmr_01_72_chh);
%exportBED(dmr_1_72_chh);

/* DAR counts:
   Hyper all, exon, intron, promoter, downstream; 0.1G, 1G
   Hypo all, exon, intron, promoter, downstream; 0.1G, 1G

 */




/*****************************************************************************************/

data up_dar_01_72 up_dar_1_72 dn_dar_01_72 dn_dar_1_72;
     set arabMAP.results_by_dar_annot;
     length feature $32.;
     feature = scan(annotation, 1, " ");
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;

     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
run;

proc freq data=up_dar_01_72; tables feature; run;
proc freq data=dn_dar_01_72; tables feature; run;
proc freq data=up_dar_1_72; tables feature; run;
proc freq data=dn_dar_1_72; tables feature; run;

/* 
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   13000       5764    10414   3488    6762        39428
0.1Gy   Hypo    428         48      94      16      45          631
1Gy     Hyper   5500        2100    3670    1314    2413        15003
1Gy     Hypo    2173        627     1175    263     736         4975

*/

/* count by gene */

data up_dar_01_72_gn up_dar_1_72_gn dn_dar_01_72_gn dn_dar_1_72_gn;
     set arabMAP.results_by_dar_annot;
     if geneID = "" then delete;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_72_gn;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_72_gn;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72_gn;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72_gn;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72_gn;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72_gn;

     keep geneID;
run;

proc sort  data=up_dar_01_72_gn nodup; by geneID; run;
proc sort data=dn_dar_01_72_gn nodup; by geneID; run;
proc sort  data=up_dar_1_72_gn nodup; by geneID; run;
proc sort  data=dn_dar_1_72_gn nodup; by geneID; run;

data up_01_1_gn;
  merge up_dar_01_72_gn (in=in1) up_dar_1_72_gn (in=in2);
  by geneID;
  if in1 then dar_up_01=1; else dar_up_01=0;
  if in2 then dar_up_1=1; else dar_up_1=0;
run;

data dn_01_1_gn;
  merge dn_dar_01_72_gn (in=in1) dn_dar_1_72_gn (in=in2);
  by geneID;
  if in1 then dar_dn_01=1; else dar_dn_01=0;
  if in2 then dar_dn_1=1; else dar_dn_1=0;
run;


data up_dn_01_gn;
  merge up_dar_01_72_gn (in=in1) dn_dar_01_72_gn (in=in2);
  by geneID;
  if in1 then dar_up_01=1; else dar_up_01=0;
  if in2 then dar_dn_01=1; else dar_dn_01=0;
run;

data up_dn_1_gn;
  merge up_dar_1_72_gn (in=in1) dn_dar_1_72_gn (in=in2);
  by geneID;
  if in1 then dar_up_1=1; else dar_up_1=0;
  if in2 then dar_dn_1=1; else dar_dn_1=0;
run;

proc freq data=up_01_1_gn; tables dar_up_01*dar_up_1; run;
proc freq data=dn_01_1_gn; tables dar_dn_01*dar_dn_1; run;
proc freq data=up_dn_01_gn; tables dar_up_01*dar_dn_01; run;
proc freq data=up_dn_1_gn; tables dar_up_1*dar_dn_1; run;


/* 

0.1Gy hyper-DAG 175370
0.1Gy hypo-DAG  491
1Gy hyper-DAG   8898
1Gy hypo-DAG    3458

  dar_up_01     dar_up_1

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    607 |    607
           |   0.00 |   3.35 |   3.35
           |   0.00 | 100.00 |
           |   0.00 |   6.82 |
  ---------+--------+--------+
         1 |   9246 |   8291 |  17537
           |  50.96 |  45.70 |  96.65
           |  52.72 |  47.28 |
           | 100.00 |  93.18 |
  ---------+--------+--------+
  Total        9246     8898    18144
              50.96    49.04   100.00



     Table of dar_dn_01 by dar_dn_1

  dar_dn_01     dar_dn_1

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   3153 |   3153
           |   0.00 |  86.53 |  86.53
           |   0.00 | 100.00 |
           |   0.00 |  91.18 |
  ---------+--------+--------+
         1 |    186 |    305 |    491
           |   5.10 |   8.37 |  13.47
           |  37.88 |  62.12 |
           | 100.00 |   8.82 |
  ---------+--------+--------+
  Total         186     3458     3644
               5.10    94.90   100.00




      dar_up_01     dar_dn_01

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |      0 |     57 |     57
               |   0.00 |   0.32 |   0.32
               |   0.00 | 100.00 |
               |   0.00 |  11.61 |
      ---------+--------+--------+
             1 |  17103 |    434 |  17537
               |  97.21 |   2.47 |  99.68
               |  97.53 |   2.47 |
               | 100.00 |  88.39 |
      ---------+--------+--------+
      Total       17103      491    17594
                  97.21     2.79   100.00

                 The SAS System

 dar_up_1     dar_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |   1825 |   1825
          |   0.00 |  17.02 |  17.02
          |   0.00 | 100.00 |
          |   0.00 |  52.78 |
 ---------+--------+--------+
        1 |   7265 |   1633 |   8898
          |  67.75 |  15.23 |  82.98
          |  81.65 |  18.35 |
          | 100.00 |  47.22 |
 ---------+--------+--------+
 Total        7265     3458    10723
             67.75    32.25   100.00

*/



data dar_up_dn_all_gn;
  merge  up_dar_01_72_gn (in=in1) up_dar_1_72_gn (in=in2) dn_dar_01_72_gn (in=in3) dn_dar_1_72_gn (in=in4);
  by geneID;
  if in1 then dar_up_01=1; else dar_up_01=0;
  if in1 then dar_up_1=1; else dar_up_1=0;
  if in2 then dar_dn_01=1; else dar_dn_01=0;
  if in2 then dar_dn_1=1; else dar_dn_1=0;
run;

proc freq data=dar_up_dn_all_gn noprint;
  tables dar_up_01*dar_up_1*dar_dn_01*dar_dn_1 / out=dag_compare;
proc print data=dag_compare;
run;

/*
  dar_                 dar_
 up_01    dar_up_1    dn_01    dar_dn_1    COUNT

   0          0         0          0         856
   0          0         1          1         607    <- only hypoaccessible
   1          1         0          0        9246    <- only hyperaccessible
   1          1         1          1        8291    <- weirdly both???


*/



/* overlapping DARs now

need to merge region here!!

0.1Gy vs 1Gy hyper (DAR)
0.1Gy vs 1Gy hypo (DAR)
0.1Gy hyper vs hypo (DAR)
0.1Gy hyper vs hypo (DAR)

*/



%macro commonDAR(dataA, dataB, outName);

data stack_&outName.;
  set &dataA. (in=in1) &dataB. (in=in2);
  length comp $20.;
  if in1 then comp="&dataA.";
  if in2 then comp="&dataB.";
  keep comp site_type chr  region_start region_stop ;
run;

proc sort data=stack_&outName. nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_&outName._make_super;
  retain superregion_num;
  set stack_&outName.;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_&outName._super1;
   set stack_&outName._make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_&outName._super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_&outName._super1 out=stack_&outName._super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;

data stack_&outName._super_sbys2;
  set stack_&outName._super_sbys;
  if &dataA.=. then &dataA.=0;
  if &dataB.=. then &dataB.=0;
run;

proc freq data=stack_&outName._super_sbys2;
  tables &dataA.*&dataB. ;
run;
%mend;

%commonDAR(up_dar_01_72, up_dar_1_72, up_dar_72);
%commonDAR(dn_dar_01_72, dn_dar_1_72, dn_dar_72);
%commonDAR(up_dar_01_72, dn_dar_01_72, updn_01_dar_72);
%commonDAR(up_dar_1_72, dn_dar_1_72, updn_1_dar_72);


/*
 up_dar_01_72     up_dar_1_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |   3611 |   3611
          |   0.00 |   8.40 |   8.40
          |   0.00 | 100.00 |
          |   0.00 |  24.16 |
 ---------+--------+--------+
        1 |  28062 |  11335 |  39397
          |  65.25 |  26.36 |  91.60
          |  71.23 |  28.77 |
          | 100.00 |  75.84 |
 ---------+--------+--------+
 Total       28062    14946    43008
             65.25    34.75   100.00


 dn_dar_01_72     dn_dar_1_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |   4899 |   4899
          |   0.00 |  88.59 |  88.59
          |   0.00 | 100.00 |
          |   0.00 |  98.49 |
 ---------+--------+--------+
        1 |    556 |     75 |    631
          |  10.05 |   1.36 |  11.41
          |  88.11 |  11.89 |
          | 100.00 |   1.51 |
 ---------+--------+--------+
 Total         556     4974     5530
             10.05    89.95   100.00


  up_dar_01_72     dn_dar_01_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    631 |    631
           |   0.00 |   1.58 |   1.58
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |  39428 |      0 |  39428
           |  98.42 |   0.00 |  98.42
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total       39428      631    40059
              98.42     1.58   100.00


  up_dar_1_72     dn_dar_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   4974 |   4974
           |   0.00 |  24.90 |  24.90
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |  15003 |      0 |  15003
           |  75.10 |   0.00 |  75.10
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total       15003     4974    19977
              75.10    24.90   100.00

*/

/* 4 way venn of regions */





data stack_all_Dar;
  set up_dar_01_72 (in=in1) dn_dar_01_72 (in=in2) up_dar_1_72 (in=in3) dn_dar_1_72 (in=in4);
  length comp $20.;
  if in1 then comp="up_dar_01_72";
  if in2 then comp="dn_dar_01_72";
  if in3 then comp="up_dar_1_72";
  if in4 then comp="dn_dar_1_72";
  keep comp site_type chr  region_start region_stop ;
run;

proc sort data=stack_all_Dar nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_all_Dar_make_super;
  retain superregion_num;
  set stack_all_Dar;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_all_Dar_super1;
   set stack_all_Dar_make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_all_Dar_super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_all_Dar_super1 out=stack_all_Dar_super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;


data stack_all_Dar_super_sbys2;
  set stack_all_Dar_super_sbys;
  if up_dar_01_72=. then up_dar_01_72=0;
  if dn_dar_01_72=. then dn_dar_01_72=0;
  if up_dar_1_72=. then up_dar_1_72=0;
  if dn_dar_1_72=. then dn_dar_1_72=0;
run;

proc freq data=stack_all_Dar_super_sbys2 noprint;
  tables up_dar_01_72*up_dar_1_72*dn_dar_01_72*dn_dar_1_72 / out=dar_compare;
run;
proc print data=dar_compare;
run;


/*
 up_dar_    up_dar_    dn_dar_    dn_dar_
  01_72       1_72      01_72       1_72     COUNT    PERCENT

    0          0          0          1        4470     9.2916
    0          0          1          0         541     1.1246
    0          0          1          1          70     0.1455
    0          1          0          0        3615     7.5143
    0          1          1          0          15     0.0312
    1          0          0          0       27683    57.5434
    1          0          0          1         393     0.8169
    1          0          1          1           1     0.0021
    1          1          0          0       11277    23.4410
    1          1          0          1          39     0.0811
    1          1          1          0           3     0.0062
    1          1          1          1           1     0.0021

Conclusion on DARs:
when you irradiate plants with low-dose gamma radiation, accessibilty decreases in relatively few regions consistently 
Instead, there is a substantially more regions with increased chromatin accessibility

In other words,the response is to open up a lot of chromatin and express a bunch of genes (or ready/prime the expression
of a bunch of genes */


/* PRep and export data for the following LINE PLOTS:

(1) Average GC accessibility in hyper/hypo DARs
(2) TSS accessibility plots

*/





* GC methylation data;
data gc_data;
  set arabMAP.methylation_data_gc;
  where flag_normalized=1;
  keep chr stop_pos treatment units rep total_C_norm methyl_C_norm perc_methyl_norm;
run;

proc sort data=gc_data;
  by chr stop_pos treatment units rep ;
proc means data=gc_data noprint;
  by chr stop_pos treatment units  ;
  var total_C_norm methyl_C_norm perc_methyl_norm;
  output out=gc_data2 sum(total_C_norm)=total_C sum(methyl_C_norm)=methyl_C mean(perc_methyl_norm)=perc_methyl;
run;

data gc_data3;
  set gc_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=gc_data3 out=gc_sbys10;
  where total_C >= 10;
  by chr stop_pos;
  id treatment units;
  var perc_methyl2;
run;

data gc_sbys10_2;
  set gc_sbys10;
  if _01Gy0U ne  . and _01Gy100U ne . then _01Gy_100U_0U=_01Gy100U - _01Gy0U;
  else _01Gy_100U_0U=.;

  if _1Gy0U ne . and _1Gy100U ne . then _1Gy_100U_0U=_1Gy100U - _1Gy0U;
  else _1Gy_100U_0U=.;

  if _0Gy0U ne . and _0Gy100U ne . then _0Gy_100U_0U=_0Gy100U - _0Gy0U;
  else _0Gy_100U_0U=.;

  if _01Gy_100U_0U ne . and _0Gy_100U_0U ne . then _01Gy_100U_0U_common=_01Gy_100U_0U; else _01Gy_100U_0U_common=.;
  if _1Gy_100U_0U ne . and _0Gy_100U_0U ne . then _1Gy_100U_0U_common=_1Gy_100U_0U; else _1Gy_100U_0U_common=.;
  if (_1Gy_100U_0U ne . or _01Gy_100U_0U ne .) and _0Gy_100U_0U ne . then _0Gy_100U_0U_common=_0Gy_100U_0U; else _0Gy_100U_0U_common=.;

  if _01Gy_100U_0U_common ne . and _1Gy_100U_0U_common ne .  then do;
    _01Gy_100U_0U_common_all = _01Gy_100U_0U_common;
    _1Gy_100U_0U_common_all = _1Gy_100U_0U_common;
    _0Gy_100U_0U_common_all = _0Gy_100U_0U_common;
    end;
  else do;
   _01Gy_100U_0U_common_all = .;
   _1Gy_100U_0U_common_all = .;
   _0Gy_100U_0U_common_all = . ;
    end;
  rename stop_pos=pos;
run;


/* Get DARs */


data up_dar_01_1kb up_dar_1_1kb dn_dar_01_1kb dn_dar_1_1kb;
     set arabMAP.results_by_dar_annot;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";


     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;

     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;


     keep comparison chr dar_Center plot_start plot_stop ;
run;

proc sort data=up_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;





%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=gc_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) gc_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dar_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_100U_0U_common _1Gy_100U_0U_common  _0Gy_100U_0U_common
       _01Gy_100U_0U_common_all _1Gy_100U_0U_common_all  _0Gy_100U_0U_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean( _0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean( _0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy_100U_0U_common mean_1Gy_100U_0U_common  mean_0Gy_100U_0U_common
       mean_01Gy_100U_0U_common_all mean_1Gy_100U_0U_common_all  mean_0Gy_100U_0U_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(mean_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean(mean_0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(mean_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(mean_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean(mean_0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dar_01_1kb);
%mergeMETH(dn_dar_01_1kb);
%mergeMETH(up_dar_1_1kb);
%mergeMETH(dn_dar_1_1kb);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/GC_accessibility_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hyper_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hyper_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hyper_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hyper_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hyper_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hyper_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hyper_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hyper_DAR_1kb_1Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hypo_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hypo_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hypo_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hypo_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hypo_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hypo_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hypo_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hypo_DAR_1kb_1Gy_0Gy_common_binned);

/* As above, but now for TSSs */


data site2promoter;
  set arabMAP.results_by_dac_annot;
  length gene_id $20.;
  if abs(distance_to_tss) > 999 then delete;
  if count(nearest_promoterID, "-T1") > 0 then gene_ID=compress(upcase(tranwrd(Nearest_PromoterID,"-T1","")));
  else gene_ID=compress(upcase(scan(Nearest_PromoterID,1,".")));
  keep gene_ID chr start_pos stop_pos strand distance_to_tss nearest_promoterID;
  rename stop_pos=pos nearest_promoterID=transcript_id;
run;

proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/Arabidopsis_thaliana.TAIR10.37.gtf"
  out=gtf dbms=tab replace;
  guessingrows=max;
  getnames=no;
run;

data gene;
  set gtf;
  where VAR3 = "gene";
  length gene_id $15.;
  gene_id=compress(tranwrd(tranwrd(scan(VAR9,2," "), ";", ""), '"', ''));
  keep gene_id VAR7; 
  rename VAR7=strand;
run;

proc sort data=gene nodup;
  by gene_id;
proc sort data=site2promoter nodup;
    by gene_id;
run;

data site2promoter2 no_strand no_site;
  merge site2promoter (in=in1) gene (in=in2);
  by gene_id;
  if in1 and in2 then output site2promoter2;
  else if in1 then output no_strand;
  else output no_site;
run;

/* distance to TSS is relative to strand of transcript, so I don't need to flip anything!!! */

proc sort data=site2promoter2 nodup;
  by chr pos;
proc sort data=gc_sbys10_2;
  by chr pos;
run;

data site2promoter2_w_meth;
  merge site2promoter2 (in=in1) gc_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data site2promoter2_w_meth2;
  set site2promoter2_w_meth;
  grouped_pos=int(distance_to_TSS/10) * 10;
  if distance_to_TSS < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=site2promoter2_w_meth2;
  by distance_to_TSS grouped_pos2 ;
run;


proc means data=site2promoter2_w_meth2 noprint;
  by distance_to_TSS grouped_pos2  ;
  var  _01Gy_100U_0U_common _1Gy_100U_0U_common  _0Gy_100U_0U_common
       _01Gy_100U_0U_common_all _1Gy_100U_0U_common_all  _0Gy_100U_0U_common_all ;
  output out=mean_diff_tss
  mean(_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean( _0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean( _0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;

proc sort data=mean_diff_tss ;
  by grouped_pos2;
run;


proc means data=mean_diff_tss noprint;
  by grouped_pos2  ;
  var  mean_01Gy_100U_0U_common mean_1Gy_100U_0U_common  mean_0Gy_100U_0U_common
       mean_01Gy_100U_0U_common_all mean_1Gy_100U_0U_common_all  mean_0Gy_100U_0U_common_all ;
  output out=mean_diff_tss_2
  mean(mean_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(mean_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean(mean_0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(mean_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(mean_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean(mean_0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;


data mean_diff_tss_1;
  set mean_diff_tss;
  drop _TYPE_ _FREQ_;
  keep distance_To_tss mean_: ;
  rename distance_to_tss=pos;
run;

data mean_diff_tss_3;
  set mean_diff_tss_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;


/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, var3, outName);

data export;
  retain pos &var1. &var2. &var3.;
  set &inData.;
  keep pos &var1. &var2. &var3.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/GC_accessibility_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_tss_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, mean_1Gy_100U_0U_common, GC_acc_TSS_1kb);
%exportLine( mean_diff_tss_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, mean_1Gy_100U_0U_common, GC_acc_TSS_1kb_binned);

%exportLine( mean_diff_tss_1, mean_0Gy_100U_0U_common_all, mean_01Gy_100U_0U_common_all, mean_1Gy_100U_0U_common_all, GC_acc_TSS_1kb_common);
%exportLine( mean_diff_tss_3, mean_0Gy_100U_0U_common_all, mean_01Gy_100U_0U_common_all, mean_1Gy_100U_0U_common_all, GC_acc_TSS_1kb_common_binned);



/*************************************************************************************/


/* Methylation plots:
   same as GC plots (se we can reuse the same code!)
   but also need to add 
   methylation in DMRs by site type (i.e. CG for CG DMRs, CHG for CHG DMRs, CHH for CHH DMRs) not just DARs (but code is good)

   Make a big-ass macro to do everything for ONE site type at a time, then kick off CG, CHG, CHH

*/


%macro methDataGen(siteType);


/* DMR counts:
   Hyper all, exon, intron, promoter, downstream; 0.1G, 1G
   Hypo all, exon, intron, promoter, downstream; 0.1G, 1G

 */

data up_dmr_01_72 up_dmr_1_72 dn_dmr_01_72 dn_dmr_1_72;
     set arabMAP.results_by_dmr_annot;
     length feature $32.;
     where site_type="&siteType.";
     feature = scan(annotation, 1, " ");
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;

     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
run;

proc sort data=up_dmr_01_72 nodup; by _all_; run;
proc sort data=dn_dmr_01_72 nodup; by _all_; run;
proc sort data=up_dmr_1_72 nodup; by _all_; run;
proc sort data=dn_dmr_1_72 nodup; by _all_; run;


proc freq data=up_dmr_01_72; tables feature; run;
proc freq data=dn_dmr_01_72; tables feature; run;
proc freq data=up_dmr_1_72; tables feature; run;
proc freq data=dn_dmr_1_72; tables feature; run;

/* 

CG sites:
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   5           9       10      0       3           27
0.1Gy   Hypo    21          20      11      6       54          112
1Gy     Hyper   0           1       1       0       0           2
1Gy     Hypo    2           4       1       1       12          20



CHG sites:
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   9           12      13      0       3           37
0.1Gy   Hypo    258         29      4       9       25          325
1Gy     Hyper   4           2       1       0       0           7
1Gy     Hypo    2           0       0       0       1           3


CHH sites:
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   737         119     126     31      139         1152
0.1Gy   Hypo    16996       2510    531     603     4705        25345
1Gy     Hyper   7203        978     271     229     1762        10443
1Gy     Hypo    2372        386     96      93      617         3564


*/

/* count by gene */

data up_dmr_01_72_gn up_dmr_1_72_gn dn_dmr_01_72_gn dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="&siteType.";
     if geneID = "" then delete;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_gn;
     keep geneID;
run;

proc sort  data=up_dmr_01_72_gn nodup; by geneID; run;
proc sort data=dn_dmr_01_72_gn nodup; by geneID; run;
proc sort  data=up_dmr_1_72_gn nodup; by geneID; run;
proc sort  data=dn_dmr_1_72_gn nodup; by geneID; run;

data up_01_1_gn;
  merge up_dmr_01_72_gn (in=in1) up_dmr_1_72_gn (in=in2);
  by geneID;
  if in1 then dmr_up_01=1; else dmr_up_01=0;
  if in2 then dmr_up_1=1; else dmr_up_1=0;
run;

data dn_01_1_gn;
  merge dn_dmr_01_72_gn (in=in1) dn_dmr_1_72_gn (in=in2);
  by geneID;
  if in1 then dmr_dn_01=1; else dmr_dn_01=0;
  if in2 then dmr_dn_1=1; else dmr_dn_1=0;
run;


data up_dn_01_gn;
  merge up_dmr_01_72_gn (in=in1) dn_dmr_01_72_gn (in=in2);
  by geneID;
  if in1 then dmr_up_01=1; else dmr_up_01=0;
  if in2 then dmr_dn_01=1; else dmr_dn_01=0;
run;

data up_dn_1_gn;
  merge up_dmr_1_72_gn (in=in1) dn_dmr_1_72_gn (in=in2);
  by geneID;
  if in1 then dmr_up_1=1; else dmr_up_1=0;
  if in2 then dmr_dn_1=1; else dmr_dn_1=0;
run;

proc freq data=up_01_1_gn; tables dmr_up_01*dmr_up_1; run;
proc freq data=dn_01_1_gn; tables dmr_dn_01*dmr_dn_1; run;
proc freq data=up_dn_01_gn; tables dmr_up_01*dmr_dn_01; run;
proc freq data=up_dn_1_gn; tables dmr_up_1*dmr_dn_1; run;


/* 
CG:
0.1Gy hyper-DMG 25
0.1Gy hypo-DMG  92
1Gy hyper-DMG   2
1Gy hypo-DMG    20


       Table of dmr_up_01 by dmr_up_1

    dmr_up_01     dmr_up_1

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |      0 |      1 |      1
             |   0.00 |   3.85 |   3.85
             |   0.00 | 100.00 |
             |   0.00 |  50.00 |
    ---------+--------+--------+
           1 |     24 |      1 |     25
             |  92.31 |   3.85 |  96.15
             |  96.00 |   4.00 |
             | 100.00 |  50.00 |
    ---------+--------+--------+
    Total          24        2       26
                92.31     7.69   100.00


    Table of dmr_dn_01 by dmr_dn_1

 dmr_dn_01     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |      6 |      6
          |   0.00 |   6.12 |   6.12
          |   0.00 | 100.00 |
          |   0.00 |  30.00 |
 ---------+--------+--------+
        1 |     78 |     14 |     92
          |  79.59 |  14.29 |  93.88
          |  84.78 |  15.22 |
          | 100.00 |  70.00 |
 ---------+--------+--------+
 Total          78       20       98
             79.59    20.41   100.00

            The SAS System


   dmr_up_01     dmr_dn_01

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |     91 |     91
            |   0.00 |  78.45 |  78.45
            |   0.00 | 100.00 |
            |   0.00 |  98.91 |
   ---------+--------+--------+
          1 |     24 |      1 |     25
            |  20.69 |   0.86 |  21.55
            |  96.00 |   4.00 |
            | 100.00 |   1.09 |
   ---------+--------+--------+
   Total          24       92      116
               20.69    79.31   100.00



 dmr_up_1     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |     20 |     20
          |   0.00 |  90.91 |  90.91
          |   0.00 | 100.00 |
          |   0.00 | 100.00 |
 ---------+--------+--------+
        1 |      2 |      0 |      2
          |   9.09 |   0.00 |   9.09
          | 100.00 |   0.00 |
          | 100.00 |   0.00 |
 ---------+--------+--------+
 Total           2       20       22
              9.09    90.91   100.00




CHG:
0.1Gy hyper-DMG 
0.1Gy hypo-DMG  
1Gy hyper-DMG   
1Gy hypo-DMG    


     Table of dmr_up_01 by dmr_up_1

  dmr_up_01     dmr_up_1

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      6 |      6
           |   0.00 |  14.29 |  14.29
           |   0.00 | 100.00 |
           |   0.00 |  85.71 |
  ---------+--------+--------+
         1 |     35 |      1 |     36
           |  83.33 |   2.38 |  85.71
           |  97.22 |   2.78 |
           | 100.00 |  14.29 |
  ---------+--------+--------+
  Total          35        7       42
              83.33    16.67   100.00

             The SAS System


   Table of dmr_dn_01 by dmr_dn_1

dmr_dn_01     dmr_dn_1

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |      0 |      1 |      1
         |   0.00 |   0.37 |   0.37
         |   0.00 | 100.00 |
         |   0.00 |  33.33 |
---------+--------+--------+
       1 |    269 |      2 |    271
         |  98.90 |   0.74 |  99.63
         |  99.26 |   0.74 |
         | 100.00 |  66.67 |
---------+--------+--------+
Total         269        3      272
            98.90     1.10   100.00

   Table of dmr_up_01 by dmr_dn_01

 dmr_up_01     dmr_dn_01

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    266 |    266
          |   0.00 |  88.08 |  88.08
          |   0.00 | 100.00 |
          |   0.00 |  98.15 |
 ---------+--------+--------+
        1 |     31 |      5 |     36
          |  10.26 |   1.66 |  11.92
          |  86.11 |  13.89 |
          | 100.00 |   1.85 |
 ---------+--------+--------+
 Total          31      271      302
             10.26    89.74   100.00

     Table of dmr_up_1 by dmr_dn_1

  dmr_up_1     dmr_dn_1

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      3 |      3
           |   0.00 |  30.00 |  30.00
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |      7 |      0 |      7
           |  70.00 |   0.00 |  70.00
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total           7        3       10
              70.00    30.00   100.00

             The SAS System



CHH:
0.1Gy hyper-DMG 
0.1Gy hypo-DMG  
1Gy hyper-DMG   
1Gy hypo-DMG    

      Table of dmr_up_01 by dmr_up_1

   dmr_up_01     dmr_up_1

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |   3134 |   3134
            |   0.00 |  79.16 |  79.16
            |   0.00 | 100.00 |
            |   0.00 |  83.71 |
   ---------+--------+--------+
          1 |    215 |    610 |    825
            |   5.43 |  15.41 |  20.84
            |  26.06 |  73.94 |
            | 100.00 |  16.29 |
   ---------+--------+--------+
   Total         215     3744     3959
                5.43    94.57   100.00

              The SAS System

            The FREQ Procedure

    Table of dmr_dn_01 by dmr_dn_1

 dmr_dn_01     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    119 |    119
          |   0.00 |   1.84 |   1.84
          |   0.00 | 100.00 |
          |   0.00 |   5.62 |
 ---------+--------+--------+
        1 |   4366 |   1997 |   6363
          |  67.36 |  30.81 |  98.16
          |  68.62 |  31.38 |
          | 100.00 |  94.38 |
 ---------+--------+--------+
 Total        4366     2116     6482
             67.36    32.64   100.00

    Table of dmr_up_01 by dmr_dn_01

  dmr_up_01     dmr_dn_01

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   5735 |   5735
           |   0.00 |  87.42 |  87.42
           |   0.00 | 100.00 |
           |   0.00 |  90.13 |
  ---------+--------+--------+
         1 |    197 |    628 |    825
           |   3.00 |   9.57 |  12.58
           |  23.88 |  76.12 |
           | 100.00 |   9.87 |
  ---------+--------+--------+
  Total         197     6363     6560
               3.00    97.00   100.00

    Table of dmr_up_1 by dmr_dn_1

 dmr_up_1     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    794 |    794
          |   0.00 |  17.50 |  17.50
          |   0.00 | 100.00 |
          |   0.00 |  37.52 |
 ---------+--------+--------+
        1 |   2422 |   1322 |   3744
          |  53.37 |  29.13 |  82.50
          |  64.69 |  35.31 |
          | 100.00 |  62.48 |
 ---------+--------+--------+
 Total        2422     2116     4538
             53.37    46.63   100.00

*/



data dmr_up_dn_all_gn;
  merge  up_dmr_01_72_gn (in=in1) up_dmr_1_72_gn (in=in2) dn_dmr_01_72_gn (in=in3) dn_dmr_1_72_gn (in=in4);
  by geneID;
  if in1 then dmr_up_01=1; else dmr_up_01=0;
  if in1 then dmr_up_1=1; else dmr_up_1=0;
  if in2 then dmr_dn_01=1; else dmr_dn_01=0;
  if in2 then dmr_dn_1=1; else dmr_dn_1=0;
run;

proc freq data=dmr_up_dn_all_gn noprint;
  tables dmr_up_01*dmr_up_1*dmr_dn_01*dmr_dn_1 / out=dag_compare;
proc print data=dag_compare;
run;

/*
CG:

    dmr_                 dmr_
   up_01    dmr_up_1    dn_01    dmr_dn_1    COUNT

     0          0         0          0         97
     0          0         1          1          1
     1          1         0          0         24
     1          1         1          1          1





CHG:

   dmr_                 dmr_
  up_01    dmr_up_1    dn_01    dmr_dn_1    COUNT

    0          0         0          0        267
    0          0         1          1          6
    1          1         0          0         35
    1          1         1          1          1



CHH:

     dmr_                 dmr_
    up_01    dmr_up_1    dn_01    dmr_dn_1    COUNT

      0          0         0          0        3268
      0          0         1          1        3134
      1          1         0          0         215
      1          1         1          1         610



*/



/* overlapping DMRs now

need to merge region here!!

0.1Gy vs 1Gy hyper (DMR)
0.1Gy vs 1Gy hypo (DMR)
0.1Gy hyper vs hypo (DMR)
0.1Gy hyper vs hypo (DMR)

*/



%macro commonDMR(dataA, dataB, outName);

data stack_&outName.;
  set &dataA. (in=in1) &dataB. (in=in2);
  length comp $20.;
  if in1 then comp="&dataA.";
  if in2 then comp="&dataB.";
  keep comp site_type chr  region_start region_stop ;
run;

proc sort data=stack_&outName. nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_&outName._make_super;
  retain superregion_num;
  set stack_&outName.;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_&outName._super1;
   set stack_&outName._make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_&outName._super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_&outName._super1 out=stack_&outName._super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;

data stack_&outName._super_sbys2;
  set stack_&outName._super_sbys;
  if &dataA.=. then &dataA.=0;
  if &dataB.=. then &dataB.=0;
run;

proc freq data=stack_&outName._super_sbys2;
  tables &dataA.*&dataB. ;
run;
%mend;

%commonDMR(up_dmr_01_72, up_dmr_1_72, up_dmr_72);
%commonDMR(dn_dmr_01_72, dn_dmr_1_72, dn_dmr_72);
%commonDMR(up_dmr_01_72, dn_dmr_01_72, updn_01_dmr_72);
%commonDMR(up_dmr_1_72, dn_dmr_1_72, updn_1_dmr_72);


/*


CG:
   Table of up_dmr_01_72 by up_dmr_1_72

   up_dmr_01_72     up_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |      1 |      1
            |   0.00 |   3.57 |   3.57
            |   0.00 | 100.00 |
            |   0.00 |  50.00 |
   ---------+--------+--------+
          1 |     26 |      1 |     27
            |  92.86 |   3.57 |  96.43
            |  96.30 |   3.70 |
            | 100.00 |  50.00 |
   ---------+--------+--------+
   Total          26        2       28
               92.86     7.14   100.00

              The SAS System


  dn_dmr_01_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      6 |      6
           |   0.00 |   5.17 |   5.17
           |   0.00 | 100.00 |
           |   0.00 |  30.00 |
  ---------+--------+--------+
         1 |     96 |     14 |    110
           |  82.76 |  12.07 |  94.83
           |  87.27 |  12.73 |
           | 100.00 |  70.00 |
  ---------+--------+--------+
  Total          96       20      116
              82.76    17.24   100.00


Table of up_dmr_01_72 by dn_dmr_01_72

 up_dmr_01_72     dn_dmr_01_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    112 |    112
          |   0.00 |  80.58 |  80.58
          |   0.00 | 100.00 |
          |   0.00 | 100.00 |
 ---------+--------+--------+
        1 |     27 |      0 |     27
          |  19.42 |   0.00 |  19.42
          | 100.00 |   0.00 |
          | 100.00 |   0.00 |
 ---------+--------+--------+
 Total          27      112      139
             19.42    80.58   100.00


  Table of up_dmr_1_72 by dn_dmr_1_72

  up_dmr_1_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |     20 |     20
           |   0.00 |  90.91 |  90.91
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |      2 |      0 |      2
           |   9.09 |   0.00 |   9.09
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total           2       20       22
               9.09    90.91   100.00



CHG:
   Table of up_dmr_01_72 by up_dmr_1_72

   up_dmr_01_72     up_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |      6 |      6
            |   0.00 |  13.95 |  13.95
            |   0.00 | 100.00 |
            |   0.00 |  85.71 |
   ---------+--------+--------+
          1 |     36 |      1 |     37
            |  83.72 |   2.33 |  86.05
            |  97.30 |   2.70 |
            | 100.00 |  14.29 |
   ---------+--------+--------+
   Total          36        7       43
               83.72    16.28   100.00

 Table of dn_dmr_01_72 by dn_dmr_1_72

 dn_dmr_01_72     dn_dmr_1_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |      2 |      2
          |   0.00 |   0.61 |   0.61
          |   0.00 | 100.00 |
          |   0.00 |  66.67 |
 ---------+--------+--------+
        1 |    324 |      1 |    325
          |  99.08 |   0.31 |  99.39
          |  99.69 |   0.31 |
          | 100.00 |  33.33 |
 ---------+--------+--------+
 Total         324        3      327
             99.08     0.92   100.00

            The SAS System

 Table of up_dmr_01_72 by dn_dmr_01_72

  up_dmr_01_72     dn_dmr_01_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    325 |    325
           |   0.00 |  89.78 |  89.78
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |     37 |      0 |     37
           |  10.22 |   0.00 |  10.22
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total          37      325      362
              10.22    89.78   100.00

  Table of up_dmr_1_72 by dn_dmr_1_72

  up_dmr_1_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      3 |      3
           |   0.00 |  30.00 |  30.00
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |      7 |      0 |      7
           |  70.00 |   0.00 |  70.00
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total           7        3       10
              70.00    30.00   100.00

CHH:



   Table of up_dmr_01_72 by up_dmr_1_72

   up_dmr_01_72     up_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |  10118 |  10118
            |   0.00 |  89.80 |  89.80
            |   0.00 | 100.00 |
            |   0.00 |  96.96 |
   ---------+--------+--------+
          1 |    832 |    317 |   1149
            |   7.38 |   2.81 |  10.20
            |  72.41 |  27.59 |
            | 100.00 |   3.04 |
   ---------+--------+--------+
   Total         832    10435    11267
                7.38    92.62   100.00




  dn_dmr_01_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    876 |    876
           |   0.00 |   3.37 |   3.37
           |   0.00 | 100.00 |
           |   0.00 |  24.85 |
  ---------+--------+--------+
         1 |  22439 |   2649 |  25088
           |  86.42 |  10.20 |  96.63
           |  89.44 |  10.56 |
           | 100.00 |  75.15 |
  ---------+--------+--------+
  Total       22439     3525    25964
              86.42    13.58   100.00


Table of up_dmr_01_72 by dn_dmr_01_72

 up_dmr_01_72     dn_dmr_01_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |  25345 |  25345
          |   0.00 |  95.65 |  95.65
          |   0.00 | 100.00 |
          |   0.00 | 100.00 |
 ---------+--------+--------+
        1 |   1152 |      0 |   1152
          |   4.35 |   0.00 |   4.35
          | 100.00 |   0.00 |
          | 100.00 |   0.00 |
 ---------+--------+--------+
 Total        1152    25345    26497
              4.35    95.65   100.00

 Table of up_dmr_1_72 by dn_dmr_1_72

    up_dmr_1_72     dn_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |   3564 |   3564
            |   0.00 |  25.44 |  25.44
            |   0.00 | 100.00 |
            |   0.00 | 100.00 |
   ---------+--------+--------+
          1 |  10443 |      0 |  10443
            |  74.56 |   0.00 |  74.56
            | 100.00 |   0.00 |
            | 100.00 |   0.00 |
   ---------+--------+--------+
   Total       10443     3564    14007
               74.56    25.44   100.00

*/

/* 4 way venn of regions */


data stack_all_dmr;
  set up_dmr_01_72 (in=in1) dn_dmr_01_72 (in=in2) up_dmr_1_72 (in=in3) dn_dmr_1_72 (in=in4);
  length comp $20.;
  if in1 then comp="up_dmr_01_72";
  if in2 then comp="dn_dmr_01_72";
  if in3 then comp="up_dmr_1_72";
  if in4 then comp="dn_dmr_1_72";
  keep comp site_type chr  region_start region_stop ;
run;

proc sort data=stack_all_dmr nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_all_dmr_make_super;
  retain superregion_num;
  set stack_all_dmr;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_all_dmr_super1;
   set stack_all_dmr_make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_all_dmr_super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_all_dmr_super1 out=stack_all_dmr_super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;


data stack_all_dmr_super_sbys2;
  set stack_all_dmr_super_sbys;
  if up_dmr_01_72=. then up_dmr_01_72=0;
  if dn_dmr_01_72=. then dn_dmr_01_72=0;
  if up_dmr_1_72=. then up_dmr_1_72=0;
  if dn_dmr_1_72=. then dn_dmr_1_72=0;
run;

proc freq data=stack_all_dmr_super_sbys2 noprint;
  tables up_dmr_01_72*up_dmr_1_72*dn_dmr_01_72*dn_dmr_1_72 / out=dmr_compare;
run;
proc print data=dmr_compare;
run;


/*
CG:

up_dmr_    up_dmr_    dn_dmr_    dn_dmr_
 01_72       1_72      01_72       1_72     COUNT

   0          0          0          1          6
   0          0          1          0         96
   0          0          1          1         14 <- hypo
   0          1          0          0          1
   1          0          0          0         26
   1          1          0          0          1 <- hyper

                    The SAS System

CHG:


 up_dmr_    up_dmr_    dn_dmr_    dn_dmr_
  01_72       1_72      01_72       1_72     COUNT

    0          0          0          1          2
    0          0          1          0        324
    0          0          1          1          1 <- hypo
    0          1          0          0          6
    1          0          0          0         36
    1          1          0          0          1 <- hyper


CHH:


 up_dmr_    up_dmr_    dn_dmr_    dn_dmr_
  01_72       1_72      01_72       1_72     COUNT

    0          0          0          1         876
    0          0          1          0       19628
    0          0          1          1        2532   <- hypomethylated-CHH
    0          1          0          0        7300
    0          1          1          0        2585
    0          1          1          1          96
    1          0          0          0         814
    1          0          0          1          14
    1          0          1          1           5
    1          1          0          0         295   <- hypermethylated-CHH
    1          1          1          0          18
    1          1          1          1           2



 */


/* PRep and export data for the following LINE PLOTS:

(1) Average accessibility in hyper/hypo DARs
(1) Average accessibility in hyper/hypo DMRs
(2) TSS accessibility plots

*/


/* methylation in DARs */


* GC methylation data;
data meth_data;
  set arabMAP.methylation_data_cg_chg_chh;
  where site_type="&siteType.";
  keep chr stop_pos treatment rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by chr stop_pos treatment  rep ;
proc means data=meth_data noprint;
  by chr stop_pos treatment   ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;

data meth_data3;
  set meth_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=meth_data3 out=meth_sbys10;
  where total_C >= 10;
  by chr stop_pos;
  id treatment ;
  var perc_methyl2;
run;

data meth_sbys10_2;
  set meth_sbys10;
  if _01Gy ne . and _0Gy ne . then _01Gy_common=_01Gy; else _01Gy_common=.;
  if _1Gy ne . and _0Gy ne . then _1Gy_common=_1Gy; else _1Gy_common=.;
  if (_1Gy ne . or _01Gy ne .) and _0Gy ne . then _0Gy_common=_0Gy; else _0Gy_common=.;

  if _01Gy_common ne . and _1Gy_common ne .  then do;
    _01Gy_common_all = _01Gy_common;
    _1Gy_common_all = _1Gy_common;
    _0Gy_common_all = _0Gy_common;
    end;
  else do;
   _01Gy_common_all = .;
   _1Gy_common_all = .;
   _0Gy_common_all = . ;
    end;
  rename stop_pos=pos;
run;


/* Get DARs */


data up_dar_01_1kb up_dar_1_1kb dn_dar_01_1kb dn_dar_1_1kb;
     set arabMAP.results_by_dar_annot;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";

     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;

     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;


     keep comparison chr dar_Center plot_start plot_stop ;
run;

proc sort data=up_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;





%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=meth_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dar_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dar_01_1kb);
%mergeMETH(dn_dar_01_1kb);
%mergeMETH(up_dar_1_1kb);
%mergeMETH(dn_dar_1_1kb);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hyper_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hyper_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hyper_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hyper_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hyper_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hyper_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hyper_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hyper_DAR_1kb_1Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hypo_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hypo_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hypo_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hypo_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hypo_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hypo_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hypo_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hypo_DAR_1kb_1Gy_0Gy_common_binned);



/* Do the same but for DMRs */

/* Get DMRs */




data up_dmr_01_1kb up_dmr_1_1kb dn_dmr_01_1kb dn_dmr_1_1kb;
     set arabMAP.results_by_dmr_annot;
     length feature $32.;
     where site_type="&siteType.";
     feature = scan(annotation, 1, " ");

     dmr_center=int((region_start + region_stop) / 2);
     plot_start=dmr_center - 999;
     plot_stop=dmr_center + 999;

     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_1kb;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_1kb;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_1kb;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_1kb;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_1kb;

     keep comparison chr dmr_Center plot_start plot_stop ;
run;

proc sort data=up_dmr_01_1kb nodup;  by chr dmr_center plot_start plot_stop; run;
proc sort data=dn_dmr_01_1kb nodup;  by chr dmr_center plot_start plot_stop; run;
proc sort data=up_dmr_1_1kb nodup;  by chr dmr_center plot_start plot_stop; run;
proc sort data=dn_dmr_1_1kb nodup;  by chr dmr_center plot_start plot_stop; run;





%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dmr_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=meth_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dmr_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dmr_01_1kb);
%mergeMETH(dn_dmr_01_1kb);
%mergeMETH(up_dmr_1_1kb);
%mergeMETH(dn_dmr_1_1kb);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hyper_DMR_1kb_01Gy_0Gy);
%exportLine( mean_diff_up_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hyper_DMR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_up_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hyper_DMR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_up_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hyper_DMR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_up_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hyper_DMR_1kb_1Gy_0Gy);
%exportLine( mean_diff_up_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hyper_DMR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_up_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hyper_DMR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_up_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hyper_DMR_1kb_1Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hypo_DMR_1kb_01Gy_0Gy);
%exportLine( mean_diff_dn_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hypo_DMR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_dn_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hypo_DMR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_dn_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hypo_DMR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hypo_DMR_1kb_1Gy_0Gy);
%exportLine( mean_diff_dn_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hypo_DMR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_dn_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hypo_DMR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_dn_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hypo_DMR_1kb_1Gy_0Gy_common_binned);


/* As above, but now for TSSs */


data site2promoter;
  set arabMAP.results_by_dmc_annot;
  where site_type="&siteType.";
  length gene_id $20.;
  if abs(distance_to_tss) > 999 then delete;
  if count(nearest_promoterID, "-T1") > 0 then gene_ID=compress(upcase(tranwrd(Nearest_PromoterID,"-T1","")));
  else gene_ID=compress(upcase(scan(Nearest_PromoterID,1,".")));
  keep gene_ID chr start_pos stop_pos strand distance_to_tss nearest_promoterID;
  rename stop_pos=pos nearest_promoterID=transcript_id;
run;

proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/Arabidopsis_thaliana.TAIR10.37.gtf"
  out=gtf dbms=tab replace;
  guessingrows=max;
  getnames=no;
run;

data gene;
  set gtf;
  where VAR3 = "gene";
  length gene_id $15.;
  gene_id=compress(tranwrd(tranwrd(scan(VAR9,2," "), ";", ""), '"', ''));
  keep gene_id VAR7; 
  rename VAR7=strand;
run;

proc sort data=gene nodup;
  by gene_id;
proc sort data=site2promoter nodup;
    by gene_id;
run;

data site2promoter2 no_strand no_site;
  merge site2promoter (in=in1) gene (in=in2);
  by gene_id;
  if in1 and in2 then output site2promoter2;
  else if in1 then output no_strand;
  else output no_site;
run;

/* distance to TSS is relative to strand of transcript, so I don't need to flip anything!!! */

proc sort data=site2promoter2 nodup;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data site2promoter2_w_meth;
  merge site2promoter2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data site2promoter2_w_meth2;
  set site2promoter2_w_meth;
  grouped_pos=int(distance_to_TSS/10) * 10;
  if distance_to_TSS < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=site2promoter2_w_meth2;
  by distance_to_TSS grouped_pos2 ;
run;


proc means data=site2promoter2_w_meth2 noprint;
  by distance_to_TSS grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_tss
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_tss ;
  by grouped_pos2;
run;


proc means data=mean_diff_tss noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_tss_2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_tss_1;
  set mean_diff_tss;
  drop _TYPE_ _FREQ_;
  keep distance_To_tss mean_: ;
  rename distance_to_tss=pos;
run;

data mean_diff_tss_3;
  set mean_diff_tss_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;


/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, var3, outName);

data export;
  retain pos &var1. &var2. &var3.;
  set &inData.;
  keep pos &var1. &var2. &var3.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_tss_1, mean_0Gy_common, mean_01Gy_common, mean_1Gy_common, TSS_1kb);
%exportLine( mean_diff_tss_3, mean_0Gy_common, mean_01Gy_common, mean_1Gy_common, TSS_1kb_binned);

%exportLine( mean_diff_tss_1, mean_0Gy_common_all, mean_01Gy_common_all, mean_1Gy_common_all, TSS_1kb_common);
%exportLine( mean_diff_tss_3, mean_0Gy_common_all, mean_01Gy_common_all, mean_1Gy_common_all, TSS_1kb_common_binned);


%mend;

%methDataGen(CG);
%methDataGen(CHG);
%methDataGen(CHH);



/* scatterplots :
  (1) rep to rep concordance count shared/unique, pearson correlation on perc methyl (all, shared)
  (2) HCG (mean) shared only, pearson correlation
  */

%macro concordance(siteType,treatment,units);

%if &siteType.=GC %then %do;

data meth_data_rep1;
   set arabMAP.methylation_data_gc;
   where rep=1 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   keep chr stop_pos perc_methyl_norm total_C;
   rename stop_pos=pos perc_methyl_norm=perc_methyl_&treatment._&units._1 total_C=total_C_1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_gc;
   where rep=2 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   keep chr stop_pos perc_methyl_norm total_C;
   rename stop_pos=pos perc_methyl_norm=perc_methyl_&treatment._&units._2 total_C=total_C_2;
run;

%end;

%else %do;

data meth_data_rep1;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=1 and site_type="&siteType." and treatment="&treatment."  ;
   keep chr stop_pos perc_methyl total_C;
   rename stop_pos=pos perc_methyl=perc_methyl_&treatment._&units._1 total_C=total_C_1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=2 and site_type="&siteType." and treatment="&treatment."  ;
   keep chr stop_pos perc_methyl total_C;
   rename stop_pos=pos perc_methyl=perc_methyl_&treatment._&units._2 total_C=total_C_2;
run;

%end;

proc sort data=meth_data_rep1;
  by chr pos;
proc sort data=meth_data_rep2;
  by chr pos;
run;

data meth_data_compare;
   merge meth_data_rep1 (in=in1) meth_Data_rep2 (in=in2);
  by chr pos;
  if in1 then flag_rep1_&treatment._&units.=1;
  else do; flag_rep1_&treatment._&units.=0; perc_methyl_&treatment._&units._1=0; total_c_1=0; end;
  if in2 then flag_rep2_&treatment._&units.=1;
  else do; flag_rep2_&treatment._&units.=0; perc_methyl_&treatment._&units._2=0; total_c_1=0; end;
run;

data meth_data_compare2;
  set meth_data_compare;
  max_C = total_C_1 + total_C_2;
  if max_C >= 10 then flag_keep=1; else flag_keep=0;
  run;

proc freq data=meth_data_compare2;
   where flag_keep=1;
   tables flag_rep1_&treatment._&units. * flag_rep2_&treatment._&units. ; 
run;

proc corr data=meth_data_compare2 pearson;
  where total_C_1 > 0 and total_C_2 > 0;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;

proc corr data=meth_data_compare2 pearson;
  where flag_keep=1 and flag_rep1_&treatment._&units. =1 and flag_rep2_&treatment._&units. = 1 ;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;

data export_Data;
   set meth_data_compare2;
   keep chr pos perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_&siteType._&treatment._&units..csv"
     dbms=csv replace;
run;

%mend;


%concordance(CG, 0Gy, 0U);
/*
      flag_rep1_0Gy_0U
             flag_rep2_0Gy_0U

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       1|  Total
   ---------+--------+
          0 |   1128 |   1128
            |   0.03 |   0.03
            | 100.00 |
            |   0.03 |
   ---------+--------+
          1 |3924221 |3924221
            |  99.97 |  99.97
            | 100.00 |
            |  99.97 |
   ---------+--------+
   Total     3925349  3925349
              100.00   100.00

         The SAS System
  
Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

perc_methyl_0Gy_0U_1     5154586       0.30002       0.40726       1546473             0       1.00000
perc_methyl_0Gy_0U_2     5154586       0.30033       0.40696       1548099             0       1.00000


                            Pearson Correlation Coefficients, N = 5154586
                                      Prob > |r| under H0: Rho=0

                                                 perc_methyl_      perc_methyl_
                                                     0Gy_0U_1          0Gy_0U_2

                       perc_methyl_0Gy_0U_1           1.00000           0.94429
                                                                         <.0001

                       perc_methyl_0Gy_0U_2           0.94429           1.00000
                                                       <.0001                            

   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_0Gy_0U_1     3924221       0.33743       0.41563       1324132             0       1.00000
   perc_methyl_0Gy_0U_2     3924221       0.33727       0.41543       1323525             0       1.00000


                               Pearson Correlation Coefficients, N = 3924221
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        0Gy_0U_1          0Gy_0U_2

                          perc_methyl_0Gy_0U_1           1.00000           0.96658
                                                                            <.0001

                          perc_methyl_0Gy_0U_2           0.96658           1.00000
                                                          <.0001


*/
%concordance(CG, 01Gy, 0U);
/*
 flag_rep1_01Gy_0U
           flag_rep2_01Gy_0U

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        0 |    331 |    331
          |   0.01 |   0.01
          | 100.00 |
          |   0.01 |
 ---------+--------+
        1 |4489853 |4489853
          |  99.99 |  99.99
          | 100.00 |
          |  99.99 |
 ---------+--------+
 Total     4490184  4490184
            100.00   100.00



    Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_01Gy_0U_1     5316165       0.29626       0.39809       1574957             0       1.00000
    perc_methyl_01Gy_0U_2     5316165       0.29402       0.39736       1563068             0       1.00000


                                Pearson Correlation Coefficients, N = 5316165
                                          Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                         01Gy_0U_1         01Gy_0U_2

                           perc_methyl_01Gy_0U_1           1.00000           0.94806
                                                                              <.0001

                           perc_methyl_01Gy_0U_2           0.94806           1.00000
                                                            <.0001



   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1     4489853       0.31312       0.40114       1405852             0       1.00000
   perc_methyl_01Gy_0U_2     4489853       0.31086       0.40033       1395720             0       1.00000


                               Pearson Correlation Coefficients, N = 4489853
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.96354
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.96354           1.00000
                                                           <.0001
*/
%concordance(CG, 1Gy, 0U);
/*

    flag_rep1_1Gy_0U
              flag_rep2_1Gy_0U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           0 |   7961 |   7961
             |   0.22 |   0.22
             | 100.00 |
             |   0.22 |
    ---------+--------+
           1 |3560334 |3560334
             |  99.78 |  99.78
             | 100.00 |
             |  99.78 |
    ---------+--------+
    Total     3568295  3568295
               100.00   100.00


     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_1Gy_0U_1     5081067       0.29735       0.41020       1510831             0       1.00000
     perc_methyl_1Gy_0U_2     5081067       0.30495       0.40880       1549447             0       1.00000


                                 Pearson Correlation Coefficients, N = 5081067
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          1Gy_0U_1          1Gy_0U_2

                            perc_methyl_1Gy_0U_1           1.00000           0.92965
                                                                              <.0001

                            perc_methyl_1Gy_0U_2           0.92965           1.00000
                                                            <.0001

                                                 The SAS System                                          19

       Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

       perc_methyl_1Gy_0U_1     3560334       0.34248       0.41945       1219338             0       1.00000
       perc_methyl_1Gy_0U_2     3560334       0.34610       0.41779       1232244             0       1.00000


                                   Pearson Correlation Coefficients, N = 3560334
                                             Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                            1Gy_0U_1          1Gy_0U_2

                              perc_methyl_1Gy_0U_1           1.00000           0.96308
                                                                                <.0001

                              perc_methyl_1Gy_0U_2           0.96308           1.00000
                                                              <.0001



*/

%concordance(CHG, 0Gy, 0U);
/*

    flag_rep1_0Gy_0U
              flag_rep2_0Gy_0U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           0 |   1230 |   1230
             |   0.03 |   0.03
             | 100.00 |
             |   0.03 |
    ---------+--------+
           1 |4273540 |4273540
             |  99.97 |  99.97
             | 100.00 |
             |  99.97 |
    ---------+--------+
    Total     4274770  4274770
               100.00   100.00
      Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

      perc_methyl_0Gy_0U_1     5650661       0.12442       0.24155        703079             0       1.00000
      perc_methyl_0Gy_0U_2     5650661       0.12518       0.24309        707376             0       1.00000


                                  Pearson Correlation Coefficients, N = 5650661
                                            Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                           0Gy_0U_1          0Gy_0U_2

                             perc_methyl_0Gy_0U_1           1.00000           0.84020
                                                                               <.0001

                             perc_methyl_0Gy_0U_2           0.84020           1.00000
                                                             <.0001

                                                  The SAS System                                          19:4


                                                   Simple Statistics

         Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

         perc_methyl_0Gy_0U_1     4273540       0.14149       0.24878        604670             0       1.00000
         perc_methyl_0Gy_0U_2     4273540       0.14179       0.25037        605957             0       1.00000


                                     Pearson Correlation Coefficients, N = 4273540
                                               Prob > |r| under H0: Rho=0

                                                          perc_methyl_      perc_methyl_
                                                              0Gy_0U_1          0Gy_0U_2

                                perc_methyl_0Gy_0U_1           1.00000           0.90438
                                                                                  <.0001

                                perc_methyl_0Gy_0U_2           0.90438           1.00000
                                                                <.0001

                                                     The SAS System                                          19:40

*/
%concordance(CHG, 01Gy, 0U);
/*
      flag_rep1_01Gy_0U
                flag_rep2_01Gy_0U

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       1|  Total
      ---------+--------+
             0 |    331 |    331
               |   0.01 |   0.01
               | 100.00 |
               |   0.01 |
      ---------+--------+
             1 |4946525 |4946525
               |  99.99 |  99.99
               | 100.00 |
               |  99.99 |
      ---------+--------+
      Total     4946856  4946856
                 100.00   100.00



                                              Simple Statistics

   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1     5837991       0.12251       0.22717        715213             0       1.00000
   perc_methyl_01Gy_0U_2     5837991       0.12044       0.22753        703117             0       1.00000


                               Pearson Correlation Coefficients, N = 5837991
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.84203
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.84203           1.00000
                                                           <.0001

                                               Simple Statistics

    Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_01Gy_0U_1     4946525       0.12825       0.22627        634400             0       1.00000
    perc_methyl_01Gy_0U_2     4946525       0.12620       0.22650        624268             0       1.00000


                                Pearson Correlation Coefficients, N = 4946525
                                          Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                         01Gy_0U_1         01Gy_0U_2

                           perc_methyl_01Gy_0U_1           1.00000           0.88876
                                                                              <.0001

                           perc_methyl_01Gy_0U_2           0.88876           1.00000
                                                            <.0001

                                                 The SAS System                                          19:40


*/
%concordance(CHG, 1Gy, 0U);
/*
  flag_rep1_1Gy_0U
            flag_rep2_1Gy_0U

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       1|  Total
  ---------+--------+
         0 |   7888 |   7888
           |   0.20 |   0.20
           | 100.00 |
           |   0.20 |
  ---------+--------+
         1 |3874859 |3874859
           |  99.80 |  99.80
           | 100.00 |
           |  99.80 |
  ---------+--------+
  Total     3882747  3882747
             100.00   100.00

        The SAS System


                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_1Gy_0U_1     5575382       0.12117       0.24744        675590             0       1.00000
  perc_methyl_1Gy_0U_2     5575382       0.13336       0.25278        743518             0       1.00000


                              Pearson Correlation Coefficients, N = 5575382
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       1Gy_0U_1          1Gy_0U_2

                         perc_methyl_1Gy_0U_1           1.00000           0.79521
                                                                           <.0001

                         perc_methyl_1Gy_0U_2           0.79521           1.00000
                                                         <.0001

                                                Simple Statistics

   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_1Gy_0U_1     3874859       0.14311       0.25683        554515             0       1.00000
   perc_methyl_1Gy_0U_2     3874859       0.15000       0.25710        581233             0       1.00000


                               Pearson Correlation Coefficients, N = 3874859
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        1Gy_0U_1          1Gy_0U_2

                          perc_methyl_1Gy_0U_1           1.00000           0.88591
                                                                            <.0001

                          perc_methyl_1Gy_0U_2           0.88591           1.00000
                                                          <.0001                                           The SAS System                                          19:40 S

*/

%concordance(CHH, 0Gy, 0U);
/*

  flag_rep1_0Gy_0U
            flag_rep2_0Gy_0U

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       1|  Total
  ---------+--------+
         0 |  10311 |  10311
           |   0.06 |   0.06
           | 100.00 |
           |   0.06 |
  ---------+--------+
         1 |1.829E7 |1.829E7
           |  99.94 |  99.94
           | 100.00 |
           |  99.94 |
  ---------+--------+
  Total     1.831E7  1.831E7
             100.00   100.00

                                              Simple Statistics

    Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_0Gy_0U_1    27234750       0.06308       0.14888       1717888             0       1.00000
    perc_methyl_0Gy_0U_2    27234750       0.06332       0.14614       1724404             0       1.00000


                               Pearson Correlation Coefficients, N = 27234750
                                          Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                         0Gy_0U_1          0Gy_0U_2

                           perc_methyl_0Gy_0U_1           1.00000           0.53149
                                                                             <.0001

                           perc_methyl_0Gy_0U_2           0.53149           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 S



                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_0Gy_0U_1    18294791       0.07021       0.14067       1284486             0       1.00000
  perc_methyl_0Gy_0U_2    18294791       0.06908       0.13745       1263805             0       1.00000


                             Pearson Correlation Coefficients, N = 18294791
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       0Gy_0U_1          0Gy_0U_2

                         perc_methyl_0Gy_0U_1           1.00000           0.72832
                                                                           <.0001

                         perc_methyl_0Gy_0U_2           0.72832           1.00000
                                                         <.0001

*/
%concordance(CHH, 01Gy, 0U);
/*

       flag_rep1_01Gy_0U
                 flag_rep2_01Gy_0U

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       1|  Total
       ---------+--------+
              0 |   3249 |   3249
                |   0.01 |   0.01
                | 100.00 |
                |   0.01 |
       ---------+--------+
              1 | 2.37E7 | 2.37E7
                |  99.99 |  99.99
                | 100.00 |
                |  99.99 |
       ---------+--------+
       Total      2.37E7   2.37E7
                  100.00   100.00

             The SAS System


                                                   Simple Statistics

        Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

        perc_methyl_01Gy_0U_1    29288053       0.06168       0.12649       1806543             0       1.00000
        perc_methyl_01Gy_0U_2    29288053       0.05897       0.12353       1727175             0       1.00000


                                    Pearson Correlation Coefficients, N = 29288053
                                              Prob > |r| under H0: Rho=0

                                                          perc_methyl_      perc_methyl_
                                                             01Gy_0U_1         01Gy_0U_2

                               perc_methyl_01Gy_0U_1           1.00000           0.51221
                                                                                  <.0001

                               perc_methyl_01Gy_0U_2           0.51221           1.00000
                                                                <.0001

                                                     The SAS System                                          19:

                                                   The CORR Procedure

                                                 Simple Statistics

      Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

      perc_methyl_01Gy_0U_1    23696673       0.06363       0.11653       1507771             0       1.00000
      perc_methyl_01Gy_0U_2    23696673       0.06093       0.11436       1443834             0       1.00000


                                  Pearson Correlation Coefficients, N = 23696673
                                            Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                           01Gy_0U_1         01Gy_0U_2

                             perc_methyl_01Gy_0U_1           1.00000           0.64156
                                                                                <.0001

                             perc_methyl_01Gy_0U_2           0.64156           1.00000
                                                              <.0001

                                                   The SAS System                                          19:40 S

*/
%concordance(CHH, 1Gy, 0U);
/*

      flag_rep1_1Gy_0U
                flag_rep2_1Gy_0U

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       1|  Total
      ---------+--------+
             0 |  76558 |  76558
               |   0.44 |   0.44
               | 100.00 |
               |   0.44 |
      ---------+--------+
             1 | 1.73E7 | 1.73E7
               |  99.56 |  99.56
               | 100.00 |
               |  99.56 |
      ---------+--------+
      Total     1.738E7  1.738E7
                 100.00   100.00

            The SAS System



                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_1Gy_0U_1    27286581       0.05769       0.14796       1574213             0       1.00000
  perc_methyl_1Gy_0U_2    27286581       0.09236       0.19675       2520142             0       1.00000


                             Pearson Correlation Coefficients, N = 27286581
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       1Gy_0U_1          1Gy_0U_2

                         perc_methyl_1Gy_0U_1           1.00000           0.41093
                                                                           <.0001

                         perc_methyl_1Gy_0U_2           0.41093           1.00000
                                                         <.0001

                                              Simple Statistics

    Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_0U_1    17299403       0.06694       0.14378       1158074             0       1.00000
    perc_methyl_1Gy_0U_2    17299403       0.08553       0.16447       1479686             0       1.00000


                               Pearson Correlation Coefficients, N = 17299403
                                          Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                         1Gy_0U_1          1Gy_0U_2

                           perc_methyl_1Gy_0U_1           1.00000           0.64831
                                                                             <.0001

                           perc_methyl_1Gy_0U_2           0.64831           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 S

*/

%concordance(GC, 0Gy, 0U);
/*

 flag_rep1_0Gy_0U
           flag_rep2_0Gy_0U

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        1 |4996963 |4996963
          | 100.00 | 100.00
          | 100.00 |
          | 100.00 |
 ---------+--------+
 Total     4996963  4996963
            100.00   100.00

       The SAS System



          Pearson Correlation Coefficients, N = 4996963
                    Prob > |r| under H0: Rho=0

                               perc_methyl_      perc_methyl_
                                   0Gy_0U_1          0Gy_0U_2

     perc_methyl_0Gy_0U_1           1.00000           0.90848
                                                       <.0001

     perc_methyl_0Gy_0U_2           0.90848           1.00000
                                     <.0001

                          The SAS System                                          19:4
*/
%concordance(GC, 01Gy, 0U);
/*

      flag_rep1_01Gy_0U
                flag_rep2_01Gy_0U

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       1|  Total
      ---------+--------+
             1 |5803543 |5803543
               | 100.00 | 100.00
               | 100.00 |
               | 100.00 |
      ---------+--------+
      Total     5803543  5803543
                 100.00   100.00

                                                Simple Statistics

     Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_01Gy_0U_1     5803543       0.09914       0.20332        575387             0       0.96874
     perc_methyl_01Gy_0U_2     5803543       0.10323       0.21691        599089             0       1.03335


                                 Pearson Correlation Coefficients, N = 5803543
                                           Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                          01Gy_0U_1         01Gy_0U_2

                            perc_methyl_01Gy_0U_1           1.00000           0.89538
                                                                               <.0001

                            perc_methyl_01Gy_0U_2           0.89538           1.00000
                                                             <.0001

                                                  The SAS System                                          19:40

*/
%concordance(GC, 1Gy, 0U);
/*

    flag_rep1_1Gy_0U
              flag_rep2_1Gy_0U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           1 |4539526 |4539526
             | 100.00 | 100.00
             | 100.00 |
             | 100.00 |
    ---------+--------+
    Total     4539526  4539526
               100.00   100.00


   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_1Gy_0U_1     4539526       0.11949       0.25791        542407             0       1.10532
   perc_methyl_1Gy_0U_2     4539526       0.10815       0.21407        490956             0       0.91300


                               Pearson Correlation Coefficients, N = 4539526
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        1Gy_0U_1          1Gy_0U_2

                          perc_methyl_1Gy_0U_1           1.00000           0.88470
                                                                            <.0001

                          perc_methyl_1Gy_0U_2           0.88470           1.00000
                                                          <.0001



                                           Simple Statistics

 Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

 perc_methyl_1Gy_0U_1     4539526       0.11949       0.25791        542407             0       1.10532
 perc_methyl_1Gy_0U_2     4539526       0.10815       0.21407        490956             0       0.91300


                             Pearson Correlation Coefficients, N = 4539526
                                       Prob > |r| under H0: Rho=0

                                                  perc_methyl_      perc_methyl_
                                                      1Gy_0U_1          1Gy_0U_2

                        perc_methyl_1Gy_0U_1           1.00000           0.88470
                                                                          <.0001

                        perc_methyl_1Gy_0U_2           0.88470           1.00000
                                                        <.0001

                                             The SAS System                                          19:40 Sunda
*/

%concordance(GC, 0Gy, 100U);
/*

         flag_rep1_0Gy_100U
                   flag_rep2_0Gy_100U

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       1|  Total
         ---------+--------+
                1 |4712204 |4712204
                  | 100.00 | 100.00
                  | 100.00 |
                  | 100.00 |
         ---------+--------+
         Total     4712204  4712204
                    100.00   100.00



                                                Simple Statistics

     Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_100U_1     4712204       0.29215       0.21867       1376651             0       0.99418
     perc_methyl_0Gy_100U_2     4712204       0.29735       0.23239       1401161             0       1.00588


                                  Pearson Correlation Coefficients, N = 4712204
                                            Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                          0Gy_100U_1        0Gy_100U_2

                            perc_methyl_0Gy_100U_1           1.00000           0.60706
                                                                                <.0001

                            perc_methyl_0Gy_100U_2           0.60706           1.00000
                                                              <.0001

*/
%concordance(GC, 01Gy, 100U);
/*

    flag_rep1_01Gy_100U
              flag_rep2_01Gy_100U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           1 |5988381 |5988381
             | 100.00 | 100.00
             | 100.00 |
             | 100.00 |
    ---------+--------+
    Total     5988381  5988381
               100.00   100.00

          The SAS System

                                               Simple Statistics

   Variable                          N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_100U_1     5988381       0.44072       0.20608       2639227             0       1.06212
   perc_methyl_01Gy_100U_2     5988381       0.44296       0.18222       2652589             0       0.97123


                                Pearson Correlation Coefficients, N = 5988381
                                          Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                        01Gy_100U_1       01Gy_100U_2

                          perc_methyl_01Gy_100U_1           1.00000           0.46823
                                                                               <.0001

                          perc_methyl_01Gy_100U_2           0.46823           1.00000
                                                             <.0001

*/
%concordance(GC, 1Gy, 100U);
/*
          flag_rep1_1Gy_100U
                    flag_rep2_1Gy_100U

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       1|  Total
          ---------+--------+
                 1 |4703751 |4703751
                   | 100.00 | 100.00
                   | 100.00 |
                   | 100.00 |
          ---------+--------+
          Total     4703751  4703751
                     100.00   100.00

                                               Simple Statistics

    Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_100U_1     4703751       0.34516       0.21316       1623528             0       0.97086
    perc_methyl_1Gy_100U_2     4703751       0.34052       0.22484       1601728             0       1.03094


                                 Pearson Correlation Coefficients, N = 4703751
                                           Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                         1Gy_100U_1        1Gy_100U_2

                           perc_methyl_1Gy_100U_1           1.00000           0.55288
                                                                               <.0001

                           perc_methyl_1Gy_100U_2           0.55288           1.00000
                                                             <.0001



*/




/* 100X sites only and bin into 100bp windows */



%macro concordance(siteType,treatment,units);

%if &siteType.=GC %then %do;

data meth_data_rep1;
   set arabMAP.methylation_data_gc;
   where rep=1 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads; 
   rename  perc_methyl2=perc_methyl_&treatment._&units._1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_gc;
   where rep=2 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads;
   rename perc_methyl2=perc_methyl_&treatment._&units._2;
run;

%end;

%else %do;

data meth_data_rep1;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=1 and site_type="&siteType." and treatment="&treatment."  ;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads;
   rename perc_methyl2=perc_methyl_&treatment._&units._1 ;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=2 and site_type="&siteType." and treatment="&treatment."  ;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads;
   rename perc_methyl2=perc_methyl_&treatment._&units._2 ;
run;

%end;



proc sort data=meth_data_rep1;
  by chr chr_bin;
proc sort data=meth_data_rep2;
  by chr chr_bin;
run;

proc means data=meth_data_rep1 noprint;
 by chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl_&treatment._&units._1;
  output out=mean_methyl_by_bin_rep1 (drop=_TYPE_ _FREQ_) mean=;
run;

proc means data=meth_data_rep2 noprint;
 by chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl_&treatment._&units._2;
  output out=mean_methyl_by_bin_rep2 (drop=_TYPE_ _FREQ_) mean=;
run;

proc sort data=mean_methyl_by_bin_rep1;
  by chr chr_bin;
proc sort data=mean_methyl_by_bin_rep2;
  by chr chr_bin;
run;


data meth_data_compare;
   merge mean_methyl_by_bin_rep1 (in=in1) mean_methyl_by_bin_rep2 (in=in2);
  by chr chr_bin;
  if in1 then flag_rep1_&treatment._&units.=1;
  else do; flag_rep1_&treatment._&units.=0; perc_methyl_&treatment._&units._1=0; total_c_1=0; end;
  if in2 then flag_rep2_&treatment._&units.=1;
  else do; flag_rep2_&treatment._&units.=0; perc_methyl_&treatment._&units._2=0; total_c_1=0; end;
run;

proc freq data=meth_data_compare;
   tables flag_rep1_&treatment._&units. * flag_rep2_&treatment._&units. ; 
run;

proc corr data=meth_data_compare pearson;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;


proc corr data=meth_data_compare pearson;
  where flag_rep1_&treatment._&units. =1 and flag_rep2_&treatment._&units. = 1 ;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;

data export_Data;
   set meth_data_compare;
   keep chr pos chr_bin perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_&siteType._&treatment._&units._100X_100bp_binned.csv"
     dbms=csv replace;
run;

%mend;


%concordance(CG, 0Gy, 0U);
/*


                                               Simple Statistics

     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_0U_1       13009      50.14782      39.86885        652373             0     100.00000
     perc_methyl_0Gy_0U_2       13009      50.15524      39.67361        652469             0     100.00000


                                  Pearson Correlation Coefficients, N = 13009
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          0Gy_0U_1          0Gy_0U_2

                            perc_methyl_0Gy_0U_1           1.00000           0.99253
                                                                              <.0001

                            perc_methyl_0Gy_0U_2           0.99253           1.00000
                                                            <.0001


*/


%concordance(CG, 01Gy, 0U);
/*

                                               Simple Statistics

    Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_01Gy_0U_1        2636      34.14883      31.56339         90016             0      95.40990
    perc_methyl_01Gy_0U_2        2636      33.51767      31.61295         88353             0      95.80665


                                  Pearson Correlation Coefficients, N = 2636
                                          Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                         01Gy_0U_1         01Gy_0U_2

                           perc_methyl_01Gy_0U_1           1.00000           0.99318
                                                                              <.0001

                           perc_methyl_01Gy_0U_2           0.99318           1.00000
                                                            <.0001




*/


%concordance(CG, 1Gy, 0U);
/*


   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_1Gy_0U_1        4195      40.39137      37.47797        169442             0     100.00000
   perc_methyl_1Gy_0U_2        4195      41.21943      37.27794        172915             0      97.34594


                                Pearson Correlation Coefficients, N = 4195
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        1Gy_0U_1          1Gy_0U_2

                          perc_methyl_1Gy_0U_1           1.00000           0.99080
                                                                            <.0001

                          perc_methyl_1Gy_0U_2           0.99080           1.00000
                                                          <.0001




*/


%concordance(CHG, 0Gy, 0U);
/*



                                               Simple Statistics

     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_0U_1       12966      27.31133      28.11818        354119             0     100.00000
     perc_methyl_0Gy_0U_2       12966      27.37780      27.99949        354981             0      96.26168


                                  Pearson Correlation Coefficients, N = 12966
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          0Gy_0U_1          0Gy_0U_2

                            perc_methyl_0Gy_0U_1           1.00000           0.98143
                                                                              <.0001

                            perc_methyl_0Gy_0U_2           0.98143           1.00000
                                                            <.0001

                                                 The SAS System                                          19:40

*/


%concordance(CHG, 01Gy, 0U);
/*

                                             Simple Statistics

  Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_01Gy_0U_1        2427      16.73066      17.18876         40605             0      95.12805
  perc_methyl_01Gy_0U_2        2427      16.41853      17.29947         39848             0      98.71795


                                Pearson Correlation Coefficients, N = 2427
                                        Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                       01Gy_0U_1         01Gy_0U_2

                         perc_methyl_01Gy_0U_1           1.00000           0.97850
                                                                            <.0001

                         perc_methyl_01Gy_0U_2           0.97850           1.00000
                                                          <.0001

                                               The SAS System                                          19:40 Su



*/


%concordance(CHG, 1Gy, 0U);
/*


                                               Simple Statistics

     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_1Gy_0U_1        3872      18.80897      22.31672         72828             0      93.16239
     perc_methyl_1Gy_0U_2        3872      19.33267      21.92126         74856             0      93.03798


                                  Pearson Correlation Coefficients, N = 3872
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          1Gy_0U_1          1Gy_0U_2

                            perc_methyl_1Gy_0U_1           1.00000           0.96834
                                                                              <.0001

                            perc_methyl_1Gy_0U_2           0.96834           1.00000
                                                            <.0001

                                                 The SAS System                                          19:40 Su


*/


%concordance(CHH, 0Gy, 0U);
/*


                                          Simple Statistics

Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

perc_methyl_0Gy_0U_1       16905      11.81273      12.42799        199694             0      86.13861
perc_methyl_0Gy_0U_2       16905      11.73605      12.32780        198398             0      88.11881


                             Pearson Correlation Coefficients, N = 16905
                                      Prob > |r| under H0: Rho=0

                                                 perc_methyl_      perc_methyl_
                                                     0Gy_0U_1          0Gy_0U_2

                       perc_methyl_0Gy_0U_1           1.00000           0.91475
                                                                         <.0001

                       perc_methyl_0Gy_0U_2           0.91475           1.00000
                                                       <.0001



*/


%concordance(CHH, 01Gy, 0U);
/*

                                              Simple Statistics

   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1        3125       8.53391       8.06098         26668       0.88496      96.10853
   perc_methyl_01Gy_0U_2        3125       8.38016       8.32528         26188             0      96.25674


                                 Pearson Correlation Coefficients, N = 3125
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.96061
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.96061           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 S



*/


%concordance(CHH, 1Gy, 0U);
/*


                                              Simple Statistics

    Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_0U_1        5329       9.39021      10.68359         50040             0      84.61539
    perc_methyl_1Gy_0U_2        5329       9.83599      10.07930         52416             0      79.41176


                                 Pearson Correlation Coefficients, N = 5329
                                          Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                         1Gy_0U_1          1Gy_0U_2

                           perc_methyl_1Gy_0U_1           1.00000           0.85068
                                                                             <.0001

                           perc_methyl_1Gy_0U_2           0.85068           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 Sun


*/


%concordance(GC, 0Gy, 0U);
/*


                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_0Gy_0U_1       13935      20.65310      23.55993        287801             0      93.50084
  perc_methyl_0Gy_0U_2       13935      22.98440      27.00505        320288             0     104.75641


                               Pearson Correlation Coefficients, N = 13935
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       0Gy_0U_1          0Gy_0U_2

                         perc_methyl_0Gy_0U_1           1.00000           0.94378
                                                                           <.0001

                         perc_methyl_0Gy_0U_2           0.94378           1.00000
                                                         <.0001

                                              The SAS System                                          19:40


*/


%concordance(GC, 01Gy, 0U);
/*


                                              Simple Statistics

   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1        2797      12.56083      13.82068         35133             0      92.79118
   perc_methyl_01Gy_0U_2        2797      13.15119      15.19590         36784             0     100.11650


                                 Pearson Correlation Coefficients, N = 2797
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.96533
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.96533           1.00000
                                                           <.0001


*/


%concordance(GC, 1Gy, 0U);
/*


                                          Simple Statistics

Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

perc_methyl_1Gy_0U_1        4416      17.78988      24.49459         78560             0     109.49917
perc_methyl_1Gy_0U_2        4416      14.69169      18.10247         64879             0      89.10914


                             Pearson Correlation Coefficients, N = 4416
                                      Prob > |r| under H0: Rho=0

                                                 perc_methyl_      perc_methyl_
                                                     1Gy_0U_1          1Gy_0U_2

                       perc_methyl_1Gy_0U_1           1.00000           0.92891
                                                                         <.0001

                       perc_methyl_1Gy_0U_2           0.92891           1.00000
                                                       <.0001


*/


%concordance(GC, 0Gy, 100U);
/*

                                                Simple Statistics

     Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_100U_1        7702      37.89239      21.87415        291847             0      97.44328
     perc_methyl_0Gy_100U_2        7702      38.89964      22.09051        299605             0      98.89780


                                   Pearson Correlation Coefficients, N = 7702
                                            Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                          0Gy_100U_1        0Gy_100U_2

                            perc_methyl_0Gy_100U_1           1.00000           0.93661
                                                                                <.0001

                            perc_methyl_0Gy_100U_2           0.93661           1.00000
                                                              <.0001



*/


%concordance(GC, 01Gy, 100U);
/*


                                             Simple Statistics

 Variable                          N          Mean       Std Dev           Sum       Minimum       Maximum

 perc_methyl_01Gy_100U_1        3406      63.06014      20.89049        214783             0     101.01176
 perc_methyl_01Gy_100U_2        3406      62.37673      19.94840        212455       0.39725      92.14254


                                Pearson Correlation Coefficients, N = 3406
                                        Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                      01Gy_100U_1       01Gy_100U_2

                        perc_methyl_01Gy_100U_1           1.00000           0.97556
                                                                             <.0001

                        perc_methyl_01Gy_100U_2           0.97556           1.00000
                                                           <.0001

*/


%concordance(GC, 1Gy, 100U);

/*

                                               Simple Statistics

    Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_100U_1        8295      40.27541      18.96855        334085             0      95.28848
    perc_methyl_1Gy_100U_2        8295      41.41053      20.28028        343500             0     102.18158


                                  Pearson Correlation Coefficients, N = 8295
                                           Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                         1Gy_100U_1        1Gy_100U_2

                           perc_methyl_1Gy_100U_1           1.00000           0.93477
                                                                               <.0001

                           perc_methyl_1Gy_100U_2           0.93477           1.00000
                                                             <.0001





*/





