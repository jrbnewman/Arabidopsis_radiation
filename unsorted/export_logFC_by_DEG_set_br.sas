
/* scatter plot of dose by FC for any significant  gene */

data de_results;
  set brassRNA.results_by_gene;
  keep gene_id flag_Mock_1hr_on flag_Mock_72hr_on
  flag_01gy_1hr_on  flag_01gy_72hr_on
  flag_1gy_1hr_on  flag_1gy_72hr_on
  flag_1_4cgy_1hr_on  flag_1_4cgy_72hr_on
  mean_cpm_: 
  fdr_: 
  flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_72h_fdr05
  flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_72h_fdr05
  flag_1_4cgy_v_Mock_1h_fdr05 flag_1_4cgy_v_Mock_72h_fdr05 ;
run;

proc contents data=de_results;run;quit;


data de_results2;
  set de_results;
  log2fc_01gy_Mock_1h = log2(mean_cpm_01gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_01gy_Mock_72h = log2(mean_cpm_01gy_72h ) - log2(mean_cpm_Mock_72h );
  log2fc_1gy_Mock_1h = log2(mean_cpm_1gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_1gy_Mock_72h = log2(mean_cpm_1gy_72h ) - log2(mean_cpm_Mock_72h );
  log2fc_1_4cgy_Mock_1h = log2(mean_cpm_1_4cgy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_1_4cgy_Mock_72h = log2(mean_cpm_1_4cgy_72h ) - log2(mean_cpm_Mock_72h );
run;



data up_01_1_fc1  up_01_72_fc1
     dn_01_1_fc1  dn_01_72_fc1
     up_1_1_fc1    up_1_72_fc1
     dn_1_1_fc1    dn_1_72_fc1
     up_1p4_1_fc1    up_1p4_72_fc1
     dn_1p4_1_fc1    dn_1p4_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;

if ((flag_1_4cgy_v_mock_1h_fdr05=1 and abs(log2fc_1_4cgy_Mock_1h) >= 1) or (flag_1_4cgy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1_4cgy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1_4cgy_1h > mean_cpm_Mock_1h) then output up_1p4_1_fc1;
if ((flag_1_4cgy_v_mock_1h_fdr05=1 and abs(log2fc_1_4cgy_Mock_1h) >= 1) or (flag_1_4cgy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1_4cgy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1_4cgy_1h < mean_cpm_Mock_1h) then output dn_1p4_1_fc1;
if ((flag_1_4cgy_v_mock_72h_fdr05=1 and abs(log2fc_1_4cgy_Mock_72h) >= 1) or (flag_1_4cgy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1_4cgy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1_4cgy_72h > mean_cpm_Mock_72h) then output up_1p4_72_fc1;
if ((flag_1_4cgy_v_mock_72h_fdr05=1 and abs(log2fc_1_4cgy_Mock_72h) >= 1) or (flag_1_4cgy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1_4cgy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1_4cgy_72h < mean_cpm_Mock_72h) then output dn_1p4_72_fc1;

keep gene_id;
run;


proc sort data=up_01_1_fc1 nodup; by gene_id; run;
proc sort data=up_01_72_fc1 nodup; by gene_id; run;
proc sort data=dn_01_1_fc1 nodup; by gene_id; run;
proc sort data=dn_01_72_fc1 nodup; by gene_id; run;
proc sort data=up_1_1_fc1 nodup; by gene_id; run;
proc sort data=up_1_72_fc1 nodup; by gene_id; run;
proc sort data=dn_1_1_fc1 nodup; by gene_id; run;
proc sort data=dn_1_72_fc1 nodup; by gene_id; run;
proc sort data=up_1p4_1_fc1 nodup; by gene_id; run;
proc sort data=up_1p4_72_fc1 nodup; by gene_id; run;
proc sort data=dn_1p4_1_fc1 nodup; by gene_id; run;
proc sort data=dn_1p4_72_fc1 nodup; by gene_id; run;


proc sort data=de_Results2;
  by gene_id;
run;

data data_up_01_1;
  merge up_01_1_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_01gy_Mock_1h log2fc_01gy_Mock_72h;
run;

data data_up_01_72;
  merge up_01_72_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_01gy_Mock_1h log2fc_01gy_Mock_72h;
run;

data data_dn_01_1;
  merge dn_01_1_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_01gy_Mock_1h log2fc_01gy_Mock_72h;
run;

data data_dn_01_72;
  merge dn_01_72_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_01gy_Mock_1h log2fc_01gy_Mock_72h;
run;


data data_up_1_1;
  merge up_1_1_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1gy_Mock_1h log2fc_1gy_Mock_72h;
run;

data data_up_1_72;
  merge up_1_72_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1gy_Mock_1h log2fc_1gy_Mock_72h;
run;

data data_dn_1_1;
  merge dn_1_1_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1gy_Mock_1h log2fc_1gy_Mock_72h;
run;

data data_dn_1_72;
  merge dn_1_72_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1gy_Mock_1h log2fc_1gy_Mock_72h;
run;



data data_up_1p4_1;
  merge up_1p4_1_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1_4cgy_Mock_1h log2fc_1_4cgy_Mock_72h;
run;

data data_up_1p4_72;
  merge up_1p4_72_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1_4cgy_Mock_1h log2fc_1_4cgy_Mock_72h;
run;

data data_dn_1p4_1;
  merge dn_1p4_1_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1_4cgy_Mock_1h log2fc_1_4cgy_Mock_72h;
run;

data data_dn_1p4_72;
  merge dn_1p4_72_fc1 (in=in1) de_Results2 (in=in2);
  by gene_id;
  if in1 and in2;
  keep gene_id log2fc_1_4cgy_Mock_1h log2fc_1_4cgy_Mock_72h;
run;


proc export data=data_up_01_1 outfile="!HOME/concannon/DTRA/Br_DEG_up_10cGy_1h.csv" dbms=csv replace; run;
proc export data=data_dn_01_1 outfile="!HOME/concannon/DTRA/Br_DEG_dn_10cGy_1h.csv" dbms=csv replace; run;
proc export data=data_up_01_72 outfile="!HOME/concannon/DTRA/Br_DEG_up_10cGy_72h.csv" dbms=csv replace; run;
proc export data=data_dn_01_72 outfile="!HOME/concannon/DTRA/Br_DEG_dn_10cGy_72h.csv" dbms=csv replace; run;

proc export data=data_up_1_1 outfile="!HOME/concannon/DTRA/Br_DEG_up_100cGy_1h.csv" dbms=csv replace; run;
proc export data=data_dn_1_1 outfile="!HOME/concannon/DTRA/Br_DEG_dn_100cGy_1h.csv" dbms=csv replace; run;
proc export data=data_up_1_72 outfile="!HOME/concannon/DTRA/Br_DEG_up_100cGy_72h.csv" dbms=csv replace; run;
proc export data=data_dn_1_72 outfile="!HOME/concannon/DTRA/Br_DEG_dn_100cGy_72h.csv" dbms=csv replace; run;

proc export data=data_up_1p4_1 outfile="!HOME/concannon/DTRA/Br_DEG_up_1p4cGy_1h.csv" dbms=csv replace; run;
proc export data=data_dn_1p4_1 outfile="!HOME/concannon/DTRA/Br_DEG_dn_1p4cGy_1h.csv" dbms=csv replace; run;
proc export data=data_up_1p4_72 outfile="!HOME/concannon/DTRA/Br_DEG_up_1p4cGy_72h.csv" dbms=csv replace; run;
proc export data=data_dn_1p4_72 outfile="!HOME/concannon/DTRA/Br_DEG_dn_1p4cGy_72h.csv" dbms=csv replace; run;




