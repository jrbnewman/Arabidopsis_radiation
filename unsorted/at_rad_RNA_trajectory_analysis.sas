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



data up_01_1_fc1  up_01_72_fc1
     dn_01_1_fc1  dn_01_72_fc1
     up_1_1_fc1    up_1_72_fc1
     dn_1_1_fc1    dn_1_72_fc1
     up_b_1_fc1    up_b_72_fc1
     dn_b_1_fc1    dn_b_72_fc1
     up_a_1_fc1    up_a_72_fc1
     dn_a_1_fc1    dn_a_72_fc1;

     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;

if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) 
and ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_b_1_fc1;

if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h)
and ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_b_1_fc1;

if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) 
and ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_b_72_fc1;

if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h)
and ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_b_72_fc1;

if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) 
or ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_a_1_fc1;

if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h)
or ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_a_1_fc1;

if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) 
or ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_a_72_fc1;

if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h)
or ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_a_72_fc1;

keep gene_id log2fc_: ;
run;


%macro prepFCdata(datain,geneSet);

proc transpose data=&datain. out=&datain._sbys;
   by gene_id;
   var log2fc_: ;
run;

data counts4model_&datain.;
   set &datain._sbys;
   format gene_set $20.;
   format comparison $20.;
   format hour best12.;
   gene_set = compress("&geneSet.");
   comparison = catt(scan(_NAME_,2,"_"),"_vs_0cGy");
   hour = tranwrd(scan(_NAME_,4,"_"),"h","") + 0 ;
   FC = 2 ** COL1;
   rename COL1=log2FC;
run;

%mend;


%prepFCdata(up_01_1_fc1,10cGy_upreg_1h);
%prepFCdata(up_01_72_fc1,10cGy_upreg_72h);
%prepFCdata(dn_01_1_fc1,10cGy_dnreg_1h);
%prepFCdata(dn_01_72_fc1,10cGy_dnreg_72h);
%prepFCdata(up_1_1_fc1,100cGy_upreg_1h);
%prepFCdata(up_1_72_fc1,100cGy_upreg_72h);
%prepFCdata(dn_1_1_fc1,100cGy_dnreg_1h);
%prepFCdata(dn_1_72_fc1,100cGy_dnreg_72h);
*%prepFCdata(up_b_1_fc1,both_upreg_1h);
*%prepFCdata(up_b_72_fc1,both_upreg_72h);
*%prepFCdata(dn_b_1_fc1,both_dnreg_1h);
*%prepFCdata(dn_b_72_fc1,both_dnreg_72h);
*%prepFCdata(up_a_1_fc1,any_upreg_1h);
*%prepFCdata(up_a_72_fc1,any_upreg_72h);
*%prepFCdata(dn_a_1_fc1,any_dnreg_1h);
*%prepFCdata(dn_a_72_fc1,any_dnreg_72h);

data counts_for_model;
  set counts4model_: ;
run;

/*transpose for repeated measures */

proc sort data=counts_for_model;
  by gene_set gene_id comparison hour;
proc transpose data=counts_for_model out=counts_for_model_RM;
  by gene_set gene_id comparison;
  id hour;
  var FC;
run;

proc sgplot data=counts_for_model;
  by gene_set;
  histogram log2FC ;
  density log2FC ;
run;


proc sort data=counts_for_model;
   by gene_set hour comparison;
run;
proc sort data=counts_for_model_RM;
   by gene_set comparison;
run;



proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison|hour / ddfm=kr;
ods output tests3=model1 ;
run;
quit;

proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison hour;
  model log2FC = comparison|hour / ddfm=kr;
ods output tests3=model2 ;
run;
quit;


proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison / ddfm=kr;
  random hour;
ods output tests3=model3 ;
run;
quit;

proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison  / solution;
  random intercept / subject=hour;
ods output tests3=model4 ;
run;



proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison  / solution;
  random _residual_ / subject=hour;
ods output tests3=model4a ;
run;


proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison hour / solution;
  random intercept / subject=hour;
ods output tests3=model5 ;
run;

proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison / solution;
  random intercept / subject=hour;
  random intercept / subject=hour(comparison);
ods output tests3=model6 ;
run;

proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison hour / solution;
  random intercept / subject=hour;
  random intercept / subject=hour(comparison);
ods output tests3=model7 ;
run;

proc glimmix data=counts_for_model ;
  by gene_set;
  class comparison;
  model log2FC = comparison|hour / solution;
  random intercept / subject=hour;
  random intercept / subject=hour(comparison);
ods output tests3=model8 ;
run;






