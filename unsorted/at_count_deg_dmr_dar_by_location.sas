/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

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


data  up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
      dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1
      no_01_1_fc1 no_01_3_fc1 no_01_24_fc1 no_01_72_fc1
      up_1_1_fc1 up_1_3_fc1 up_1_24_fc1 up_1_72_fc1
      dn_1_1_fc1 dn_1_3_fc1 dn_1_24_fc1 dn_1_72_fc1
      no_1_1_fc1 no_1_3_fc1 no_1_24_fc1 no_1_72_fc1
      off_01_1_fc1 off_01_3_fc1 off_01_24_fc1 off_01_72_fc1
      off_1_1_fc1 off_1_3_fc1 off_1_24_fc1 off_1_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
else if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
else if (mean_cpm_01gy_1h > 0 and mean_cpm_Mock_1h > 0) then output no_01_1_fc1;
else output off_01_1_fc1;

if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3_fc1;
else if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3_fc1;
else if (mean_cpm_01gy_3h > 0 and mean_cpm_Mock_3h > 0) then output no_01_3_fc1;
else output off_01_3_fc1;

if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24_fc1;
else if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24_fc1;
else if (mean_cpm_01gy_24h > 0 and mean_cpm_Mock_24h > 0) then output no_01_24_fc1;
else output off_01_24_fc1;

if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
else if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;
else if (mean_cpm_01gy_72h > 0 and mean_cpm_Mock_72h > 0) then output no_01_72_fc1;
else output off_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
else if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
else if (mean_cpm_1gy_1h > 0 and mean_cpm_Mock_1h > 0) then output no_1_1_fc1;
else output off_1_1_fc1;

if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3_fc1;
else if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3_fc1;
else if (mean_cpm_1gy_3h > 0 and mean_cpm_Mock_3h > 0) then output no_1_3_fc1;
else output off_1_3_fc1;

if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24_fc1;
else if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24_fc1;
else if (mean_cpm_1gy_24h > 0 and mean_cpm_Mock_24h > 0) then output no_1_24_fc1;
else output off_1_24_fc1;

if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
else if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;
else if (mean_cpm_1gy_72h > 0 and mean_cpm_Mock_72h > 0) then output no_1_72_fc1;
else output off_1_72_fc1;

keep gene_id;
rename gene_id=geneID;
run;


data any_up_01;   set up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1; run;
data any_dn_01;   set dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1; run;
data any_no_01;   set no_01_1_fc1 no_01_3_fc1 no_01_24_fc1 no_01_72_fc1; run;
data any_up_1;  set  up_1_1_fc1 up_1_3_fc1 up_1_24_fc1 up_1_72_fc1 ; run;
data any_dn_1;  set dn_1_1_fc1 dn_1_3_fc1 dn_1_24_fc1 dn_1_72_fc1 ; run;
data any_no_1;  set no_1_1_fc1 no_1_3_fc1 no_1_24_fc1 no_1_72_fc1 ; run;

data any_off_01;   set off_01_1_fc1 off_01_3_fc1 off_01_24_fc1 off_01_72_fc1; run;
data any_off_1;  set off_1_1_fc1 off_1_3_fc1 off_1_24_fc1 off_1_72_fc1 ; run;



proc sort data=any_up_01 nodup;
by geneID;
proc sort data=any_dn_01 nodup;
by geneID;
proc sort data=any_no_01 nodup;
by geneID;
proc sort data=any_up_1 nodup;
by geneID;
proc sort data=any_dn_1 nodup;
by geneID;
proc sort data=any_no_1 nodup;
by geneID;
proc sort data=any_off_01 nodup;
by geneID;
proc sort data=any_off_1 nodup;
by geneID;

run;


data any_no_01a;
  merge any_no_01 (in=in1) any_up_01 (in=in2) any_dn_01 (in=in3);
  by geneID;
  if in2 or in3 then delete;
run;

data any_no_1a;
  merge any_no_1 (in=in1) any_up_1 (in=in2) any_dn_1 (in=in3);
  by geneID;
  if in2 or in3 then delete;
run;



data any_off_01a;
  merge any_off_01 (in=in1) any_up_01 (in=in2) any_dn_01 (in=in3) any_no_01 (in=in4);
  by geneID;
  if in2 or in3 or in4 then delete;
run;

data any_off_1a;
  merge any_off_1 (in=in1) any_up_1 (in=in2) any_dn_1 (in=in3) any_no_1 (in=in4);
  by geneID;
  if in2 or in3 or in4 then delete;
run;


/* Get DMRs/DARs by gene  and feature annotation */



data up_dar_01_72_gn up_dar_1_72_gn dn_dar_01_72_gn dn_dar_1_72_gn;
     set arabMAP.results_by_dar_annot;
     if geneID = "" or annotation="Intergenic" then delete;
     length feature $12.;
     feature = compress(scan(annotation, 1, " "));
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

     keep geneID feature;
run;





data cg_up_dmr_01_72_gn cg_up_dmr_1_72_gn cg_dn_dmr_01_72_gn cg_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CG";
     length feature $12.;
     feature = compress(scan(annotation, 1, " "));
     if geneID = "" or annotation="Intergenic" then delete;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72_gn;
     keep geneID feature;
run;



data chg_up_dmr_01_72_gn chg_up_dmr_1_72_gn chg_dn_dmr_01_72_gn chg_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHG";
     length feature $12.;
     feature = compress(scan(annotation, 1, " "));
     if geneID = "" or annotation="Intergenic" then delete;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72_gn;
     keep geneID feature;
run;


data chh_up_dmr_01_72_gn chh_up_dmr_1_72_gn chh_dn_dmr_01_72_gn chh_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHH";
     length feature $12.;
     feature = compress(scan(annotation, 1, " "));
     if geneID = "" or annotation="Intergenic" then delete;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_gn;
     keep geneID feature;
run;


proc sort data=up_dar_01_72_gn nodup; by geneID feature; run;
proc sort data=up_dar_1_72_gn nodup; by geneID feature; run;
proc sort data=dn_dar_01_72_gn nodup; by geneID feature; run;
proc sort data=dn_dar_1_72_gn nodup; by geneID feature; run;
proc sort data=cg_up_dmr_01_72_gn nodup; by geneID feature; run;
proc sort data=cg_up_dmr_1_72_gn nodup; by geneID feature; run;
proc sort data=cg_dn_dmr_01_72_gn nodup; by geneID feature; run;
proc sort data=cg_dn_dmr_1_72_gn nodup; by geneID feature; run;
proc sort data=chg_up_dmr_01_72_gn nodup; by geneID feature; run;
proc sort data=chg_up_dmr_1_72_gn nodup; by geneID feature; run;
proc sort data=chg_dn_dmr_01_72_gn nodup; by geneID feature; run;
proc sort data=chg_dn_dmr_1_72_gn nodup; by geneID feature; run;
proc sort data=chh_up_dmr_01_72_gn nodup; by geneID feature; run;
proc sort data=chh_up_dmr_1_72_gn nodup; by geneID feature; run;
proc sort data=chh_dn_dmr_01_72_gn nodup; by geneID feature; run;
proc sort data=chh_dn_dmr_1_72_gn nodup; by geneID feature; run;

proc sort data=any_up_01 nodup; by geneID ; run;
proc sort data=up_01_72_fc1 nodup; by geneID ; run;
proc sort data=up_01_24_fc1 nodup; by geneID ; run;
proc sort data=up_01_3_fc1 nodup; by geneID ; run;
proc sort data=up_01_1_fc1 nodup; by geneID ; run;

proc sort data=any_dn_01 nodup; by geneID ; run;
proc sort data=dn_01_72_fc1 nodup; by geneID ; run;
proc sort data=dn_01_24_fc1 nodup; by geneID ; run;
proc sort data=dn_01_3_fc1 nodup; by geneID ; run;
proc sort data=dn_01_1_fc1 nodup; by geneID ; run;

proc sort data=any_no_01a nodup; by geneID ; run;
proc sort data=no_01_72_fc1 nodup; by geneID ; run;
proc sort data=no_01_24_fc1 nodup; by geneID ; run;
proc sort data=no_01_3_fc1 nodup; by geneID ; run;
proc sort data=no_01_1_fc1 nodup; by geneID ; run;


proc sort data=any_up_1 nodup; by geneID ; run;
proc sort data=up_1_72_fc1 nodup; by geneID ; run;
proc sort data=up_1_24_fc1 nodup; by geneID ; run;
proc sort data=up_1_3_fc1 nodup; by geneID ; run;
proc sort data=up_1_1_fc1 nodup; by geneID ; run;

proc sort data=any_dn_1 nodup; by geneID ; run;
proc sort data=dn_1_72_fc1 nodup; by geneID ; run;
proc sort data=dn_1_24_fc1 nodup; by geneID ; run;
proc sort data=dn_1_3_fc1 nodup; by geneID ; run;
proc sort data=dn_1_1_fc1 nodup; by geneID ; run;

proc sort data=any_no_1a nodup; by geneID; run;
proc sort data=no_1_72_fc1 nodup; by geneID; run;
proc sort data=no_1_24_fc1 nodup; by geneID; run;
proc sort data=no_1_3_fc1 nodup; by geneID ; run;
proc sort data=no_1_1_fc1 nodup; by geneID ; run;







%macro mergeAndCount(inRNA, inMETH);

proc sort data=&inRNA.;
   by geneID;
proc sort data=&inMETH.;
   by geneID;
run;


data check;
    merge &inRNA. (in=in1) &inMETH. (in=in2);
    by geneID;
    if not in2 then feature="none";
    if in1;
run;

proc freq data=check;
   tables feature;
run;

%mend;
      
%mergeAndCount(up_01_1_fc1, up_dar_01_72_gn); 
%mergeAndCount(up_01_3_fc1, up_dar_01_72_gn); 
%mergeAndCount(up_01_24_fc1, up_dar_01_72_gn); 
%mergeAndCount(up_01_72_fc1, up_dar_01_72_gn); 
%mergeAndCount(any_up_01, up_dar_01_72_gn); 
      
%mergeAndCount(up_01_1_fc1, dn_dar_01_72_gn); 
%mergeAndCount(up_01_3_fc1, dn_dar_01_72_gn); 
%mergeAndCount(up_01_24_fc1, dn_dar_01_72_gn); 
%mergeAndCount(up_01_72_fc1, dn_dar_01_72_gn); 
%mergeAndCount(any_up_01, dn_dar_01_72_gn); 

%mergeAndCount(up_01_1_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(up_01_3_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(up_01_24_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(up_01_72_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(any_up_01, cg_up_dmr_01_72_gn); 
      
%mergeAndCount(up_01_1_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_3_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_24_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_72_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(any_up_01, cg_dn_dmr_01_72_gn); 

%mergeAndCount(up_01_1_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(up_01_3_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(up_01_24_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(up_01_72_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(any_up_01, chg_up_dmr_01_72_gn); 
      
%mergeAndCount(up_01_1_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_3_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_24_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_72_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(any_up_01, chg_dn_dmr_01_72_gn); 

%mergeAndCount(up_01_1_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(up_01_3_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(up_01_24_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(up_01_72_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(any_up_01, chh_up_dmr_01_72_gn); 
      
%mergeAndCount(up_01_1_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_3_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_24_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(up_01_72_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(any_up_01, chh_dn_dmr_01_72_gn); 

     
      
%mergeAndCount(dn_01_1_fc1, up_dar_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, up_dar_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, up_dar_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, up_dar_01_72_gn); 
%mergeAndCount(any_dn_01, up_dar_01_72_gn); 
      
%mergeAndCount(dn_01_1_fc1, dn_dar_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, dn_dar_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, dn_dar_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, dn_dar_01_72_gn); 
%mergeAndCount(any_dn_01, dn_dar_01_72_gn); 

%mergeAndCount(dn_01_1_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(any_dn_01, cg_up_dmr_01_72_gn); 
      
%mergeAndCount(dn_01_1_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(any_dn_01, cg_dn_dmr_01_72_gn); 

%mergeAndCount(dn_01_1_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(any_dn_01, chg_up_dmr_01_72_gn); 
      
%mergeAndCount(dn_01_1_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(any_dn_01, chg_dn_dmr_01_72_gn); 

%mergeAndCount(dn_01_1_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(any_dn_01, chh_up_dmr_01_72_gn); 
      
%mergeAndCount(dn_01_1_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_3_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_24_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(dn_01_72_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(any_dn_01, chh_dn_dmr_01_72_gn); 

     

%mergeAndCount(no_01_1_fc1, up_dar_01_72_gn); 
%mergeAndCount(no_01_3_fc1, up_dar_01_72_gn); 
%mergeAndCount(no_01_24_fc1, up_dar_01_72_gn); 
%mergeAndCount(no_01_72_fc1, up_dar_01_72_gn); 
%mergeAndCount(any_no_01a, up_dar_01_72_gn); 
      
%mergeAndCount(no_01_1_fc1, dn_dar_01_72_gn); 
%mergeAndCount(no_01_3_fc1, dn_dar_01_72_gn); 
%mergeAndCount(no_01_24_fc1, dn_dar_01_72_gn); 
%mergeAndCount(no_01_72_fc1, dn_dar_01_72_gn); 
%mergeAndCount(any_no_01a, dn_dar_01_72_gn); 

%mergeAndCount(no_01_1_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(no_01_3_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(no_01_24_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(no_01_72_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(any_no_01a, cg_up_dmr_01_72_gn); 
      
%mergeAndCount(no_01_1_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_3_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_24_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_72_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(any_no_01a, cg_dn_dmr_01_72_gn); 

%mergeAndCount(no_01_1_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(no_01_3_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(no_01_24_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(no_01_72_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(any_no_01a, chg_up_dmr_01_72_gn); 
      
%mergeAndCount(no_01_1_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_3_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_24_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_72_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(any_no_01a, chg_dn_dmr_01_72_gn); 

%mergeAndCount(no_01_1_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(no_01_3_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(no_01_24_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(no_01_72_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(any_no_01a, chh_up_dmr_01_72_gn); 
      
%mergeAndCount(no_01_1_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_3_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_24_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(no_01_72_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(any_no_01a, chh_dn_dmr_01_72_gn); 

     



      
%mergeAndCount(up_1_1_fc1, up_dar_1_72_gn); 
%mergeAndCount(up_1_3_fc1, up_dar_1_72_gn); 
%mergeAndCount(up_1_24_fc1, up_dar_1_72_gn); 
%mergeAndCount(up_1_72_fc1, up_dar_1_72_gn); 
%mergeAndCount(any_up_1, up_dar_1_72_gn); 
      
%mergeAndCount(up_1_1_fc1, dn_dar_1_72_gn); 
%mergeAndCount(up_1_3_fc1, dn_dar_1_72_gn); 
%mergeAndCount(up_1_24_fc1, dn_dar_1_72_gn); 
%mergeAndCount(up_1_72_fc1, dn_dar_1_72_gn); 
%mergeAndCount(any_up_1, dn_dar_1_72_gn); 

%mergeAndCount(up_1_1_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(up_1_3_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(up_1_24_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(up_1_72_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(any_up_1, cg_up_dmr_1_72_gn); 
      
%mergeAndCount(up_1_1_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_3_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_24_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_72_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(any_up_1, cg_dn_dmr_1_72_gn); 

%mergeAndCount(up_1_1_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(up_1_3_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(up_1_24_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(up_1_72_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(any_up_1, chg_up_dmr_1_72_gn); 
      
%mergeAndCount(up_1_1_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_3_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_24_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_72_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(any_up_1, chg_dn_dmr_1_72_gn); 

%mergeAndCount(up_1_1_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(up_1_3_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(up_1_24_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(up_1_72_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(any_up_1, chh_up_dmr_1_72_gn); 
      
%mergeAndCount(up_1_1_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_3_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_24_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(up_1_72_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(any_up_1, chh_dn_dmr_1_72_gn); 

     
      
%mergeAndCount(dn_1_1_fc1, up_dar_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, up_dar_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, up_dar_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, up_dar_1_72_gn); 
%mergeAndCount(any_dn_1, up_dar_1_72_gn); 
      
%mergeAndCount(dn_1_1_fc1, dn_dar_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, dn_dar_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, dn_dar_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, dn_dar_1_72_gn); 
%mergeAndCount(any_dn_1, dn_dar_1_72_gn); 

%mergeAndCount(dn_1_1_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(any_dn_1, cg_up_dmr_1_72_gn); 
      
%mergeAndCount(dn_1_1_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(any_dn_1, cg_dn_dmr_1_72_gn); 

%mergeAndCount(dn_1_1_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(any_dn_1, chg_up_dmr_1_72_gn); 
      
%mergeAndCount(dn_1_1_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(any_dn_1, chg_dn_dmr_1_72_gn); 

%mergeAndCount(dn_1_1_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(any_dn_1, chh_up_dmr_1_72_gn); 
      
%mergeAndCount(dn_1_1_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_3_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_24_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(dn_1_72_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(any_dn_1, chh_dn_dmr_1_72_gn); 

     

%mergeAndCount(no_1_1_fc1, up_dar_1_72_gn); 
%mergeAndCount(no_1_3_fc1, up_dar_1_72_gn); 
%mergeAndCount(no_1_24_fc1, up_dar_1_72_gn); 
%mergeAndCount(no_1_72_fc1, up_dar_1_72_gn); 
%mergeAndCount(any_no_1a, up_dar_1_72_gn); 
      
%mergeAndCount(no_1_1_fc1, dn_dar_1_72_gn); 
%mergeAndCount(no_1_3_fc1, dn_dar_1_72_gn); 
%mergeAndCount(no_1_24_fc1, dn_dar_1_72_gn); 
%mergeAndCount(no_1_72_fc1, dn_dar_1_72_gn); 
%mergeAndCount(any_no_1a, dn_dar_1_72_gn); 

%mergeAndCount(no_1_1_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(no_1_3_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(no_1_24_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(no_1_72_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(any_no_1a, cg_up_dmr_1_72_gn); 
      
%mergeAndCount(no_1_1_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_3_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_24_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_72_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(any_no_1a, cg_dn_dmr_1_72_gn); 

%mergeAndCount(no_1_1_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(no_1_3_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(no_1_24_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(no_1_72_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(any_no_1a, chg_up_dmr_1_72_gn); 
      
%mergeAndCount(no_1_1_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_3_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_24_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_72_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(any_no_1a, chg_dn_dmr_1_72_gn); 

%mergeAndCount(no_1_1_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(no_1_3_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(no_1_24_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(no_1_72_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(any_no_1a, chh_up_dmr_1_72_gn); 
      
%mergeAndCount(no_1_1_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_3_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_24_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(no_1_72_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(any_no_1a, chh_dn_dmr_1_72_gn); 

     









%mergeAndCount(off_01_1_fc1, up_dar_01_72_gn); 
%mergeAndCount(off_01_3_fc1, up_dar_01_72_gn); 
%mergeAndCount(off_01_24_fc1, up_dar_01_72_gn); 
%mergeAndCount(off_01_72_fc1, up_dar_01_72_gn); 
%mergeAndCount(any_off_01a, up_dar_01_72_gn); 
      
%mergeAndCount(off_01_1_fc1, dn_dar_01_72_gn); 
%mergeAndCount(off_01_3_fc1, dn_dar_01_72_gn); 
%mergeAndCount(off_01_24_fc1, dn_dar_01_72_gn); 
%mergeAndCount(off_01_72_fc1, dn_dar_01_72_gn); 
%mergeAndCount(any_off_01a, dn_dar_01_72_gn); 

%mergeAndCount(off_01_1_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(off_01_3_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(off_01_24_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(off_01_72_fc1, cg_up_dmr_01_72_gn); 
%mergeAndCount(any_off_01a, cg_up_dmr_01_72_gn); 
      
%mergeAndCount(off_01_1_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_3_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_24_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_72_fc1, cg_dn_dmr_01_72_gn); 
%mergeAndCount(any_off_01a, cg_dn_dmr_01_72_gn); 

%mergeAndCount(off_01_1_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(off_01_3_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(off_01_24_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(off_01_72_fc1, chg_up_dmr_01_72_gn); 
%mergeAndCount(any_off_01a, chg_up_dmr_01_72_gn); 
      
%mergeAndCount(off_01_1_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_3_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_24_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_72_fc1, chg_dn_dmr_01_72_gn); 
%mergeAndCount(any_off_01a, chg_dn_dmr_01_72_gn); 

%mergeAndCount(off_01_1_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(off_01_3_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(off_01_24_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(off_01_72_fc1, chh_up_dmr_01_72_gn); 
%mergeAndCount(any_off_01a, chh_up_dmr_01_72_gn); 
      
%mergeAndCount(off_01_1_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_3_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_24_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(off_01_72_fc1, chh_dn_dmr_01_72_gn); 
%mergeAndCount(any_off_01a, chh_dn_dmr_01_72_gn); 

     



%mergeAndCount(off_1_1_fc1, up_dar_1_72_gn); 
%mergeAndCount(off_1_3_fc1, up_dar_1_72_gn); 
%mergeAndCount(off_1_24_fc1, up_dar_1_72_gn); 
%mergeAndCount(off_1_72_fc1, up_dar_1_72_gn); 
%mergeAndCount(any_off_1a, up_dar_1_72_gn); 
      
%mergeAndCount(off_1_1_fc1, dn_dar_1_72_gn); 
%mergeAndCount(off_1_3_fc1, dn_dar_1_72_gn); 
%mergeAndCount(off_1_24_fc1, dn_dar_1_72_gn); 
%mergeAndCount(off_1_72_fc1, dn_dar_1_72_gn); 
%mergeAndCount(any_off_1a, dn_dar_1_72_gn); 

%mergeAndCount(off_1_1_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(off_1_3_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(off_1_24_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(off_1_72_fc1, cg_up_dmr_1_72_gn); 
%mergeAndCount(any_off_1a, cg_up_dmr_1_72_gn); 
      
%mergeAndCount(off_1_1_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_3_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_24_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_72_fc1, cg_dn_dmr_1_72_gn); 
%mergeAndCount(any_off_1a, cg_dn_dmr_1_72_gn); 

%mergeAndCount(off_1_1_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(off_1_3_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(off_1_24_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(off_1_72_fc1, chg_up_dmr_1_72_gn); 
%mergeAndCount(any_off_1a, chg_up_dmr_1_72_gn); 
      
%mergeAndCount(off_1_1_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_3_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_24_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_72_fc1, chg_dn_dmr_1_72_gn); 
%mergeAndCount(any_off_1a, chg_dn_dmr_1_72_gn); 

%mergeAndCount(off_1_1_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(off_1_3_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(off_1_24_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(off_1_72_fc1, chh_up_dmr_1_72_gn); 
%mergeAndCount(any_off_1a, chh_up_dmr_1_72_gn); 
      
%mergeAndCount(off_1_1_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_3_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_24_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(off_1_72_fc1, chh_dn_dmr_1_72_gn); 
%mergeAndCount(any_off_1a, chh_dn_dmr_1_72_gn); 

     

