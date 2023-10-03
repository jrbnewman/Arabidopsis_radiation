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


data up_01_any_fc1 up_01_72_fc1
     dn_01_any_fc1 dn_01_72_fc1
     up_1_any_fc1  up_1_72_fc1
     dn_1_any_fc1  dn_1_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_any_fc1;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_any_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_any_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_any_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_any_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_any_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_any_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_any_fc1;



if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_any_fc1;
if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_any_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_any_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_any_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_any_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_any_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_any_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_any_fc1;
keep gene_id;
rename gene_id=geneID;
run;


data access_genes;
  set arabMAP.results_by_dac_annot;
  if annotation="Intergenic" or geneID="" then delete;
  keep geneID;
run;

data cg_meth_genes chg_meth_genes chh_meth_genes ;
  set arabMAP.results_by_dmc_annot;
  if annotation="Intergenic" or geneID="" then delete;
  if site_type="CG" then output cg_meth_genes;
  if site_type="CHG" then output chg_meth_genes;
  if site_type="CHH" then output chh_meth_genes;
  keep geneID;
run;

proc sort data=access_genes nodup;
  by geneID;
proc sort data=cg_meth_genes nodup;
  by geneID;
proc sort data=chg_meth_genes nodup;
  by geneID;
proc sort data=chh_meth_genes nodup;
  by geneID;
run;




data up_dar_01_72_gn up_dar_1_72_gn dn_dar_01_72_gn dn_dar_1_72_gn;
     set arabMAP.results_by_dar_annot;
     if geneID = "" or annotation="Intergenic" then delete;
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





data cg_up_dmr_01_72_gn cg_up_dmr_1_72_gn cg_dn_dmr_01_72_gn cg_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CG";
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
     keep geneID;
run;



data chg_up_dmr_01_72_gn chg_up_dmr_1_72_gn chg_dn_dmr_01_72_gn chg_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHG";
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
     keep geneID;
run;


data chh_up_dmr_01_72_gn chh_up_dmr_1_72_gn chh_dn_dmr_01_72_gn chh_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHH";
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
     keep geneID;
run;


proc sort data=up_dar_01_72_gn nodup; by geneID; run;
proc sort data=up_dar_1_72_gn nodup; by geneID; run;
proc sort data=dn_dar_01_72_gn nodup; by geneID; run;
proc sort data=dn_dar_1_72_gn nodup; by geneID; run;
proc sort data=cg_up_dmr_01_72_gn nodup; by geneID; run;
proc sort data=cg_up_dmr_1_72_gn nodup; by geneID; run;
proc sort data=cg_dn_dmr_01_72_gn nodup; by geneID; run;
proc sort data=cg_dn_dmr_1_72_gn nodup; by geneID; run;
proc sort data=chg_up_dmr_01_72_gn nodup; by geneID; run;
proc sort data=chg_up_dmr_1_72_gn nodup; by geneID; run;
proc sort data=chg_dn_dmr_01_72_gn nodup; by geneID; run;
proc sort data=chg_dn_dmr_1_72_gn nodup; by geneID; run;
proc sort data=chh_up_dmr_01_72_gn nodup; by geneID; run;
proc sort data=chh_up_dmr_1_72_gn nodup; by geneID; run;
proc sort data=chh_dn_dmr_01_72_gn nodup; by geneID; run;
proc sort data=chh_dn_dmr_1_72_gn nodup; by geneID; run;
proc sort data=up_01_any_fc1 nodup; by geneID; run;
proc sort data=up_01_72_fc1 nodup; by geneID; run;
proc sort data=dn_01_any_fc1 nodup; by geneID; run;
proc sort data=dn_01_72_fc1 nodup; by geneID; run;
proc sort data=up_1_any_fc1 nodup; by geneID; run;
proc sort data=up_1_72_fc1 nodup; by geneID; run;
proc sort data=dn_1_any_fc1 nodup; by geneID; run;
proc sort data=dn_1_72_fc1 nodup; by geneID; run;


data sig_table;
  merge up_01_any_fc1 (in=in1) dn_01_any_fc1 (in=in2) up_01_72_fc1 (in=in3) dn_01_72_fc1 (in=in4)
        up_1_any_fc1 (in=in5) dn_1_any_fc1 (in=in6) up_1_72_fc1 (in=in7) dn_1_72_fc1 (in=in8)
        up_dar_01_72_gn (in=in9) dn_dar_01_72_gn (in=in10) up_dar_1_72_gn (in=in11) dn_dar_1_72_gn (in=in12)
        cg_up_dmr_01_72_gn (in=in13) cg_dn_dmr_01_72_gn (in=in14) cg_up_dmr_1_72_gn (in=in15) cg_dn_dmr_1_72_gn (in=in16)
        chg_up_dmr_01_72_gn (in=in17) chg_dn_dmr_01_72_gn (in=in18) chg_up_dmr_1_72_gn (in=in19) chg_dn_dmr_1_72_gn (in=in20)
        chh_up_dmr_01_72_gn (in=in21) chh_dn_dmr_01_72_gn (in=in22) chh_up_dmr_1_72_gn (in=in23) chh_dn_dmr_1_72_gn (in=in24);
   by geneID;
   if in1 then flag_up_01_deg_any=1; else flag_up_01_deg_any=0;
   if in2 then flag_dn_01_deg_any=1; else flag_dn_01_deg_any=0;
   if in3 then flag_up_01_deg_72=1; else flag_up_01_deg_72=0;
   if in4 then flag_dn_01_deg_72=1; else flag_dn_01_deg_72=0;
   if in5 then flag_up_1_deg_any=1; else flag_up_1_deg_any=0;
   if in6 then flag_dn_1_deg_any=1; else flag_dn_1_deg_any=0;
   if in7 then flag_up_1_deg_72=1; else flag_up_1_deg_72=0;
   if in8 then flag_dn_1_deg_72=1; else flag_dn_1_deg_72=0;

   if in9 then flag_up_dar_01=1; else flag_up_dar_01=0;
   if in10 then flag_dn_dar_01=1; else flag_dn_dar_01=0;
   if in11 then flag_up_dar_1=1; else flag_up_dar_1=0;
   if in12 then flag_dn_dar_1=1; else flag_dn_dar_1=0;

   if in13 then flag_cg_up_dmr_01=1; else flag_cg_up_dmr_01=0;
   if in14 then flag_cg_dn_dmr_01=1; else flag_cg_dn_dmr_01=0;
   if in15 then flag_cg_up_dmr_1=1; else flag_cg_up_dmr_1=0;
   if in16 then flag_cg_dn_dmr_1=1; else flag_cg_dn_dmr_1=0;

   if in17 then flag_chg_up_dmr_01=1; else flag_chg_up_dmr_01=0;
   if in18 then flag_chg_dn_dmr_01=1; else flag_chg_dn_dmr_01=0;
   if in19 then flag_chg_up_dmr_1=1; else flag_chg_up_dmr_1=0;
   if in20 then flag_chg_dn_dmr_1=1; else flag_chg_dn_dmr_1=0;

   if in21 then flag_chh_up_dmr_01=1; else flag_chh_up_dmr_01=0;
   if in22 then flag_chh_dn_dmr_01=1; else flag_chh_dn_dmr_01=0;
   if in23 then flag_chh_up_dmr_1=1; else flag_chh_up_dmr_1=0;
   if in24 then flag_chh_dn_dmr_1=1; else flag_chh_dn_dmr_1=0;
run;




libname at "!HOME/concannon/useful_arabidopsis_data/TAIR10/sas_data";




data at2go;
   set at.gene2go;
   keep gene_ID GO_number_biopro_cat;
   rename gene_id=geneID;
run;


proc sort data=at2go;
  by geneID;
proc sort data=sig_table;
  by geneID;
run;

data sig_table_w_go;
  merge sig_table (in=in1) at2go (in=in2);
   by geneID;
run;

data sig_table_W_go2;
  set sig_table_w_go;
  array change _numeric_;
  do over change;
  if change=. then change=0;
  end;
run;



proc freq data=sig_table_w_go2 noprint;
  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_up_dar_01*flag_dn_dar_01 / out=deg_dar_01_any;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_up_dar_1*flag_dn_dar_1 / out=deg_dar_1_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_up_dar_01*flag_dn_dar_01 / out=deg_dar_01_72;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_up_dar_1*flag_dn_dar_1 / out=deg_dar_1_72;

  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_cg_up_dmr_01*flag_cg_dn_dmr_01 / out=deg_dmr_cg_01_any;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_cg_up_dmr_1*flag_cg_dn_dmr_1 / out=deg_dmr_cg_1_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_cg_up_dmr_01*flag_cg_dn_dmr_01 / out=deg_dmr_cg_01_72;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_cg_up_dmr_1*flag_cg_dn_dmr_1 / out=deg_dmr_cg_1_72;

  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_chg_up_dmr_01*flag_chg_dn_dmr_01 / out=deg_dmr_chg_01_any;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_chg_up_dmr_1*flag_chg_dn_dmr_1 / out=deg_dmr_chg_1_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_chg_up_dmr_01*flag_chg_dn_dmr_01 / out=deg_dmr_chg_01_72;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_chg_up_dmr_1*flag_chg_dn_dmr_1 / out=deg_dmr_chg_1_72;

  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_chh_up_dmr_01*flag_chh_dn_dmr_01 / out=deg_dmr_chh_01_any;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_chh_up_dmr_1*flag_chh_dn_dmr_1 / out=deg_dmr_chh_1_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_chh_up_dmr_01*flag_chh_dn_dmr_01 / out=deg_dmr_chh_01_72;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_chh_up_dmr_1*flag_chh_dn_dmr_1 / out=deg_dmr_chh_1_72;

  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_up_dar_01*flag_dn_dar_01*flag_cg_up_dmr_01*flag_cg_dn_dmr_01 / out=deg_dar_dmr_cg_01_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_up_dar_01*flag_dn_dar_01*flag_cg_up_dmr_01*flag_cg_dn_dmr_01 / out=deg_dar_dmr_cg_01_72;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_up_dar_1*flag_dn_dar_1*flag_cg_up_dmr_1*flag_cg_dn_dmr_1 / out=deg_dar_dmr_cg_1_any;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_up_dar_1*flag_dn_dar_1*flag_cg_up_dmr_1*flag_cg_dn_dmr_1 / out=deg_dar_dmr_cg_1_72;

  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_up_dar_01*flag_dn_dar_01*flag_chg_up_dmr_01*flag_chg_dn_dmr_01 / out=deg_dar_dmr_chg_01_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_up_dar_01*flag_dn_dar_01*flag_chg_up_dmr_01*flag_chg_dn_dmr_01 / out=deg_dar_dmr_chg_01_72;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_up_dar_1*flag_dn_dar_1*flag_chg_up_dmr_1*flag_chg_dn_dmr_1 / out=deg_dar_dmr_chg_1_any;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_up_dar_1*flag_dn_dar_1*flag_chg_up_dmr_1*flag_chg_dn_dmr_1 / out=deg_dar_dmr_chg_1_72;

  tables flag_up_01_deg_any*flag_dn_01_deg_any*flag_up_dar_01*flag_dn_dar_01*flag_chh_up_dmr_01*flag_chh_dn_dmr_01 / out=deg_dar_dmr_chh_01_any;
  tables flag_up_01_deg_72*flag_dn_01_deg_72*flag_up_dar_01*flag_dn_dar_01*flag_chh_up_dmr_01*flag_chh_dn_dmr_01 / out=deg_dar_dmr_chh_01_72;
  tables flag_up_1_deg_any*flag_dn_1_deg_any*flag_up_dar_1*flag_dn_dar_1*flag_chh_up_dmr_1*flag_chh_dn_dmr_1 / out=deg_dar_dmr_chh_1_any;
  tables flag_up_1_deg_72*flag_dn_1_deg_72*flag_up_dar_1*flag_dn_dar_1*flag_chh_up_dmr_1*flag_chh_dn_dmr_1 / out=deg_dar_dmr_chh_1_72;
run;

proc print data=deg_dar_01_any;
proc print data=deg_dar_1_any;
proc print data=deg_dar_01_72;
proc print data=deg_dar_1_72;
proc print data=deg_dmr_cg_01_any;
proc print data=deg_dmr_cg_1_any;
proc print data=deg_dmr_cg_01_72;
proc print data=deg_dmr_cg_1_72;
proc print data=deg_dmr_chg_01_any;
proc print data=deg_dmr_chg_1_any;
proc print data=deg_dmr_chg_01_72;
proc print data=deg_dmr_chg_1_72;
proc print data=deg_dmr_chh_01_any;
proc print data=deg_dmr_chh_1_any;
proc print data=deg_dmr_chh_01_72;
proc print data=deg_dmr_chh_1_72;
proc print data=deg_dar_dmr_cg_01_any;
proc print data=deg_dar_dmr_cg_01_72;
proc print data=deg_dar_dmr_cg_1_any;
proc print data=deg_dar_dmr_cg_1_72;
proc print data=deg_dar_dmr_chg_01_any;
proc print data=deg_dar_dmr_chg_01_72;
proc print data=deg_dar_dmr_chg_1_any;
proc print data=deg_dar_dmr_chg_1_72;
proc print data=deg_dar_dmr_chh_01_any;
proc print data=deg_dar_dmr_chh_01_72;
proc print data=deg_dar_dmr_chh_1_any;
proc print data=deg_dar_dmr_chh_1_72;
run;

/*

  flag_up_    flag_dn_
   01_deg_     01_deg_    flag_up_    flag_dn_
     any         any       dar_01      dar_01     COUNT

      0           0           0           0       14219
      0           0           0           1          52
      0           0           1           0       15928
      0           0           1           1         136
      0           1           0           0          77
      0           1           0           1           1
      0           1           1           0          58
      0           1           1           1           1
      1           0           0           0         240
      1           0           0           1           3
      1           0           1           0         219
      1           1           1           0           2


   flag_up_    flag_dn_
    1_deg_      1_deg_     flag_up_    flag_dn_
      any         any        dar_1       dar_1     COUNT

       0           0           0           0       21021
       0           0           0           1        1748
       0           0           1           0        6791
       0           0           1           1         840
       0           1           0           0          80
       0           1           0           1           6
       0           1           1           0          17
       0           1           1           1           1
       1           0           0           0         332
       1           0           0           1          16
       1           0           1           0          75
       1           0           1           1           9




 flag_up_    flag_dn_
  01_deg_     01_deg_    flag_up_    flag_dn_
    72          72        dar_01      dar_01     COUNT

     0           0           0           0       14501
     0           0           0           1          56
     0           0           1           0       16176
     0           0           1           1         137
     0           1           0           0          12
     0           1           1           0           3
     1           0           0           0          23
     1           0           1           0          28



  flag_up_    flag_dn_    flag_up_    flag_dn_
  1_deg_72    1_deg_72      dar_1       dar_1     COUNT

      0           0           0           0       21402
      0           0           0           1        1768
      0           0           1           0        6873
      0           0           1           1         849
      0           1           0           0           4
      0           1           1           0           1
      1           0           0           0          27
      1           0           0           1           2
      1           0           1           0           9
      1           0           1           1           1


 flag_up_    flag_dn_    flag_cg_    flag_cg_
  01_deg_     01_deg_     up_dmr_     dn_dmr_
    any         any         01          01       COUNT

     0           0           0           0       30239
     0           0           0           1          75
     0           0           1           0          20
     0           0           1           1           1
     0           1           0           0         136
     0           1           0           1           1
     1           0           0           0         462
     1           1           0           0           2


  flag_up_    flag_dn_
   1_deg_      1_deg_     flag_cg_    flag_cg_
     any         any      up_dmr_1    dn_dmr_1    COUNT

      0           0           0           0       30380
      0           0           0           1          18
      0           0           1           0           2
      0           1           0           0         104
      1           0           0           0         432

 flag_up_    flag_dn_    flag_cg_    flag_cg_
  01_deg_     01_deg_     up_dmr_     dn_dmr_
    72          72          01          01       COUNT

     0           0           0           0       30774
     0           0           0           1          75
     0           0           1           0          20
     0           0           1           1           1
     0           1           0           0          14
     0           1           0           1           1
     1           0           0           0          51


   flag_up_    flag_dn_    flag_cg_    flag_cg_
   1_deg_72    1_deg_72    up_dmr_1    dn_dmr_1    COUNT

       0           0           0           0       30872
       0           0           0           1          18
       0           0           1           0           2
       0           1           0           0           5
       1           0           0           0          39


 flag_up_    flag_dn_     flag_      flag_
  01_deg_     01_deg_    chg_up_    chg_dn_
    any         any       dmr_01     dmr_01    COUNT

     0           0          0          0       30244
     0           0          0          1          64
     0           0          1          0          26
     0           0          1          1           1
     0           1          0          0         137
     1           0          0          0         462
     1           1          0          0           2


 flag_up_    flag_dn_     flag_      flag_
  1_deg_      1_deg_     chg_up_    chg_dn_
    any         any       dmr_1      dmr_1     COUNT

     0           0          0          0       30396
     0           0          0          1           1
     0           0          1          0           3
     0           1          0          0         104
     1           0          0          0         432


flag_up_    flag_dn_     flag_      flag_
 01_deg_     01_deg_    chg_up_    chg_dn_
   72          72        dmr_01     dmr_01    COUNT

    0           0          0          0       30779
    0           0          0          1          64
    0           0          1          0          26
    0           0          1          1           1
    0           1          0          0          15
    1           0          0          0          51


                          flag_      flag_
 flag_up_    flag_dn_    chg_up_    chg_dn_
 1_deg_72    1_deg_72     dmr_1      dmr_1     COUNT

     0           0          0          0       30888
     0           0          0          1           1
     0           0          1          0           3
     0           1          0          0           5
     1           0          0          0          39


 flag_up_    flag_dn_     flag_      flag_
  01_deg_     01_deg_    chh_up_    chh_dn_
    any         any       dmr_01     dmr_01    COUNT

     0           0          0          0       25291
     0           0          0          1        4679
     0           0          1          0         208
     0           0          1          1         157
     0           1          0          0         120
     0           1          0          1          17
     1           0          0          0         412
     1           0          0          1          49
     1           0          1          0           1
     1           1          0          1           2

 flag_up_    flag_dn_     flag_      flag_
  1_deg_      1_deg_     chh_up_    chh_dn_
    any         any       dmr_1      dmr_1     COUNT

     0           0          0          0       27348
     0           0          0          1         695
     0           0          1          0        1984
     0           0          1          1         373
     0           1          0          0          96
     0           1          0          1           2
     0           1          1          0           5
     0           1          1          1           1
     1           0          0          0         401
     1           0          0          1           3
     1           0          1          0          27
     1           0          1          1           1

  flag_up_    flag_dn_     flag_      flag_
   01_deg_     01_deg_    chh_up_    chh_dn_
     72          72        dmr_01     dmr_01    COUNT

      0           0          0          0       25764
      0           0          0          1        4740
      0           0          1          0         209
      0           0          1          1         157
      0           1          0          0          12
      0           1          0          1           3
      1           0          0          0          47
      1           0          0          1           4


                           flag_      flag_
  flag_up_    flag_dn_    chh_up_    chh_dn_
  1_deg_72    1_deg_72     dmr_1      dmr_1     COUNT

      0           0          0          0       27804
      0           0          0          1         699
      0           0          1          0        2014
      0           0          1          1         375
      0           1          0          0           4
      0           1          0          1           1
      1           0          0          0          37
      1           0          1          0           2


flag_up_    flag_dn_                            flag_cg_    flag_cg_
 01_deg_     01_deg_    flag_up_    flag_dn_     up_dmr_     dn_dmr_
   any         any       dar_01      dar_01        01          01       COUNT

    0           0           0           0           0           0       14173
    0           0           0           0           0           1          40
    0           0           0           0           1           0           6
    0           0           0           1           0           0          47
    0           0           0           1           0           1           1
    0           0           0           1           1           0           4
    0           0           1           0           0           0       15891
    0           0           1           0           0           1          32
    0           0           1           0           1           0           5
    0           0           1           1           0           0         128
    0           0           1           1           0           1           2
    0           0           1           1           1           0           5
    0           0           1           1           1           1           1
    0           1           0           0           0           0          76
    0           1           0           0           0           1           1
    0           1           0           1           0           0           1
    0           1           1           0           0           0          58
    0           1           1           1           0           0           1
    1           0           0           0           0           0         240
    1           0           0           1           0           0           3
    1           0           1           0           0           0         219
    1           1           1           0           0           0           2






 flag_up_    flag_dn_                            flag_cg_    flag_cg_
  01_deg_     01_deg_    flag_up_    flag_dn_     up_dmr_     dn_dmr_
    72          72        dar_01      dar_01        01          01       COUNT

     0           0           0           0           0           0       14455
     0           0           0           0           0           1          40
     0           0           0           0           1           0           6
     0           0           0           1           0           0          51
     0           0           0           1           0           1           1
     0           0           0           1           1           0           4
     0           0           1           0           0           0       16139
     0           0           1           0           0           1          32
     0           0           1           0           1           0           5
     0           0           1           1           0           0         129
     0           0           1           1           0           1           2
     0           0           1           1           1           0           5
     0           0           1           1           1           1           1
     0           1           0           0           0           0          11
     0           1           0           0           0           1           1
     0           1           1           0           0           0           3
     1           0           0           0           0           0          23
     1           0           1           0           0           0          28

                                   The SAS System

  flag_up_    flag_dn_
   1_deg_      1_deg_     flag_up_    flag_dn_    flag_cg_    flag_cg_
     any         any        dar_1       dar_1     up_dmr_1    dn_dmr_1    COUNT

      0           0           0           0           0           0       21007
      0           0           0           0           0           1          13
      0           0           0           0           1           0           1
      0           0           0           1           0           0        1746
      0           0           0           1           0           1           1
      0           0           0           1           1           0           1
      0           0           1           0           0           0        6787
      0           0           1           0           0           1           4
      0           0           1           1           0           0         840
      0           1           0           0           0           0          80
      0           1           0           1           0           0           6
      0           1           1           0           0           0          17
      0           1           1           1           0           0           1
      1           0           0           0           0           0         332
      1           0           0           1           0           0          16
      1           0           1           0           0           0          75
      1           0           1           1           0           0           9

 flag_up_    flag_dn_    flag_up_    flag_dn_    flag_cg_    flag_cg_
 1_deg_72    1_deg_72      dar_1       dar_1     up_dmr_1    dn_dmr_1    COUNT

     0           0           0           0           0           0       21388
     0           0           0           0           0           1          13
     0           0           0           0           1           0           1
     0           0           0           1           0           0        1766
     0           0           0           1           0           1           1
     0           0           0           1           1           0           1
     0           0           1           0           0           0        6869
     0           0           1           0           0           1           4
     0           0           1           1           0           0         849
     0           1           0           0           0           0           4
     0           1           1           0           0           0           1
     1           0           0           0           0           0          27
     1           0           0           1           0           0           2
     1           0           1           0           0           0           9
     1           0           1           1           0           0           1

flag_up_    flag_dn_                             flag_      flag_
 01_deg_     01_deg_    flag_up_    flag_dn_    chg_up_    chg_dn_
   any         any       dar_01      dar_01      dmr_01     dmr_01    COUNT

    0           0           0           0          0          0       14189
    0           0           0           0          0          1          27
    0           0           0           0          1          0           3
    0           0           0           1          0          0          46
    0           0           0           1          1          0           6
    0           0           1           0          0          0       15886
    0           0           1           0          0          1          35
    0           0           1           0          1          0           7
    0           0           1           1          0          0         123
    0           0           1           1          0          1           2
    0           0           1           1          1          0          10
    0           0           1           1          1          1           1
    0           1           0           0          0          0          77
    0           1           0           1          0          0           1
    0           1           1           0          0          0          58
    0           1           1           1          0          0           1
    1           0           0           0          0          0         240
    1           0           0           1          0          0           3
    1           0           1           0          0          0         219
    1           1           1           0          0          0           2


 flag_up_    flag_dn_                             flag_      flag_
  01_deg_     01_deg_    flag_up_    flag_dn_    chg_up_    chg_dn_
    72          72        dar_01      dar_01      dmr_01     dmr_01    COUNT

     0           0           0           0          0          0       14471
     0           0           0           0          0          1          27
     0           0           0           0          1          0           3
     0           0           0           1          0          0          50
     0           0           0           1          1          0           6
     0           0           1           0          0          0       16134
     0           0           1           0          0          1          35
     0           0           1           0          1          0           7
     0           0           1           1          0          0         124
     0           0           1           1          0          1           2
     0           0           1           1          1          0          10
     0           0           1           1          1          1           1
     0           1           0           0          0          0          12
     0           1           1           0          0          0           3
     1           0           0           0          0          0          23
     1           0           1           0          0          0          28

 flag_up_    flag_dn_                             flag_      flag_
  1_deg_      1_deg_     flag_up_    flag_dn_    chg_up_    chg_dn_
    any         any        dar_1       dar_1      dmr_1      dmr_1     COUNT

     0           0           0           0          0          0       21019
     0           0           0           0          0          1           1
     0           0           0           0          1          0           1
     0           0           0           1          0          0        1747
     0           0           0           1          1          0           1
     0           0           1           0          0          0        6790
     0           0           1           0          1          0           1
     0           0           1           1          0          0         840
     0           1           0           0          0          0          80
     0           1           0           1          0          0           6
     0           1           1           0          0          0          17
     0           1           1           1          0          0           1
     1           0           0           0          0          0         332
     1           0           0           1          0          0          16
     1           0           1           0          0          0          75
     1           0           1           1          0          0           9


                                                   flag_      flag_
  flag_up_    flag_dn_    flag_up_    flag_dn_    chg_up_    chg_dn_
  1_deg_72    1_deg_72      dar_1       dar_1      dmr_1      dmr_1     COUNT

      0           0           0           0          0          0       21400
      0           0           0           0          0          1           1
      0           0           0           0          1          0           1
      0           0           0           1          0          0        1767
      0           0           0           1          1          0           1
      0           0           1           0          0          0        6872
      0           0           1           0          1          0           1
      0           0           1           1          0          0         849
      0           1           0           0          0          0           4
      0           1           1           0          0          0           1
      1           0           0           0          0          0          27
      1           0           0           1          0          0           2
      1           0           1           0          0          0           9
      1           0           1           1          0          0           1


 flag_up_    flag_dn_                             flag_      flag_
  01_deg_     01_deg_    flag_up_    flag_dn_    chh_up_    chh_dn_
    any         any       dar_01      dar_01      dmr_01     dmr_01    COUNT

     0           0           0           0          0          0       12213
     0           0           0           0          0          1        1880
     0           0           0           0          1          0          78
     0           0           0           0          1          1          48
     0           0           0           1          0          0          30
     0           0           0           1          0          1          12
     0           0           0           1          1          0           9
     0           0           0           1          1          1           1
     0           0           1           0          0          0       12980
     0           0           1           0          0          1        2745
     0           0           1           0          1          0         107
     0           0           1           0          1          1          96
     0           0           1           1          0          0          68
     0           0           1           1          0          1          42
     0           0           1           1          1          0          14
     0           0           1           1          1          1          12
     0           1           0           0          0          0          64
     0           1           0           0          0          1          13
     0           1           0           1          0          1           1
     0           1           1           0          0          0          55
     0           1           1           0          0          1           3
     0           1           1           1          0          0           1
     1           0           0           0          0          0         213
     1           0           0           0          0          1          27
     1           0           0           1          0          0           3
     1           0           1           0          0          0         196
     1           0           1           0          0          1          22
     1           0           1           0          1          0           1
     1           1           1           0          0          1           2


  flag_up_    flag_dn_                             flag_      flag_
   01_deg_     01_deg_    flag_up_    flag_dn_    chh_up_    chh_dn_
     72          72        dar_01      dar_01      dmr_01     dmr_01    COUNT

      0           0           0           0          0          0       12459
      0           0           0           0          0          1        1916
      0           0           0           0          1          0          78
      0           0           0           0          1          1          48
      0           0           0           1          0          0          33
      0           0           0           1          0          1          13
      0           0           0           1          1          0           9
      0           0           0           1          1          1           1
      0           0           1           0          0          0       13203
      0           0           1           0          0          1        2769
      0           0           1           0          1          0         108
      0           0           1           0          1          1          96
      0           0           1           1          0          0          69
      0           0           1           1          0          1          42
      0           0           1           1          1          0          14
      0           0           1           1          1          1          12
      0           1           0           0          0          0           9
      0           1           0           0          0          1           3
      0           1           1           0          0          0           3
      1           0           0           0          0          0          22
      1           0           0           0          0          1           1
      1           0           1           0          0          0          25
      1           0           1           0          0          1           3

                                   The SAS System

flag_up_    flag_dn_                             flag_      flag_
 1_deg_      1_deg_     flag_up_    flag_dn_    chh_up_    chh_dn_
   any         any        dar_1       dar_1      dmr_1      dmr_1     COUNT

    0           0           0           0          0          0       19254
    0           0           0           0          0          1         409
    0           0           0           0          1          0        1165
    0           0           0           0          1          1         193
    0           0           0           1          0          0        1502
    0           0           0           1          0          1          38
    0           0           0           1          1          0         171
    0           0           0           1          1          1          37
    0           0           1           0          0          0        5934
    0           0           1           0          0          1         220
    0           0           1           0          1          0         527
    0           0           1           0          1          1         110
    0           0           1           1          0          0         658
    0           0           1           1          0          1          28
    0           0           1           1          1          0         121
    0           0           1           1          1          1          33
    0           1           0           0          0          0          76
    0           1           0           0          0          1           2
    0           1           0           0          1          0           2
    0           1           0           1          0          0           3
    0           1           0           1          1          0           2
    0           1           0           1          1          1           1
    0           1           1           0          0          0          17
    0           1           1           1          1          0           1
    1           0           0           0          0          0         312
    1           0           0           0          1          0          19
    1           0           0           0          1          1           1
    1           0           0           1          0          0          14
    1           0           0           1          0          1           1
    1           0           0           1          1          0           1
    1           0           1           0          0          0          67
    1           0           1           0          0          1           2
    1           0           1           0          1          0           6
    1           0           1           1          0          0           8
    1           0           1           1          1          0           1

                                 The SAS System
                                                 flag_      flag_
flag_up_    flag_dn_    flag_up_    flag_dn_    chh_up_    chh_dn_
1_deg_72    1_deg_72      dar_1       dar_1      dmr_1      dmr_1     COUNT

    0           0           0           0          0          0       19613
    0           0           0           0          0          1         410
    0           0           0           0          1          0        1185
    0           0           0           0          1          1         194
    0           0           0           1          0          0        1517
    0           0           0           1          0          1          39
    0           0           0           1          1          0         174
    0           0           0           1          1          1          38
    0           0           1           0          0          0        6008
    0           0           1           0          0          1         222
    0           0           1           0          1          0         533
    0           0           1           0          1          1         110
    0           0           1           1          0          0         666
    0           0           1           1          0          1          28
    0           0           1           1          1          0         122
    0           0           1           1          1          1          33
    0           1           0           0          0          0           3
    0           1           0           0          0          1           1
    0           1           1           0          0          0           1
    1           0           0           0          0          0          26
    1           0           0           0          1          0           1
    1           0           0           1          0          0           2
    1           0           1           0          0          0           9
    1           0           1           1          1          0           1

*/

proc export data=sig_table_w_go outfile="!HOME/concannon/DTRA/at_rad_deg_dar_dmr_sig_table_w_go.csv"
   dbms=csv replace;
run;
