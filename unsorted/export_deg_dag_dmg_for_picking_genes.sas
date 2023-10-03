/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


/* Find genes that are DEG, DAG, and DMG */


/* RNAseq genes */


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

data deg_01_1 deg_01_3 deg_01_24 deg_01_72
     deg1_01_1 deg1_01_3 deg1_01_24 deg1_01_72
     deg_1_1 deg_1_3 deg_1_24 deg_1_72
     deg1_1_1 deg1_1_3 deg1_1_24 deg1_1_72;
     set de_Results2;
     if (flag_01gy_v_mock_1h_fdr05=1 or (flag_01gy_1hr_on ^= flag_Mock_1hr_on )) then output deg_01_1;
     if (flag_01gy_v_mock_3h_fdr05=1 or (flag_01gy_3hr_on ^= flag_Mock_3hr_on )) then output deg_01_3;
     if (flag_01gy_v_mock_24h_fdr05=1 or (flag_01gy_24hr_on ^= flag_Mock_24hr_on )) then output deg_01_24;
     if (flag_01gy_v_mock_72h_fdr05=1 or (flag_01gy_72hr_on ^= flag_Mock_72hr_on )) then output deg_01_72;

     if (flag_1gy_v_mock_1h_fdr05=1 or (flag_1gy_1hr_on ^= flag_Mock_1hr_on )) then output deg_1_1;
     if (flag_1gy_v_mock_3h_fdr05=1 or (flag_1gy_3hr_on ^= flag_Mock_3hr_on )) then output deg_1_3;
     if (flag_1gy_v_mock_24h_fdr05=1 or (flag_1gy_24hr_on ^= flag_Mock_24hr_on )) then output deg_1_24;
     if (flag_1gy_v_mock_72h_fdr05=1 or (flag_1gy_72hr_on ^= flag_Mock_72hr_on )) then output deg_1_72;

     if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) then output deg1_01_1;
     if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) then output deg1_01_3;
     if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) then output deg1_01_24;
     if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) then output deg1_01_72;

     if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) then output deg1_1_1;
     if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) then output deg1_1_3;
     if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) then output deg1_1_24;
     if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) then output deg1_1_72;

keep gene_id;
run;


proc sort data=deg_01_1; by gene_id;     proc sort data=deg1_01_1; by gene_id;
proc sort data=deg_01_3; by gene_id;     proc sort data=deg1_01_3; by gene_id;
proc sort data=deg_01_24; by gene_id;     proc sort data=deg1_01_24; by gene_id;
proc sort data=deg_01_72; by gene_id;     proc sort data=deg1_01_72; by gene_id;
proc sort data=deg_1_1; by gene_id;     proc sort data=deg1_1_1; by gene_id;
proc sort data=deg_1_3; by gene_id;     proc sort data=deg1_1_3; by gene_id;
proc sort data=deg_1_24; by gene_id;     proc sort data=deg1_1_24; by gene_id;
proc sort data=deg_1_72; by gene_id;     proc sort data=deg1_1_72; by gene_id;
run;

/* DMRs */


data dmg_cg_01 dmg_cg_1 dmg_chg_01 dmg_chg_1 dmg_chh_01 dmg_chh_1;
     set arabMAP.results_by_dmr_annot;
     if geneID = "" then delete;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and site_type="CG" and comparison="0Gy_vs_01G" then output dmg_cg_01;
     if num_sites_FDR05_diff_10perc>=2 and site_type="CG" and comparison="0Gy_vs_1Gy" then output dmg_cg_1;
     if num_sites_FDR05_diff_10perc>=2 and site_type="CHG" and comparison="0Gy_vs_01G" then output dmg_chg_01;
     if num_sites_FDR05_diff_10perc>=2 and site_type="CHG" and comparison="0Gy_vs_1Gy" then output dmg_chg_1;
     if num_sites_FDR05_diff_10perc>=2 and site_type="CHH" and comparison="0Gy_vs_01G" then output dmg_chh_01;
     if num_sites_FDR05_diff_10perc>=2 and site_type="CHH" and comparison="0Gy_vs_1Gy" then output dmg_chh_1;
     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and site_type="CG" and comparison="0Gy_vs_01G" then output dmg_cg_01;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and site_type="CG" and comparison="0Gy_vs_1Gy" then output dmg_cg_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and site_type="CHG" and comparison="0Gy_vs_01G" then output dmg_chg_01;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and site_type="CHG" and comparison="0Gy_vs_1Gy" then output dmg_chg_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and site_type="CHH" and comparison="0Gy_vs_01G" then output dmg_chh_01;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and site_type="CHH" and comparison="0Gy_vs_1Gy" then output dmg_chh_1;
     keep geneID;
    rename geneID=gene_id;
run;



data dag_01 dag_1;
     set arabMAP.results_by_dar_annot;
     if geneID = "" then delete;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if comparison2="0Gy_vs_01G" then output dag_01;
     if comparison2="0Gy_vs_1Gy" then output dag_1;
     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison2="0Gy_vs_01G" then output dag_01;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison2="0Gy_vs_1Gy" then output dag_1;
     keep geneID;
    rename geneID=gene_id;
run;


proc sort data=dmg_cg_01 nodup; by gene_id;   proc sort data=dmg_cg_1 nodup; by gene_id;
proc sort data=dmg_chg_01 nodup; by gene_id;   proc sort data=dmg_chg_1 nodup; by gene_id;
proc sort data=dmg_chh_01 nodup; by gene_id;   proc sort data=dmg_chh_1 nodup; by gene_id;
proc sort data=dag_01 nodup; by gene_id;      proc sort data=dag_1 nodup; by gene_id;
run;

data compare_deg_dmg_dag;
  merge   deg_01_1 (in=in1) deg_01_3 (in=in2) deg_01_24 (in=in3) deg_01_72 (in=in4) deg_1_1 (in=in5) deg_1_3 (in=in6)
          deg_1_24 (in=in7) deg_1_72 (in=in8) deg1_01_1 (in=in9) deg1_01_3 (in=in10) deg1_01_24 (in=in11) deg1_01_72 (in=in12)
          deg1_1_1 (in=in13) deg1_1_3 (in=in14) deg1_1_24 (in=in15) deg1_1_72 (in=in16) dmg_cg_01 (in=in17) dmg_cg_1 (in=in18)
          dmg_chg_01 (in=in19) dmg_chg_1 (in=in20) dmg_chh_01 (in=in21) dmg_chh_1 (in=in22) dag_01 (in=in23) dag_1 (in=in24);
  by gene_id;
  if in1 then deg_01_1=1; else deg_01_1=0;
  if in2 then deg_01_3=1; else deg_01_3=0;
  if in3 then deg_01_24=1; else deg_01_24=0;
  if in4 then deg_01_72=1; else deg_01_72=0;
  if in5 then deg_1_1=1; else deg_1_1=0;
  if in6 then deg_1_3=1; else deg_1_3=0;
  if in7 then deg_1_24=1; else deg_1_24=0;
  if in8 then deg_1_72=1; else deg_1_72=0;
  if in9 then deg1_01_1=1; else deg1_01_1=0;
  if in10 then deg1_01_3=1; else deg1_01_3=0;
  if in11 then deg1_01_24=1; else deg1_01_24=0;
  if in12 then deg1_01_72=1; else deg1_01_72=0;
  if in13 then deg1_1_1=1; else deg1_1_1=0;
  if in14 then deg1_1_3=1; else deg1_1_3=0;
  if in15 then deg1_1_24=1; else deg1_1_24=0;
  if in16 then deg1_1_72=1; else deg1_1_72=0;
  if in17 then dmg_cg_01=1; else dmg_cg_01=0;
  if in18 then dmg_cg_1=1; else dmg_cg_1=0;
  if in19 then dmg_chg_01=1; else dmg_chg_01=0;
  if in20 then dmg_chg_1=1; else dmg_chg_1=0;
  if in21 then dmg_chh_01=1; else dmg_chh_01=0;
  if in22 then dmg_chh_1=1; else dmg_chh_1=0;
  if in23 then dag_01=1; else dag_01=0;
  if in24 then dag_1=1; else dag_1=0;
run;

proc freq data=compare_deg_dmg_dag noprint;
  tables deg_01_72*dmg_cg_01*dmg_chg_01*dmg_chh_01*dag_01 /out= ctabs_01_72;
  tables deg1_01_72*dmg_cg_01*dmg_chg_01*dmg_chh_01*dag_01 /out= ctabs_01_72_fc1;
  tables deg_1_72*dmg_cg_1*dmg_chg_1*dmg_chh_1*dag_1 /out= ctabs_1_72;
  tables deg1_1_72*dmg_cg_1*dmg_chg_1*dmg_chh_1*dag_1 /out= ctabs_1_72_fc1;
run;

proc export data=compare_deg_dmg_dag outfile="!HOME/concannon/DTRA/arabidopsis_radiation_DEG_DMG_DAG.csv"
   dbms=csv replace;
run;


