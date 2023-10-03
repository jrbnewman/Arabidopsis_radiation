libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Split results summaries by chrom, site type, sig */

/*


data dmc_sig_01 dmc_sig_1
     dmc_cg_01_1 dmc_cg_01_2 dmc_cg_01_3 dmc_cg_01_4 dmc_cg_01_5
     dmc_chg_01_1 dmc_chg_01_2 dmc_chg_01_3 dmc_chg_01_4 dmc_chg_01_5
     dmc_chh_01_1 dmc_chh_01_2 dmc_chh_01_3 dmc_chh_01_4 dmc_chh_01_5

     dmc_cg_1_1 dmc_cg_1_2 dmc_cg_1_3 dmc_cg_1_4 dmc_cg_1_5
     dmc_chg_1_1 dmc_chg_1_2 dmc_chg_1_3 dmc_chg_1_4 dmc_chg_1_5
     dmc_chh_1_1 dmc_chh_1_2 dmc_chh_1_3 dmc_chh_1_4 dmc_chh_1_5;
retain  comparison site_type chr start_pos stop_pos
methyl_CTL methyl_TRT methyl_diff FET_P FET_FDR_P
flag_p05 flag_fdr05 flag_fdr05_10perc flag_meth_diff flag_meth_diff_10perc;
  set wgbslocA.results_by_dmc_annot;
  start_pos=start_pos-1;
  stop_pos=stop_pos-1;
     if flag_fdr05_10perc=1 and comparison="0Gy_vs_01G" then output dmc_sig_01;
     if flag_fdr05_10perc=1 and comparison="0Gy_vs_1Gy" then output dmc_sig_1;
  if comparison="0Gy_vs_01G" and site_type="CG" and chr="1" then output dmc_cg_01_1;
  if comparison="0Gy_vs_01G" and site_type="CG" and chr="2" then output dmc_cg_01_2;
  if comparison="0Gy_vs_01G" and site_type="CG" and chr="3" then output dmc_cg_01_3;
  if comparison="0Gy_vs_01G" and site_type="CG" and chr="4" then output dmc_cg_01_4;
  if comparison="0Gy_vs_01G" and site_type="CG" and chr="5" then output dmc_cg_01_5;

  if comparison="0Gy_vs_01G" and site_type="CHG" and chr="1" then output dmc_chg_01_1;
  if comparison="0Gy_vs_01G" and site_type="CHG" and chr="2" then output dmc_chg_01_2;
  if comparison="0Gy_vs_01G" and site_type="CHG" and chr="3" then output dmc_chg_01_3;
  if comparison="0Gy_vs_01G" and site_type="CHG" and chr="4" then output dmc_chg_01_4;
  if comparison="0Gy_vs_01G" and site_type="CHG" and chr="5" then output dmc_chg_01_5;

  if comparison="0Gy_vs_01G" and site_type="CHH" and chr="1" then output dmc_chh_01_1;
  if comparison="0Gy_vs_01G" and site_type="CHH" and chr="2" then output dmc_chh_01_2;
  if comparison="0Gy_vs_01G" and site_type="CHH" and chr="3" then output dmc_chh_01_3;
  if comparison="0Gy_vs_01G" and site_type="CHH" and chr="4" then output dmc_chh_01_4;
  if comparison="0Gy_vs_01G" and site_type="CHH" and chr="5" then output dmc_chh_01_5;

  if comparison="0Gy_vs_1Gy" and site_type="CG" and chr="1" then output dmc_cg_1_1;
  if comparison="0Gy_vs_1Gy" and site_type="CG" and chr="2" then output dmc_cg_1_2;
  if comparison="0Gy_vs_1Gy" and site_type="CG" and chr="3" then output dmc_cg_1_3;
  if comparison="0Gy_vs_1Gy" and site_type="CG" and chr="4" then output dmc_cg_1_4;
  if comparison="0Gy_vs_1Gy" and site_type="CG" and chr="5" then output dmc_cg_1_5;

  if comparison="0Gy_vs_1Gy" and site_type="CHG" and chr="1" then output dmc_chg_1_1;
  if comparison="0Gy_vs_1Gy" and site_type="CHG" and chr="2" then output dmc_chg_1_2;
  if comparison="0Gy_vs_1Gy" and site_type="CHG" and chr="3" then output dmc_chg_1_3;
  if comparison="0Gy_vs_1Gy" and site_type="CHG" and chr="4" then output dmc_chg_1_4;
  if comparison="0Gy_vs_1Gy" and site_type="CHG" and chr="5" then output dmc_chg_1_5;

  if comparison="0Gy_vs_1Gy" and site_type="CHH" and chr="1" then output dmc_chh_1_1;
  if comparison="0Gy_vs_1Gy" and site_type="CHH" and chr="2" then output dmc_chh_1_2;
  if comparison="0Gy_vs_1Gy" and site_type="CHH" and chr="3" then output dmc_chh_1_3;
  if comparison="0Gy_vs_1Gy" and site_type="CHH" and chr="4" then output dmc_chh_1_4;
  if comparison="0Gy_vs_1Gy" and site_type="CHH" and chr="5" then output dmc_chh_1_5;
   drop flag_meth_diff;
run;
*/



data dmr_sig_01 dmr_sig_1 dmr_sig_01_bin dmr_sig_1_bin
     dmr_cg_01 dmr_chg_01 dmr_chh_01
     dmr_cg_1 dmr_chg_1 dmr_chh_1;
retain comparison site_type chr region_start region_stop
mean_methyl_CTL mean_methyl_TRT mean_methyl_diff
binomial_P binomial_FDR_P 
flag_binomial_p05 flag_binomial_fdr05
num_sites_in_region num_sites_FDR05 num_sites_diff_10perc num_sites_FDR05_diff_10perc;

  set wgbslocA.results_by_dmr_annot;
  region_start=region_start-1;
  region_stop=region_stop-1;

     if num_sites_FDR05_diff_10perc>=2 and comparison="0Gy_vs_01G" then output dmr_sig_01;
     if num_sites_FDR05_diff_10perc>=2 and comparison="0Gy_vs_1Gy" then output dmr_sig_1;
  if flag_binomial_FDR05=1 and abs(mean_methyl_diff)>=0.1 and comparison="0Gy_vs_01G" then output dmr_sig_01_bin;
  if flag_binomial_FDR05=1 and abs(mean_methyl_diff)>=0.1 and comparison="0Gy_vs_1Gy" then output dmr_sig_1_bin;
  if comparison="0Gy_vs_01G" and site_type="CG" then output dmr_cg_01;
  if comparison="0Gy_vs_01G" and site_type="CHG" then output dmr_chg_01;
  if comparison="0Gy_vs_01G" and site_type="CHH" then output dmr_chh_01;

  if comparison="0Gy_vs_1Gy" and site_type="CG" then output dmr_cg_1;
  if comparison="0Gy_vs_1Gy" and site_type="CHG" then output dmr_chg_1;
  if comparison="0Gy_vs_1Gy" and site_type="CHH" then output dmr_chh_1;
drop start end  region_num;
run;



/*

data dac_sig_01 dac_sig_1
     dac_cg_01_1 dac_cg_01_2 dac_cg_01_3 dac_cg_01_4 dac_cg_01_5
     dac_cg_1_1 dac_cg_1_2 dac_cg_1_3 dac_cg_1_4 dac_cg_1_5;

retain comparison2 site_type chr start_pos stop_pos
methyl_CTL_0U methyl_CTL_100U
methyl_TRT_0U methyl_TRT_100U
methyl_diff_CTL methyl_diff_TRT methyl_diff_TRT_CTL
FET_P_CTL_100U_0U FET_P_TRT_100U_0U FET_FDR_P_CTL_100U_0U FET_FDR_P_TRT_100U_0U
flag_CTL_fdr05 flag_TRT_fdr05 flag_FDR05_CTL_and_TRT flag_FDR05_CTL_or_TRT
flag_meth_diff_20perc flag_fdr05_20perc ;
  length comparison2 $20.;

  set wgbslocA.results_by_dac_annot;
  start_pos=start_pos-1;
  stop_pos=stop_pos-1;

  if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
  if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
  if comparison2="0Gy_vs_01G" and flag_fdr05_20perc=1 then output dac_sig_01;
  if comparison2="0Gy_vs_1Gy" and flag_fdr05_20perc=1 then output dac_sig_1;
  if comparison2="0Gy_vs_01G" and chr="1" then output dac_cg_01_1;
  if comparison2="0Gy_vs_01G" and chr="2" then output dac_cg_01_2;
  if comparison2="0Gy_vs_01G" and chr="3" then output dac_cg_01_3;
  if comparison2="0Gy_vs_01G" and chr="4" then output dac_cg_01_4;
  if comparison2="0Gy_vs_01G" and chr="5" then output dac_cg_01_5;

  if comparison2="0Gy_vs_1Gy" and chr="1" then output dac_cg_1_1;
  if comparison2="0Gy_vs_1Gy" and chr="2" then output dac_cg_1_2;
  if comparison2="0Gy_vs_1Gy" and chr="3" then output dac_cg_1_3;
  if comparison2="0Gy_vs_1Gy" and chr="4" then output dac_cg_1_4;
  if comparison2="0Gy_vs_1Gy" and chr="5" then output dac_cg_1_5;
  drop comparison flag_meth_diff;
  rename comparison2=comparison;
run;
*/

data 
     dar_cg_01 
     dar_cg_1 ;
retain comparison2 site_type chr region_start region_stop
binomial_P binomial_FDR_P 
flag_binomial_p05 flag_binomial_fdr05
num_sites_in_region
num_sites_FDR05_CTL num_sites_FDR05_TRT
num_sites_FDR05_CTL_and_TRT num_sites_FDR05_CTL_or_TRT
mean_methyl_CTL_0U mean_methyl_CTL_100U
mean_methyl_TRT_0U mean_methyl_TRT_100U
mean_methyl_diff_CTL_100U_0U mean_methyl_diff_TRT_100U_0U
mean_methyl_diff_TRT_CTL
num_sites_diff_20perc num_sites_FDR05_diff_20perc;

  set wgbslocA.results_by_dar_annot;
  region_start=region_start-1;
  region_stop=region_stop-1;
  if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
  if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
  if comparison2="0Gy_vs_01G" then output dar_cg_01;
  if comparison2="0Gy_vs_1Gy" then output dar_cg_1;
  drop comparison  region_num;
  rename comparison2=comparison;
run;

/*


data dmc_chh_1_1A dmc_chh_1_1B dmc_chh_1_1C; set dmc_chh_1_1; if _N_ <= 1000000 then output dmc_chh_1_1A; else if _N_ <= 2000000 then output dmc_chh_1_1B;  else output dmc_chh_1_1C; run;
data dmc_chh_1_2A dmc_chh_1_2B dmc_chh_1_2C; set dmc_chh_1_2; if _N_ <= 1000000 then output dmc_chh_1_2A; else if _N_ <= 2000000 then output dmc_chh_1_2B;  else output dmc_chh_1_2C; run;
data dmc_chh_1_3A dmc_chh_1_3B dmc_chh_1_3C; set dmc_chh_1_3; if _N_ <= 1000000 then output dmc_chh_1_3A; else if _N_ <= 2000000 then output dmc_chh_1_3B;  else output dmc_chh_1_3C; run;
data dmc_chh_1_4A dmc_chh_1_4B dmc_chh_1_4C; set dmc_chh_1_4; if _N_ <= 1000000 then output dmc_chh_1_4A; else if _N_ <= 2000000 then output dmc_chh_1_4B;  else output dmc_chh_1_4C; run;
data dmc_chh_1_5A dmc_chh_1_5B dmc_chh_1_5C; set dmc_chh_1_5; if _N_ <= 1000000 then output dmc_chh_1_5A; else if _N_ <= 2000000 then output dmc_chh_1_5B;  else output dmc_chh_1_5C; run;


data dmc_chh_01_1A dmc_chh_01_1B dmc_chh_01_1C; set dmc_chh_01_1; if _N_ <= 1000000 then output dmc_chh_01_1A; else if _N_ <= 2000000 then output dmc_chh_01_1B;  else output dmc_chh_01_1C; run;
data dmc_chh_01_2A dmc_chh_01_2B dmc_chh_01_2C; set dmc_chh_01_2; if _N_ <= 1000000 then output dmc_chh_01_2A; else if _N_ <= 2000000 then output dmc_chh_01_2B;  else output dmc_chh_01_2C; run;
data dmc_chh_01_3A dmc_chh_01_3B dmc_chh_01_3C; set dmc_chh_01_3; if _N_ <= 1000000 then output dmc_chh_01_3A; else if _N_ <= 2000000 then output dmc_chh_01_3B;  else output dmc_chh_01_3C; run;
data dmc_chh_01_4A dmc_chh_01_4B dmc_chh_01_4C; set dmc_chh_01_4; if _N_ <= 1000000 then output dmc_chh_01_4A; else if _N_ <= 2000000 then output dmc_chh_01_4B;  else output dmc_chh_01_4C; run;
data dmc_chh_01_5A dmc_chh_01_5B dmc_chh_01_5C; set dmc_chh_01_5; if _N_ <= 1000000 then output dmc_chh_01_5A; else if _N_ <= 2000000 then output dmc_chh_01_5B;  else output dmc_chh_01_5C; run;
*/
/* Export everything */


%macro exportData(inData,outFile);

proc export data=&inData. outfile="/TB14/TB14/sandbox/dtra_sandbox/&outFile..csv" dbms=csv replace; run; 

%mend;
/*
%exportData(dac_cg_01_1,Arabidopsis_results_by_site_GC_01Gy_v_0Gy_chr1);
%exportData(dac_cg_01_2,Arabidopsis_results_by_site_GC_01Gy_v_0Gy_chr2);
%exportData(dac_cg_01_3,Arabidopsis_results_by_site_GC_01Gy_v_0Gy_chr3);
%exportData(dac_cg_01_4,Arabidopsis_results_by_site_GC_01Gy_v_0Gy_chr4);
%exportData(dac_cg_01_5,Arabidopsis_results_by_site_GC_01Gy_v_0Gy_chr5);

%exportData(dac_cg_1_1,Arabidopsis_results_by_site_GC_1Gy_v_0Gy_chr1);
%exportData(dac_cg_1_2,Arabidopsis_results_by_site_GC_1Gy_v_0Gy_chr2);
%exportData(dac_cg_1_3,Arabidopsis_results_by_site_GC_1Gy_v_0Gy_chr3);
%exportData(dac_cg_1_4,Arabidopsis_results_by_site_GC_1Gy_v_0Gy_chr4);
%exportData(dac_cg_1_5,Arabidopsis_results_by_site_GC_1Gy_v_0Gy_chr5);

%exportData(dac_sig_1,Arabidopsis_results_by_site_GC_1Gy_v_0Gy_sig_only);
%exportData(dac_sig_01,Arabidopsis_results_by_site_GC_01Gy_v_0Gy_sig_only);


%exportData(dmc_cg_01_1,Arabidopsis_results_by_site_CG_01Gy_v_0Gy_chr1);
%exportData(dmc_cg_01_2,Arabidopsis_results_by_site_CG_01Gy_v_0Gy_chr2);
%exportData(dmc_cg_01_3,Arabidopsis_results_by_site_CG_01Gy_v_0Gy_chr3);
%exportData(dmc_cg_01_4,Arabidopsis_results_by_site_CG_01Gy_v_0Gy_chr4);
%exportData(dmc_cg_01_5,Arabidopsis_results_by_site_CG_01Gy_v_0Gy_chr5);

%exportData(dmc_cg_1_1,Arabidopsis_results_by_site_CG_1Gy_v_0Gy_chr1);
%exportData(dmc_cg_1_2,Arabidopsis_results_by_site_CG_1Gy_v_0Gy_chr2);
%exportData(dmc_cg_1_3,Arabidopsis_results_by_site_CG_1Gy_v_0Gy_chr3);
%exportData(dmc_cg_1_4,Arabidopsis_results_by_site_CG_1Gy_v_0Gy_chr4);
%exportData(dmc_cg_1_5,Arabidopsis_results_by_site_CG_1Gy_v_0Gy_chr5);


%exportData(dmc_chg_01_1,Arabidopsis_results_by_site_CHG_01Gy_v_0Gy_chr1);
%exportData(dmc_chg_01_2,Arabidopsis_results_by_site_CHG_01Gy_v_0Gy_chr2);
%exportData(dmc_chg_01_3,Arabidopsis_results_by_site_CHG_01Gy_v_0Gy_chr3);
%exportData(dmc_chg_01_4,Arabidopsis_results_by_site_CHG_01Gy_v_0Gy_chr4);
%exportData(dmc_chg_01_5,Arabidopsis_results_by_site_CHG_01Gy_v_0Gy_chr5);

%exportData(dmc_chg_1_1,Arabidopsis_results_by_site_CHG_1Gy_v_0Gy_chr1);
%exportData(dmc_chg_1_2,Arabidopsis_results_by_site_CHG_1Gy_v_0Gy_chr2);
%exportData(dmc_chg_1_3,Arabidopsis_results_by_site_CHG_1Gy_v_0Gy_chr3);
%exportData(dmc_chg_1_4,Arabidopsis_results_by_site_CHG_1Gy_v_0Gy_chr4);
%exportData(dmc_chg_1_5,Arabidopsis_results_by_site_CHG_1Gy_v_0Gy_chr5);


%exportData(dmc_chh_01_1a,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr1_part1);
%exportData(dmc_chh_01_2a,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr2_part1);
%exportData(dmc_chh_01_3a,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr3_part1);
%exportData(dmc_chh_01_4a,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr4_part1);
%exportData(dmc_chh_01_5a,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr5_part1);

%exportData(dmc_chh_1_1a,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr1_part1);
%exportData(dmc_chh_1_2a,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr2_part1);
%exportData(dmc_chh_1_3a,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr3_part1);
%exportData(dmc_chh_1_4a,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr4_part1);
%exportData(dmc_chh_1_5a,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr5_part1);


%exportData(dmc_chh_01_1b,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr1_part2);
%exportData(dmc_chh_01_2b,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr2_part2);
%exportData(dmc_chh_01_3b,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr3_part2);
%exportData(dmc_chh_01_4b,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr4_part2);
%exportData(dmc_chh_01_5b,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr5_part2);

%exportData(dmc_chh_1_1b,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr1_part2);
%exportData(dmc_chh_1_2b,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr2_part2);
%exportData(dmc_chh_1_3b,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr3_part2);
%exportData(dmc_chh_1_4b,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr4_part2);
%exportData(dmc_chh_1_5b,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr5_part2);


%exportData(dmc_chh_01_1c,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr1_part3);
%exportData(dmc_chh_01_2c,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr2_part3);
%exportData(dmc_chh_01_3c,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr3_part3);
%exportData(dmc_chh_01_4c,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr4_part3);
%exportData(dmc_chh_01_5c,Arabidopsis_results_by_site_CHH_01Gy_v_0Gy_chr5_part3);

%exportData(dmc_chh_1_1c,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr1_part3);
%exportData(dmc_chh_1_2c,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr2_part3);
%exportData(dmc_chh_1_3c,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr3_part3);
%exportData(dmc_chh_1_4c,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr4_part3);
%exportData(dmc_chh_1_5c,Arabidopsis_results_by_site_CHH_1Gy_v_0Gy_chr5_part3);

%exportData(dmc_sig_01,Arabidopsis_results_by_site_01Gy_v_0Gy_sig_only);
%exportData(dmc_sig_1,Arabidopsis_results_by_site_1Gy_v_0Gy_sig_only);
*/
%exportData(dmr_sig_01,Arabidopsis_results_by_region_01Gy_v_0Gy_2_sig_sites);
%exportData(dmr_sig_1,Arabidopsis_results_by_region_1Gy_v_0Gy_2_sig_sites);
%exportData(dmr_sig_01_bin,Arabidopsis_results_by_region_01Gy_v_0Gy_binomial_sig);
%exportData(dmr_sig_1_bin,Arabidopsis_results_by_region_1Gy_v_0Gy_binomial_sig);

%exportData(dmr_cg_01,Arabidopsis_results_by_region_CG_01Gy_v_0Gy_all);
%exportData(dmr_chg_01,Arabidopsis_results_by_region_CHG_01Gy_v_0Gy_all);
%exportData(dmr_chh_01,Arabidopsis_results_by_region_CHH_01Gy_v_0Gy_all);
%exportData(dmr_cg_1,Arabidopsis_results_by_region_CG_1Gy_v_0Gy_all);
%exportData(dmr_chg_1,Arabidopsis_results_by_region_CHG_1Gy_v_0Gy_all);
%exportData(dmr_chh_1,Arabidopsis_results_by_region_CHH_1Gy_v_0Gy_all);

%exportData(dar_cg_01,Arabidopsis_results_by_region_GC_01Gy_v_0Gy_all);
%exportData(dar_cg_1,Arabidopsis_results_by_region_GC_1Gy_v_0Gy_all);

