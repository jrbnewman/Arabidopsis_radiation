libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


/* Count DMRs */

data meth2region;
  set wgbslocA.cytosine_to_meth_region_ge2_v4;
  if FET_FDR_P = .  then flag_FET_FDR05=.;
  else if FET_FDR_P < 0.05 then flag_FET_FDR05=1;
  else flag_FET_FDR05=0;
  if flag_FET_FDR05=1 and flag_meth_diff_10perc=1 then flag_FET_FDR05_methdiff_10=1;
  else  flag_FET_FDR05_methdiff_10=0;
  if flag_FET_FDR05=1 and flag_meth_diff_20perc=1 then flag_FET_FDR05_methdiff_20=1;
  else  flag_FET_FDR05_methdiff_20=0;
  keep region_num site_type chr start_pos stop_pos comparison flag_FET_FDR05 flag_FET_FDR05_methdiff_10 flag_FET_FDR05_methdiff_20
   flag_meth_diff_20perc flag_meth_diff_10perc methyl_diff;
run;

data region_flags;
   set wgbslocA.meth_region_flag_site_ge2_v4;
   keep comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
run;


proc sort data=meth2region;
  by comparison site_type chr region_num;
proc sort data=region_flags;
  by comparison site_type chr region_num;
run;

data region_w_dmc;
  merge region_flags (in=in1) meth2region (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


proc sort data=region_w_dmc ;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
proc means data=region_w_dmc noprint;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
    var flag_FET_FDR05 flag_FET_FDR05_methdiff_10 flag_FET_FDR05_methdiff_20 flag_meth_diff_20perc flag_meth_diff_10perc;
    output out=dmc_per_dmr sum=;
run;

data dmc_per_dmr2;
  set dmc_per_dmr;
  if flag_FET_FDR05_methdiff_10 >= 1 then flag_DMC10_ge1=1; else flag_DMC10_ge1=0;
  if flag_FET_FDR05_methdiff_10 >= 2 then flag_DMC10_ge2=1; else flag_DMC10_ge2=0;
  if flag_FET_FDR05_methdiff_10 >= 3 then flag_DMC10_ge3=1; else flag_DMC10_ge3=0;
  if flag_FET_FDR05_methdiff_10 >= 5 then flag_DMC10_ge5=1; else flag_DMC10_ge5=0;
  if flag_FET_FDR05_methdiff_10 >= 10 then flag_DMC10_ge10=1; else flag_DMC10_ge10=0;

  if flag_FET_FDR05_methdiff_20 >= 1 then flag_DMC20_ge1=1; else flag_DMC20_ge1=0;
  if flag_FET_FDR05_methdiff_20 >= 2 then flag_DMC20_ge2=1; else flag_DMC20_ge2=0;
  if flag_FET_FDR05_methdiff_20 >= 3 then flag_DMC20_ge3=1; else flag_DMC20_ge3=0;
  if flag_FET_FDR05_methdiff_20 >= 5 then flag_DMC20_ge5=1; else flag_DMC20_ge5=0;
  if flag_FET_FDR05_methdiff_20 >= 10 then flag_DMC20_ge10=1; else flag_DMC20_ge10=0;
run;


/* Binomial test */

data binomial;
  set wgbslocA.bin2_dmr_anova_1 wgbslocA.bin2_dmr_anova_2 wgbslocA.bin2_dmr_anova_3 wgbslocA.bin2_dmr_anova_4 wgbslocA.bin2_dmr_anova_5;
run;

proc means data=region_w_dmc noprint;
    by comparison site_type chr  region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
    var methyl_diff start_pos stop_pos;
    output out=dmr_methdiff mean(methyl_diff)=methyl_diff min(start_pos)=start max(stop_pos)=end;
run;

proc sort data=binomial;
  by comparison site_type chr region_num;
proc sort data=dmr_methdiff;
  by comparison site_type chr region_num;
run;

data binomial2region;
  merge dmr_methdiff (in=in1) binomial (in=in2);
  by comparison site_type chr region_num;
  if in1 ;
run;

proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge2 noprint;
  where flag_num_sites_ge2=1;
  by comparison site_type;
run;

proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge3 noprint;
  where flag_num_sites_ge3=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge10 noprint;
  where flag_num_sites_ge10=1;
  by comparison site_type;
run;

data flag_fdr_ge2;
set binom_fdr_ge2;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge3;
set binom_fdr_ge3;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge5;
set binom_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge10;
set binom_fdr_ge10;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


proc freq data=flag_fdr_ge2 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge3 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;

proc freq data=flag_fdr_ge5 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge10 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;

proc print data=c1;
proc print data=c2;
run;



/* linear model */

data mixed;
  set wgbslocA.mixed_dmr_anova_1 wgbslocA.mixed_dmr_anova_2 wgbslocA.mixed_dmr_anova_3 wgbslocA.mixed_dmr_anova_4 wgbslocA.mixed_dmr_anova_5;
run;

proc means data=region_w_dmc noprint;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
    var methyl_diff;
    output out=dmr_methdiff_lin mean=;
run;

proc sort data=mixed;
  by comparison site_type chr region_num;
proc sort data=dmr_methdiff_lin;
  by comparison site_type chr region_num;
run;

data mixed2region;
  merge dmr_methdiff (in=in1) mixed (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

proc multtest inpvalues(probf)=mixed2region fdr out=mixed_fdr_ge2 noprint;
  where flag_num_sites_ge2=1;
  by comparison site_type;
run;

proc multtest inpvalues(probf)=mixed2region fdr out=mixed_fdr_ge3 noprint;
  where flag_num_sites_ge3=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=mixed2region fdr out=mixed_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=mixed2region fdr out=mixed_fdr_ge10 noprint;
  where flag_num_sites_ge10=1;
  by comparison site_type;
run;

data flag_fdr_ge2_lin;
set mixed_fdr_ge2;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge3_lin;
set mixed_fdr_ge3;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge5_lin;
set mixed_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge10_lin;
set mixed_fdr_ge10;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


proc freq data=flag_fdr_ge2_lin noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge3_lin noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;

proc freq data=flag_fdr_ge5_lin noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge10_lin noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;

proc print data=c1;
proc print data=c2;
run;

data flag_fdr_ge5_2;
  set flag_fdr_ge5;
  rename flag_fdr05=flag_binomial_fdr05;
run;

proc sort data=flag_fdr_ge5_lin;
  by comparison site_type chr region_num;
proc sort data=flag_fdr_ge5_2;
  by comparison site_type chr region_num;
run;

data dmr_lin_vs_binom5;
 merge flag_fdr_ge5_lin flag_fdr_ge5_2;
  by comparison site_type chr region_num;
run;

proc freq data=dmr_lin_vs_binom5 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05*flag_binomial_fdr05 / out=c1;
run;

proc print data=c1;
run;




%macro checkDMR(numSites);

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_CG_minSites_&numSites..DMR.output.txt" out=met_01_CG dbms=tab replace; guessingrows=50000; getnames=no; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_CHG_minSites_&numSites..DMR.output.txt" out=met_01_CHG dbms=tab replace; guessingrows=50000;getnames=no; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_CHH_minSites_&numSites..DMR.output.txt" out=met_01_CHH dbms=tab replace; guessingrows=50000; getnames=no;run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_CG_minSites_&numSites..DMR.output.txt" out=met_1_CG dbms=tab replace; guessingrows=50000; getnames=no;run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_CHG_minSites_&numSites..DMR.output.txt" out=met_1_CHG dbms=tab replace; guessingrows=50000; getnames=no;run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_CHH_minSites_&numSites..DMR.output.txt" out=met_1_CHH dbms=tab replace; guessingrows=50000; getnames=no;run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CG_minSites_&numSites..DMR.results.overdispersion.txt" out=mk_1_CG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_01gy_CG_minSites_&numSites..DMR.results.overdispersion.txt" out=mk_01_CG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_&numSites..DMR.results.overdispersion.txt" out=mk_1_CHG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_01gy_CHG_minSites_&numSites..DMR.results.overdispersion.txt" out=mk_01_CHG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHH_minSites_&numSites..DMR.results.overdispersion.txt" out=mk_1_CHH dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_01gy_CHH_minSites_&numSites..DMR.results.overdispersion.txt" out=mk_01_CHH dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CG_minSites_&numSites..binomial.DMR.txt" out=ms_bin_01_CG  dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CG_minSites_&numSites..betabinomial.DMR.txt" out=ms_bb_01_CG  dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_&numSites..binomial.DMR.txt" out=ms_bin_1_CG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_&numSites..betabinomial.DMR.txt" out=ms_bb_1_CG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHG_minSites_&numSites..binomial.DMR.txt" out=ms_bin_01_CHG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHG_minSites_&numSites..betabinomial.DMR.txt" out=ms_bb_01_CHG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CHG_minSites_&numSites..binomial.DMR.txt" out=ms_bin_1_CHG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CHG_minSites_&numSites..betabinomial.DMR.txt" out=ms_bb_1_CHG dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHH_minSites_&numSites..binomial.DMR.txt" out=ms_bin_01_CHH dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHH_minSites_&numSites..betabinomial.DMR.txt" out=ms_bb_01_CHH dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CHH_minSites_&numSites..binomial.DMR.txt" out=ms_bin_1_CHH dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CHH_minSites_&numSites..betabinomial.DMR.txt" out=ms_bb_1_CHH dbms=tab replace; guessingrows=50000; run;


data met_01_cg2; set met_01_cg; var9_2=var9 * 1; drop var9; rename var9_2=var9; run;
data met_01_chg2; set met_01_chg; var9_2=var9 * 1; drop var9; rename var9_2=var9; run;
data met_01_chh2; set met_01_chh; var9_2=var9 * 1; drop var9; rename var9_2=var9; run;

data met_1_cg2; set met_1_cg; var9_2=var9 * 1; drop var9; rename var9_2=var9; run;
data met_1_chg2; set met_1_chg; var9_2=var9 * 1; drop var9; rename var9_2=var9; run;
data met_1_chh2; set met_1_chh; var9_2=var9 * 1; drop var9; rename var9_2=var9; run;



data metilene_all;
   set met_01_CG2 (in=in1) met_01_CHG2 (in=in2) met_01_CHH2 (in=in3) 
       met_1_CG2 (in=in4) met_1_CHG2 (in=in5) met_1_CHH2 (in=in6) ;
   length comparison $12.;
   length site_type $4.;
   if in1 or in2 or in3 then comparison="0Gy_vs_01G";
   if in4 or in5 or in6 then comparison="0Gy_vs_1Gy";
   if in1 or in4 then site_type="CG";
   if in2 or in5 then site_type="CHG";
   if in3 or in6 then site_type="CHH";
   chr2 = compress(VAR1);
   rename chr2=chr VAR2=start VAR3=end VAR4=metilene_q VAR5=metilene_meth_diff VAR6=num_sites VAR7=metilene_MWU_P;
run;




 proc sort data=flag_fdr_ge&numSites.;
  by comparison site_type chr start end;
 proc sort data=metilene_all;
  by comparison site_type chr  start end;
run;


data met_w_diff;
  merge flag_fdr_ge&numSites. (in=in1) metilene_all (in=in2);
  by comparison site_type chr start end;
  if in1;
  rename fdr_p=binomial_fdr_p;
run;

proc multtest inpvalues(metilene_MWU_P)=met_w_diff fdr out=met_fdr noprint;
  by comparison site_type;
run;

data met_fdr2;
  set met_fdr;
if fdr_p = . then flag_met_fdr05=.;
else if fdr_p < 0.05 then flag_met_fdr05=1;
else  flag_met_fdr05=0;

if metilene_q = . then flag_met_q05=.;
else if metilene_q < 0.05 then flag_met_q05=1;
else  flag_met_q05=0;

if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;

run;

proc freq data=met_fdr2 noprint;
  table comparison*site_type*flaG_methdiff10*flag_met_q05 / out=c1;
  table comparison*site_type*flaG_methdiff20*flag_met_q05 / out=c2;
  table comparison*site_type*flaG_methdiff10*flag_met_q05*flag_fdr05 / out=c3;
run;

proc print data=c1;
proc print data=c2;
proc print data=c3;
run;










data mk_all;
   set mk_01_CG (in=in1) mk_01_CHG (in=in2) mk_01_CHH (in=in3) 
       mk_1_CG (in=in4) mk_1_CHG (in=in5) mk_1_CHH (in=in6) ;
   length comparison $12.;
   length site_type $4.;
   if in1 or in2 or in3 then comparison="0Gy_vs_01G";
   if in4 or in5 or in6 then comparison="0Gy_vs_1Gy";
   if in1 or in4 then site_type="CG";
   if in2 or in5 then site_type="CHG";
   if in3 or in6 then site_type="CHH";
   start = start -1;
run;


 


 proc sort data=flag_fdr_ge&numSites.;
  by comparison site_type chr start end;
 proc sort data=mk_all;
  by comparison site_type chr  start end;
run;


data mk_w_diff;
  merge flag_fdr_ge&numSites. (in=in1) mk_all (in=in2);
  by comparison site_type chr start end;
  if in1 ;
  rename fdr_p=binomial_fdr_p;
run;

proc multtest inpvalues(pvalue)=mk_w_diff fdr out=met_fdr noprint;
  by comparison site_type;
run;

data met_fdr2;
  set met_fdr;
if fdr_p = . then flag_mk_fdr05=.;
else if fdr_p < 0.05 then flag_mk_fdr05=1;
else  flag_mk_fdr05=0;

if qvalue = . then flag_mk_q05=.;
else if qvalue < 0.05 then flag_mk_q05=1;
else  flag_mk_q05=0;

if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;

run;

proc freq data=met_fdr2 noprint;
  table comparison*site_type*flaG_methdiff10*flag_mk_q05 / out=c1;
  table comparison*site_type*flaG_methdiff20*flag_mk_q05 / out=c2;
  table comparison*site_type*flaG_methdiff10*flag_mk_q05*flag_fdr05 / out=c3;
run;

proc print data=c1;
proc print data=c2;
proc print data=c3;
run;



data ms_bin_01_CG2; set ms_bin_01_CG; pvalue2=pvalue * 1; fdr2=fdr*1; drop pvalue fdr; run;
data ms_bin_01_CHG2; set ms_bin_01_CHG; pvalue2=pvalue * 1; fdr2=fdr*1; drop pvalue fdr; run;
data ms_bin_01_CHH2; set ms_bin_01_CHH; pvalue2=pvalue * 1; fdr2=fdr*1; drop pvalue fdr; run;
data ms_bin_1_CG2; set ms_bin_1_CG; pvalue2=pvalue * 1; fdr2=fdr*1; drop  pvalue fdr; run;
data ms_bin_1_CHG2; set ms_bin_1_CHG; pvalue2=pvalue * 1;fdr2=fdr*1;  drop pvalue fdr; run;
data ms_bin_1_CHH2; set ms_bin_1_CHH; pvalue2=pvalue * 1; fdr2=fdr*1; drop  pvalue fdr; run;




data ms_bin_all;
   set ms_bin_01_CG2 (in=in1) ms_bin_01_CHG2 (in=in2) ms_bin_01_CHH2 (in=in3) 
       ms_bin_1_CG2 (in=in4) ms_bin_1_CHG2 (in=in5) ms_bin_1_CHH2 (in=in6) ;
   length comparison $12.;
   length site_type $4.;
   if in1 or in2 or in3 then comparison="0Gy_vs_01G";
   if in4 or in5 or in6 then comparison="0Gy_vs_1Gy";
   if in1 or in4 then site_type="CG";
   if in2 or in5 then site_type="CHG";
   if in3 or in6 then site_type="CHH";
   start = start - 1;
  rename seqnames=chr;
run;





 proc sort data=flag_fdr_ge&numSites.;
  by comparison site_type chr start end;
 proc sort data=ms_bin_all;
  by comparison site_type chr  start end;
run;


data ms_bin_w_diff;
  merge flag_fdr_ge&numSites. (in=in1) ms_bin_all (in=in2);
  by comparison site_type chr start end;
  if in1 ;
  rename fdr_p=binomial_fdr_p;
run;

proc multtest inpvalues(pvalue2)=ms_bin_w_diff fdr out=met_fdr noprint;
  by comparison site_type;
run;

data met_fdr2;
  set met_fdr;
if fdr_p = . then flag_bin_fdr05=.;
else if fdr_p < 0.05 then flag_bin_fdr05=1;
else  flag_bin_fdr05=0;


if fdr2 = . then flag_msbin_fdr05=.;
else if fdr2 < 0.05 then flag_msbin_fdr05=1;
else  flag_msbin_fdr05=0;

if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;

proc freq data=met_fdr2 noprint;
  table comparison*site_type*flaG_methdiff10*flag_msbin_fdr05 / out=c1;
  table comparison*site_type*flaG_methdiff20*flag_msbin_fdr05 / out=c2;
  table comparison*site_type*flaG_methdiff10*flag_msbin_fdr05*flag_fdr05 / out=c3;
run;

proc print data=c1;
proc print data=c2;
proc print data=c3;
run;






data ms_bb_01_CG2; set ms_bb_01_CG; pvalue2=pvalue * 1;  fdr2=fdr*1;drop pvalue fdr; run;
data ms_bb_01_CHG2; set ms_bb_01_CHG; pvalue2=pvalue * 1;  fdr2=fdr*1;drop pvalue fdr; run;
data ms_bb_01_CHH2; set ms_bb_01_CHH; pvalue2=pvalue * 1;  fdr2=fdr*1;drop pvalue fdr; run;
data ms_bb_1_CG2; set ms_bb_1_CG; pvalue2=pvalue * 1;  fdr2=fdr*1;drop  pvalue fdr; run;
data ms_bb_1_CHG2; set ms_bb_1_CHG; pvalue2=pvalue * 1; fdr2=fdr*1; drop pvalue fdr; run;
data ms_bb_1_CHH2; set ms_bb_1_CHH; pvalue2=pvalue * 1;  fdr2=fdr*1;drop  pvalue fdr; run;




data ms_bb_all;
   set ms_bb_01_CG2 (in=in1) ms_bb_01_CHG2 (in=in2) ms_bb_01_CHH2 (in=in3) 
       ms_bb_1_CG2 (in=in4) ms_bb_1_CHG2 (in=in5) ms_bb_1_CHH2 (in=in6) ;
   length comparison $12.;
   length site_type $4.;
   if in1 or in2 or in3 then comparison="0Gy_vs_01G";
   if in4 or in5 or in6 then comparison="0Gy_vs_1Gy";
   if in1 or in4 then site_type="CG";
   if in2 or in5 then site_type="CHG";
   if in3 or in6 then site_type="CHH";
   start = start - 1;
  rename seqnames=chr;
run;





 proc sort data=flag_fdr_ge&numSites.;
  by comparison site_type chr start end;
 proc sort data=ms_bb_all;
  by comparison site_type chr  start end;
run;


data ms_bb_w_diff;
  merge flag_fdr_ge&numSites. (in=in1) ms_bin_all (in=in2);
  by comparison site_type chr start end;
  if in1 ;
  rename fdr_p=binomial_fdr_p;
run;

proc multtest inpvalues(pvalue2)=ms_bb_w_diff fdr out=met_fdr noprint;
  by comparison site_type;
run;

data met_fdr2;
  set met_fdr;
if fdr_p = . then flag_bb_fdr05=.;
else if fdr_p < 0.05 then flag_bb_fdr05=1;
else  flag_bb_fdr05=0;

if fdr2 = . then flag_msbb_fdr05=.;
else if fdr2 < 0.05 then flag_msbb_fdr05=1;
else  flag_msbb_fdr05=0;


if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;

proc freq data=met_fdr2 noprint;
  table comparison*site_type*flaG_methdiff10*flag_msbb_fdr05 / out=c1;
  table comparison*site_type*flaG_methdiff20*flag_msbb_fdr05 / out=c2;
  table comparison*site_type*flaG_methdiff10*flag_msbb_fdr05*flag_fdr05 / out=c3;
run;

proc print data=c1;
proc print data=c2;
run;


%mend;


*%checkDMR(2);
*%checkDMR(3);
*%checkDMR(5);
*%checkDMR(10);

/* DARs */


proc datasets lib=work kill noprint;
run;
quit;




/* Count DMRs */

data meth2region;
  set wgbslocA.cytosine_to_acc_region_ge2_v4;
run;

data cmh_test;
  set wgbsA.all_dac_tests_01gy_compare (in=in1) wgbsA.all_dac_tests_1gy_compare (in=in2);
  length comparison $12.;
  if in1 then do;
     comparison = "01Gy_0Gy";
     cmh_fdr_p = CMH_FDR_P_01Gy_0gy_v2;
     end;
  if in2 then do; comparison = "1Gy_0Gy";
     cmh_fdr_p = CMH_FDR_P_1Gy_0gy_v2;
     end;
  keep comparison site_type chr start_pos stop_pos cmh_fdr_p;
run;

proc sort data=meth2region;
  by comparison site_type chr start_pos stop_pos;
proc sort data=cmh_test;
  by comparison site_type chr start_pos stop_pos;
run;

data meth2region2;
  merge meth2region (in=in1) cmh_test (in=in2);
  by comparison site_type chr start_pos stop_pos;
  if in1;
  if cmh_fdr_p = . then flag_CMH_FDR05=.;
  else if cmh_fdr_p < 0.05 then flag_CMH_FDR05=1;
  else flag_CMH_FDR05=0;
  if flag_meth_diff_10perc=1 and flag_CMH_FDR05=1 then flag_CHM_FDR05_10perc=1; else flag_CHM_FDR05_10perc=0;
  if flag_meth_diff_20perc=1 and flag_CMH_FDR05=1 then flag_CHM_FDR05_20perc=1; else flag_CHM_FDR05_20perc=0;
run;




data region_flags;
   set wgbslocA.acc_region_flag_site_ge2_v4;
   keep comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
run;

proc freq data=region_flags;
  tables comparison * (flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10);
run;



proc sort data=meth2region2;
  by comparison site_type chr region_num;
proc sort data=region_flags;
  by comparison site_type chr region_num;
run;

data region_w_dmc;
  merge region_flags (in=in1) meth2region2 (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


proc sort data=region_w_dmc ;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
proc means data=region_w_dmc noprint;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
    var flag_CMH_FDR05 flag_CHM_FDR05_10perc flag_CHM_FDR05_20perc flag_meth_diff_20perc flag_meth_diff_10perc;
    output out=dmc_per_dmr sum=;
run;

data dmc_per_dmr2;
  set dmc_per_dmr;
  if flag_CHM_FDR05_10perc >= 1 then flag_DMC10_ge1=1; else flag_DMC10_ge1=0;
  if flag_CHM_FDR05_10perc >= 2 then flag_DMC10_ge2=1; else flag_DMC10_ge2=0;
  if flag_CHM_FDR05_10perc >= 3 then flag_DMC10_ge3=1; else flag_DMC10_ge3=0;
  if flag_CHM_FDR05_10perc >= 5 then flag_DMC10_ge5=1; else flag_DMC10_ge5=0;
  if flag_CHM_FDR05_10perc >= 10 then flag_DMC10_ge10=1; else flag_DMC10_ge10=0;

  if flag_CHM_FDR05_20perc >= 1 then flag_DMC20_ge1=1; else flag_DMC20_ge1=0;
  if flag_CHM_FDR05_20perc >= 2 then flag_DMC20_ge2=1; else flag_DMC20_ge2=0;
  if flag_CHM_FDR05_20perc >= 3 then flag_DMC20_ge3=1; else flag_DMC20_ge3=0;
  if flag_CHM_FDR05_20perc >= 5 then flag_DMC20_ge5=1; else flag_DMC20_ge5=0;
  if flag_CHM_FDR05_20perc >= 10 then flag_DMC20_ge10=1; else flag_DMC20_ge10=0;
run;
/* Binomial test */

data binomial;
  set wgbslocA.bin2_dar_anova2_all;
  where effect = "group*units";
run;

proc means data=region_w_dmc noprint;
    by comparison site_type chr  region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
    var methyl_diff_TRT_CTL start_pos stop_pos;
    output out=dmr_methdiff mean(methyl_diff_TRT_CTL)=methyl_diff min(start_pos)=start max(stop_pos)=end;
run;

proc sort data=binomial;
  by comparison site_type chr region_num;
proc sort data=dmr_methdiff;
  by comparison site_type chr region_num;
run;

data binomial2region;
  merge dmr_methdiff (in=in1) binomial (in=in2);
  by comparison site_type chr region_num;
  if in1;
run;

proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge2 noprint;
  where flag_num_sites_ge2=1;
  by comparison site_type;
run;

proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge3 noprint;
  where flag_num_sites_ge3=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge10 noprint;
  where flag_num_sites_ge10=1;
  by comparison site_type;
run;

data flag_fdr_ge2;
set binom_fdr_ge2;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge3;
set binom_fdr_ge3;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge5;
set binom_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge10;
set binom_fdr_ge10;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;




proc freq data=flag_fdr_ge2 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge3 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;

proc freq data=flag_fdr_ge5 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge10 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;

proc print data=c1;
proc print data=c2;
run;




/* linear model */

data mixed1;
  set wgbslocA.mixedmeth_dar_anova_all;
   where effect="group*units";
run;


proc sort data=mixed1;
  by comparison site_type chr region_num;
proc sort data=dmr_methdiff;
  by comparison site_type chr region_num;
run;

data mixed2region;
  merge dmr_methdiff (in=in1) mixed1 (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge2 noprint;
  where flag_num_sites_ge2=1;
  by comparison site_type;
run;

proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge3 noprint;
  where flag_num_sites_ge3=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge10 noprint;
  where flag_num_sites_ge10=1;
  by comparison site_type;
run;

data flag_fdr_ge2_m1;
set binom_fdr_ge2;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge3_m1;
set binom_fdr_ge3;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge5_m1;
set binom_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge10_m1;
set binom_fdr_ge10;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;





proc freq data=flag_fdr_ge2_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge3_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;

proc freq data=flag_fdr_ge5_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge10_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;

proc print data=c1;
proc print data=c2;
run;


data flag_fdr_ge5_2;
  set flag_fdr_ge5;
  rename flag_fdr05=flag_binomial_fdr05;
run;

proc sort data=flag_fdr_ge5_m1;
  by comparison site_type chr region_num;
proc sort data=flag_fdr_ge5_2;
  by comparison site_type chr region_num;
run;

data dmr_lin_vs_binom5;
 merge flag_fdr_ge5_m1 flag_fdr_ge5_2;
  by comparison site_type chr region_num;
run;

proc freq data=dmr_lin_vs_binom5 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05*flag_binomial_fdr05 / out=c1;
run;

proc print data=c1;
run;



data mixed1;
  set wgbslocA.mixedacc_dar_anova_all;
   where effect="group";
run;


proc sort data=mixed1;
  by comparison site_type chr region_num;
proc sort data=dmr_methdiff;
  by comparison site_type chr region_num;
run;

data mixed2region;
  merge dmr_methdiff (in=in1) mixed1 (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge2 noprint;
  where flag_num_sites_ge2=1;
  by comparison site_type;
run;

proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge3 noprint;
  where flag_num_sites_ge3=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=mixed2region fdr out=binom_fdr_ge10 noprint;
  where flag_num_sites_ge10=1;
  by comparison site_type;
run;

data flag_fdr_ge2_m1;
set binom_fdr_ge2;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge3_m1;
set binom_fdr_ge3;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge5_m1;
set binom_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


data flag_fdr_ge10_m1;
set binom_fdr_ge10;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


proc freq data=flag_fdr_ge2_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;proc print data=c1;
proc print data=c2;
run;



proc freq data=flag_fdr_ge3_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;proc print data=c1;
proc print data=c2;
run;



proc freq data=flag_fdr_ge5_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;
proc print data=c1;
proc print data=c2;
run;


proc freq data=flag_fdr_ge10_m1 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05 / out=c1;
  tables comparison*site_type*flag_methdiff20*flag_fdr05 / out=c2;
run;

proc print data=c1;
proc print data=c2;
run;




data flag_fdr_ge5_2;
  set flag_fdr_ge5;
  rename flag_fdr05=flag_binomial_fdr05;
run;

proc sort data=flag_fdr_ge5_m1;
  by comparison site_type chr region_num;
proc sort data=flag_fdr_ge5_2;
  by comparison site_type chr region_num;
run;

data dmr_lin_vs_binom5;
 merge flag_fdr_ge5_m1 flag_fdr_ge5_2;
  by comparison site_type chr region_num;
run;

proc freq data=dmr_lin_vs_binom5 noprint;
  tables comparison*site_type*flag_methdiff10*flag_fdr05*flag_binomial_fdr05 / out=c1;
run;

proc print data=c1;
run;






%macro checkDAR(numSites);


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_GC_minSites_&numSites..DAR.output.txt" out=met_01 dbms=tab replace; guessingrows=50000; getnames=no; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_GC_minSites_&numSites..DAR.output.txt" out=met_1 dbms=tab replace; guessingrows=50000; getnames=no;run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0cGy_1cGy_GC_minSites_&numSites..DAR.results2.txt" out=mk_1 dbms=tab replace; guessingrows=50000; run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0cGy_01cGy_GC_minSites_&numSites..DAR.results2.txt" out=mk_01 dbms=tab replace; guessingrows=50000; run;

data met_01_2; set met_01; var9_2=var9*1; drop var9; rename var9_2=var9; run;
data met_1_2; set met_01; var9_2=var9*1; drop var9; rename var9_2=var9; run;


data metilene_all;
   set met_01_2 (in=in1) met_1_2 (in=in2) ;
   length comparison $12.;
   length site_type $4.;
   site_type = "GC";
   if in1 then comparison="01Gy_0Gy";
   if in2 then comparison="1Gy_0Gy";
   chr2 = compress(VAR1);
   rename chr2=chr VAR2=start VAR3=end VAR4=metilene_q VAR5=metilene_meth_diff VAR6=num_sites VAR7=metilene_MWU_P;
run;




 proc sort data=flag_fdr_ge&numSites.;
  by comparison site_type chr start end;
 proc sort data=metilene_all;
  by comparison site_type chr  start end;
run;


data met_w_diff;
  merge flag_fdr_ge&numSites. (in=in1) metilene_all (in=in2);
  by comparison site_type chr start end;
  if in1 ;
  rename fdr_p = binomial_fdr_p;
run;

proc multtest inpvalues(metilene_MWU_P)=met_w_diff fdr out=met_fdr noprint;
  by comparison site_type;
run;

data met_fdr2;
  set met_fdr;
if fdr_p = . then flag_met_fdr05=.;
else if fdr_p < 0.05 then flag_met_fdr05=1;
else  flag_met_fdr05=0;

if metilene_q = . then flag_met_q05=.;
else if metilene_q < 0.05 then flag_met_q05=1;
else  flag_met_q05=0;

if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;

proc freq data=met_fdr2 noprint;
  table comparison*site_type*flaG_methdiff10*flag_met_q05 / out=c1;
  table comparison*site_type*flaG_methdiff20*flag_met_q05 / out=c2;
  table comparison*site_type*flaG_methdiff10*flag_met_q05*flag_fdr05 / out=c3;
run;

proc print data=c1;
proc print data=c2;
proc print data=c3;
run;










data mk_all;
   set mk_01 (in=in1) mk_1 (in=in2) ;
   length comparison $12.;
   length site_type $4.;
   if in1 then comparison="01Gy_0Gy";
   if in2 then comparison="1Gy_0Gy";
   site_type="GC";
   start = start -1;
run;


 


 proc sort data=flag_fdr_ge&numSites.;
  by comparison site_type chr start end;
 proc sort data=mk_all;
  by comparison site_type chr  start end;
run;


data mk_w_diff;
  merge flag_fdr_ge&numSites. (in=in1) mk_all (in=in2);
  by comparison site_type chr start end;
  if in1 ;
  rename fdr_p = binomial_fdr_p;
run;

proc multtest inpvalues(pvalue)=mk_w_diff fdr out=met_fdr noprint;
  by comparison site_type;
run;

data met_fdr2;
  set met_fdr;
if fdr_p = . then flag_mk_fdr05=.;
else if fdr_p < 0.05 then flag_mk_fdr05=1;
else  flag_mk_fdr05=0;

if qvalue = . then flag_mk_q05=.;
else if qvalue < 0.05 then flag_mk_q05=1;
else  flag_mk_q05=0;

if  abs(methyl_diff) >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if  abs(methyl_diff) >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;

proc freq data=met_fdr2 noprint;
  table comparison*site_type*flaG_methdiff10*flag_mk_q05 / out=c1;
  table comparison*site_type*flaG_methdiff20*flag_mk_q05 / out=c2;
  table comparison*site_type*flaG_methdiff10*flag_mk_q05*flag_fdr05 / out=c3;
run;

proc print data=c1;
proc print data=c2;
proc print data=c3;
run;

%mend;

%checkDAR(2);
%checkDAR(3);
%checkDAR(5);
%checkDAR(10);
