/* Finalize DMR and DAR criteria */


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
    var flag_FET_FDR05 flag_FET_FDR05_methdiff_10 flag_FET_FDR05_methdiff_20 flag_meth_diff_20perc flag_meth_diff_10perc methyl_diff;
    output out=dmc_per_dmr sum= min(start_pos)=region_start max(stop_pos)=region_stop mean(methyl_diff)=mean_methyl_diff;
run;

data dmc_per_dmr2;
  set dmc_per_dmr;
  if flag_FET_FDR05_methdiff_10 >= 2 then flag_DMC10_ge2=1; else flag_DMC10_ge2=0;
  if flag_FET_FDR05_methdiff_10 >= 3 then flag_DMC10_ge3=1; else flag_DMC10_ge3=0;
  region_length = region_stop - region_start;
  if region_length >= 150 then flag_length_ge_150bp=1; else flag_length_ge_150bp=0;
  if region_length >= 100 then flag_length_ge_100bp=1; else flag_length_ge_100bp=0;
run;

proc sort data=dmc_per_dmr2;
  by comparison site_type region_num ; 
proc means data=dmc_per_dmr2 noprint;
  where  flag_num_sites_ge5 =1;
by comparison site_type  ; 
  var region_length;
  output out=region_length_check mean=mean stddev=stddev min=min q1=q1 median=median q3=q3 max=max;
run;


proc freq data=dmc_per_dmr2 noprint;
  where flag_num_sites_ge5=1;
  tables comparison*site_type*flag_length_ge_150bp / out=check1;
  tables comparison*site_type*flag_length_ge_100bp / out=check2;
run;

proc print data=check1;
proc print data=check2;
run;


/*


               site_    flag_length_
 comparison    type       ge_100bp      COUNT

 0Gy_vs_01G     CG            0          2762
 0Gy_vs_01G     CG            1          4247

 0Gy_vs_01G     CHG           0          4164
 0Gy_vs_01G     CHG           1          8788

 0Gy_vs_01G     CHH           0         43446
 0Gy_vs_01G     CHH           1         38612

 0Gy_vs_1Gy     CG            0          2116
 0Gy_vs_1Gy     CG            1          3305

 0Gy_vs_1Gy     CHG           0          3381
 0Gy_vs_1Gy     CHG           1          5984

 0Gy_vs_1Gy     CHH           0         48293
 0Gy_vs_1Gy     CHH           1         29763



               site_    flag_length_
 comparison    type       ge_150bp      COUNT

 0Gy_vs_01G     CG            0          4652
 0Gy_vs_01G     CG            1          2357

 0Gy_vs_01G     CHG           0          7630
 0Gy_vs_01G     CHG           1          5322

 0Gy_vs_01G     CHH           0         61658
 0Gy_vs_01G     CHH           1         20400

 0Gy_vs_1Gy     CG            0          3644
 0Gy_vs_1Gy     CG            1          1777

 0Gy_vs_1Gy     CHG           0          6082
 0Gy_vs_1Gy     CHG           1          3283

 0Gy_vs_1Gy     CHH           0         64073
 0Gy_vs_1Gy     CHH           1         13983


*/




/* Binomial test */

data binomial;
  set wgbslocA.bin2_dmr_anova_1 wgbslocA.bin2_dmr_anova_2 wgbslocA.bin2_dmr_anova_3 wgbslocA.bin2_dmr_anova_4 wgbslocA.bin2_dmr_anova_5;
run;


proc sort data=binomial;
  by comparison site_type chr region_num;
proc sort data=dmc_per_dmr2;
  by comparison site_type chr region_num;
run;

data binomial2region;
  merge dmc_per_dmr2 (in=in1) binomial (in=in2);
  by comparison site_type chr region_num;
  if in1;
run;

proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5_150bp noprint;
  where flag_num_sites_ge5=1 and flag_length_ge_150bp=1;
  by comparison site_type;
run;
proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5_100bp noprint;
  where flag_num_sites_ge5=1 and flag_length_ge_100bp=1;
  by comparison site_type;
run;


data flag_fdr_ge5;
set binom_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if mean_methyl_diff >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if mean_methyl_diff >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;






data flag_fdr_ge5_100bp;
set binom_fdr_ge5_100bp;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if mean_methyl_diff >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if mean_methyl_diff >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;



data flag_fdr_ge5_150bp;
set binom_fdr_ge5_150bp;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if mean_methyl_diff >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if mean_methyl_diff >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


proc freq data=flag_fdr_ge5 noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_150bp / out=c1;
run;

proc freq data=flag_fdr_ge5 noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_100bp / out=c1a;
run;



proc freq data=flag_fdr_ge5_100bp noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_100bp / out=c2;
run;

proc freq data=flag_fdr_ge5_150bp noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_150bp / out=c3;
run;


proc print data=c1a;
proc print data=c1;
proc print data=c2;
proc print data=c3;
run;


/* Binomial DMRs with >10% diff:

               site_    flag_length_
 comparison    type       ge_100bp      COUNT    PERCENT

 0Gy_vs_01G     CG            0           623     1.6115
 0Gy_vs_01G     CG            1           874     2.2607
 0Gy_vs_01G     CHG           0           427     1.1045
 0Gy_vs_01G     CHG           1           706     1.8262
 0Gy_vs_01G     CHH           0          4522    11.6968
 0Gy_vs_01G     CHH           1          2926     7.5685
 0Gy_vs_1Gy     CG            0            19     0.0491
 0Gy_vs_1Gy     CG            1            69     0.1785
 0Gy_vs_1Gy     CHG           0            78     0.2018
 0Gy_vs_1Gy     CHG           1           495     1.2804
 0Gy_vs_1Gy     CHH           0         16175    41.8391
 0Gy_vs_1Gy     CHH           1         11746    30.3828

                  The SAS System


               site_    flag_length_
 comparison    type       ge_150bp      COUNT    PERCENT

 0Gy_vs_01G     CG            0           986     2.5504
 0Gy_vs_01G     CG            1           511     1.3218
 0Gy_vs_01G     CHG           0           717     1.8546
 0Gy_vs_01G     CHG           1           416     1.0760
 0Gy_vs_01G     CHH           0          5999    15.5173
 0Gy_vs_01G     CHH           1          1449     3.7481
 0Gy_vs_1Gy     CG            0            39     0.1009
 0Gy_vs_1Gy     CG            1            49     0.1267
 0Gy_vs_1Gy     CHG           0           190     0.4915
 0Gy_vs_1Gy     CHG           1           383     0.9907
 0Gy_vs_1Gy     CHH           0         22090    57.1392
 0Gy_vs_1Gy     CHH           1          5831    15.0828

                  The SAS System



 Binomial DMRs with >10% diff and min 100bp in length:


              site_    flag_length_
comparison    type       ge_100bp      COUNT    PERCENT

0Gy_vs_01G     CG            1           965     4.9827
0Gy_vs_01G     CHG           1           797     4.1152
0Gy_vs_01G     CHH           1          3190    16.4713
0Gy_vs_1Gy     CG            1           560     2.8915
0Gy_vs_1Gy     CHG           1          1497     7.7296
0Gy_vs_1Gy     CHH           1         12358    63.8096




 Binomial DMRs with >10% diff and min 150bp in length:

               site_    flag_length_
 comparison    type       ge_150bp      COUNT    PERCENT

 0Gy_vs_01G     CG            1           576     5.5390
 0Gy_vs_01G     CHG           1           476     4.5774
 0Gy_vs_01G     CHH           1          1644    15.8092
 0Gy_vs_1Gy     CG            1           432     4.1542
 0Gy_vs_1Gy     CHG           1          1111    10.6837
 0Gy_vs_1Gy     CHH           1          6160    59.2365



*/


data wgbsA.results_by_dmr_5sites;
   set flag_fdr_ge5;
run;

data wgbsA.results_by_dmr_5sites_100bp;
   set flag_fdr_ge5_100bp;
run;


data wgbsA.results_by_dmr_5sites_150bp;
   set flag_fdr_ge5_150bp;
run;




/* DARs */


proc datasets lib=work kill noprint;
run;
quit;




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

proc sort data=meth2region2;
  by comparison site_type chr region_num;
proc sort data=region_flags;
  by comparison site_type chr region_num;
run;

data region_w_dac;
  merge region_flags (in=in1) meth2region2 (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


proc sort data=region_w_dac ;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
proc means data=region_w_dac noprint;
    by comparison site_type chr region_num flag_num_sites_ge2 flag_num_sites_ge3 flag_num_sites_ge5 flag_num_sites_ge10;
    var flag_CMH_FDR05 flag_CHM_FDR05_10perc flag_CHM_FDR05_20perc flag_meth_diff_20perc flag_meth_diff_10perc methyl_diff_TRT_CTL;
    output out=dac_per_dar sum= min(start_pos)=region_start max(stop_pos)=region_stop mean(methyl_diff_TRT_CTL)=mean_methyl_diff;
run;

data dac_per_dar2;
  set dac_per_dar;
  if flag_CHM_FDR05_10perc >= 2 then flag_DAC10_ge2=1; else flag_DAC10_ge2=0;
  if flag_CHM_FDR05_10perc >= 3 then flag_DAC10_ge3=1; else flag_DAC10_ge3=0;
  if flag_CHM_FDR05_20perc >= 2 then flag_DAC20_ge2=1; else flag_DAC20_ge2=0;
  if flag_CHM_FDR05_20perc >= 3 then flag_DAC20_ge3=1; else flag_DAC20_ge3=0;
  region_length = region_stop - region_start;
  if region_length >= 150 then flag_length_ge_150bp=1; else flag_length_ge_150bp=0;
  if region_length >= 100 then flag_length_ge_100bp=1; else flag_length_ge_100bp=0;
run;



proc sort data=dac_per_dar2;
  by comparison site_type region_num ; 
proc means data=dac_per_dar2 noprint;
  where  flag_num_sites_ge5 =1;
by comparison site_type  ; 
  var region_length;
  output out=region_length_check_dar mean=mean stddev=stddev min=min q1=q1 median=median q3=q3 max=max;
run;


proc freq data=dac_per_dar2 noprint;
  where flag_num_sites_ge5=1;
  tables comparison*site_type*flag_length_ge_150bp / out=check1;
  tables comparison*site_type*flag_length_ge_100bp/ out=check2;
run;

proc print data=check1;
proc print data=check2;
run;


/*

                site_    flag_length_
  comparison    type       ge_150bp      COUNT    PERCENT

   01Gy_0Gy      GC            0          2202    35.3451
   01Gy_0Gy      GC            1          2368    38.0096

   1Gy_0Gy       GC            0           897    14.3981
   1Gy_0Gy       GC            1           763    12.2472


               site_    flag_length_
 comparison    type       ge_100bp      COUNT    PERCENT

  01Gy_0Gy      GC            0           993    15.9390
  01Gy_0Gy      GC            1          3577    57.4157

  1Gy_0Gy       GC            0           432     6.9342
  1Gy_0Gy       GC            1          1228    19.7111



*/



/* Binomial test */

data binomial;
  set wgbslocA.bin2_dar_anova2_all;
  where effect = "group*units";
run;


proc sort data=binomial;
  by comparison site_type chr region_num;
proc sort data=dac_per_dar2;
  by comparison site_type chr region_num;
run;

data binomial2region;
  merge dac_per_dar2 (in=in1) binomial (in=in2);
  by comparison site_type chr region_num;
  if in1 ;
run;





proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5 noprint;
  where flag_num_sites_ge5=1;
  by comparison site_type;
run;


proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5_150bp noprint;
  where flag_num_sites_ge5=1 and flag_length_ge_150bp=1;
  by comparison site_type;
run;
proc multtest inpvalues(probf)=binomial2region fdr out=binom_fdr_ge5_100bp noprint;
  where flag_num_sites_ge5=1 and flag_length_ge_100bp=1;
  by comparison site_type;
run;


data flag_fdr_ge5;
set binom_fdr_ge5;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if mean_methyl_diff >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if mean_methyl_diff >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;






data flag_fdr_ge5_100bp;
set binom_fdr_ge5_100bp;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if mean_methyl_diff >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if mean_methyl_diff >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;



data flag_fdr_ge5_150bp;
set binom_fdr_ge5_150bp;
if fdr_p = . then flag_fdr05=.;
else if fdr_p < 0.05 then flag_fdr05=1;
else  flag_fdr05=0;
if mean_methyl_diff >= 0.1 then flaG_methdiff10=1; else flag_methdiff10=0;
if mean_methyl_diff >= 0.2 then flaG_methdiff20=1; else flag_methdiff20=0;
run;


proc freq data=flag_fdr_ge5 noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_150bp / out=c1;
run;

proc freq data=flag_fdr_ge5 noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_100bp / out=c1a;
run;



proc freq data=flag_fdr_ge5_100bp noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_100bp / out=c2;
run;

proc freq data=flag_fdr_ge5_150bp noprint;
   where flag_methdiff10=1 and flag_fdr05=1; 
  tables comparison*site_type*flag_length_ge_150bp / out=c3;
run;


proc print data=c1a;
proc print data=c1;
proc print data=c2;
proc print data=c3;
run;


/*   

Binomial DMRs with >10% diff:


                site_    flag_length_
  comparison    type       ge_100bp      COUNT

   01Gy_0Gy      GC            0           628
   01Gy_0Gy      GC            1          2850
   1Gy_0Gy       GC            0           153
   1Gy_0Gy       GC            1           619



               site_    flag_length_
 comparison    type       ge_150bp      COUNT    PERCENT

  01Gy_0Gy      GC            0          1513    35.6000
  01Gy_0Gy      GC            1          1965    46.2353
  1Gy_0Gy       GC            0           359     8.4471
  1Gy_0Gy       GC            1           413     9.7176



Binomial DMRs with >10% diff and min length of 100bp:

                site_    flag_length_
  comparison    type       ge_100bp      COUNT

   01Gy_0Gy      GC            1          2856
   1Gy_0Gy       GC            1           625





Binomial DMRs with >10% diff and min length of 150bp:



                 site_    flag_length_
   comparison    type       ge_150bp      COUNT    PERCENT

    01Gy_0Gy      GC            1          1971    82.3652
    1Gy_0Gy       GC            1           422    17.6348


*/




data wgbsA.results_by_dar_5sites;
   set flag_fdr_ge5;
run;

data wgbsA.results_by_dar_5sites_100bp;
   set flag_fdr_ge5_100bp;
run;


data wgbsA.results_by_dar_5sites_150bp;
   set flag_fdr_ge5_150bp;
run;




/* Export DARs and DMRs for HOMER annotation */

data dmr_bed;
  set wgbsA.results_by_dmr_5sites;
  length regionID $100.;
  regionID=catx("|",comparison,site_type,chr,region_num);
  keep chr region_start region_stop regionID;
run;


data dar_bed;
  set wgbsA.results_by_dar_5sites;
  length regionID $100.;
  regionID=catx("|",comparison,site_type,chr,region_num);
  keep chr region_start region_stop regionID;
run;

proc export data=dmr_bed
     outfile="/TB14/TB14/sandbox/dtra_sandbox/DMRs_min_5_sites_for_HOMER.bed"
     dbms=tab replace;
     putnames=no;
run;

proc export data=dar_bed
     outfile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER.bed"
     dbms=tab replace;
     putnames=no;
run;


