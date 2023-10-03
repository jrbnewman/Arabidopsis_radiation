libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Compare DMC and DAC: FET/CMH vs "credible difference" vs binomial */

/* Define a credible difference in methylation as >10%, >20%, >30%, >50% */

%macro prepData(trt);

data counts;
   set wgbslocA.methylation_data_cg_chg_chh;
   if chr = "Mt" then delete;
   if chr = "Pt" then delete;
   where units="0U" and treatment="&trt.";
run;  

proc sort data=counts;
   by chr start_pos stop_pos site_type treatment units;
proc means data=counts noprint;
   by chr start_pos stop_pos site_type treatment units;
   var total_C methyl_C perc_methyl;
   output out=counts2 sum(total_C)=total_C_&trt. sum(methyl_C)=methyl_C_&trt. mean(perc_methyl)=perc_methyl_&trt.;
run;

data cov;
   set wgbslocA.flag_coverage_10x_cg_chg_chh;
   where flag_&trt._coverage_ge_10x=1;
   keep site_type chr start_pos stop_pos;
run;

proc sort data=cov;
  by site_type chr start_pos stop_pos;
proc sort data=counts2;
  by site_type chr start_pos stop_pos;
run;

data counts_&trt.;
   merge counts2 (in=in1) cov (in=in2);
   by site_type chr start_pos stop_pos;
   if in1 and in2;
run;

%mend;

%prepData(0Gy);
%prepData(01Gy);
%prepData(1Gy);

proc sort data=counts_0Gy;
  by site_type chr start_pos stop_pos;
proc sort data=counts_01Gy;
  by site_type chr start_pos stop_pos;
proc sort data=counts_1Gy;
  by site_type chr start_pos stop_pos;
run;

data meth_0Gy_01Gy;
  merge counts_0Gy (in=in1) counts_01Gy (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

data meth_0Gy_1Gy;
  merge counts_0Gy (in=in1) counts_1Gy (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

data methdiff_0Gy_01Gy;
   set meth_0Gy_01Gy;
   methdiff_01Gy_0Gy = perc_methyl_01gy - perc_methyl_0gy ;
   if abs(methdiff_01Gy_0Gy) > 0.1 then flag_methdiff_gt_10perc = 1;
   else flag_methdiff_gt_10perc = 0;
   if abs(methdiff_01Gy_0Gy) > 0.2 then flag_methdiff_gt_20perc = 1;
   else flag_methdiff_gt_20perc = 0;
   if abs(methdiff_01Gy_0Gy) > 0.3 then flag_methdiff_gt_30perc = 1;
   else flag_methdiff_gt_30perc = 0;
   if abs(methdiff_01Gy_0Gy) > 0.5 then flag_methdiff_gt_50perc = 1;
   else flag_methdiff_gt_50perc = 0;
run;

data methdiff_0Gy_1Gy;
   set meth_0Gy_1Gy;
   methdiff_1Gy_0Gy = perc_methyl_1gy - perc_methyl_0gy ;
   if abs(methdiff_1Gy_0Gy) > 0.1 then flag_methdiff_gt_10perc = 1;
   else flag_methdiff_gt_10perc = 0;
   if abs(methdiff_1Gy_0Gy) > 0.2 then flag_methdiff_gt_20perc = 1;
   else flag_methdiff_gt_20perc = 0;
   if abs(methdiff_1Gy_0Gy) > 0.3 then flag_methdiff_gt_30perc = 1;
   else flag_methdiff_gt_30perc = 0;
   if abs(methdiff_1Gy_0Gy) > 0.5 then flag_methdiff_gt_50perc = 1;
   else flag_methdiff_gt_50perc = 0;
run;

proc freq data=methdiff_0Gy_01Gy;
  tables  site_type site_type*(flag_methdiff_gt_10perc flag_methdiff_gt_20perc
          flag_methdiff_gt_30perc flag_methdiff_gt_50perc);
run;


/*

                                         Cumulative    Cumulative
   site_type    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   CG            3721316       14.72       3721316        14.72
   CHG           4067017       16.09       7788333        30.81
   CHH          17486471       69.19      25274804       100.00


            flag_methdiff_gt_10perc

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
  CG       |3242342 | 478974 |3721316
           |  12.83 |   1.90 |  14.72
           |  87.13 |  12.87 |
           |  14.55 |  16.05 |
  ---------+--------+--------+
  CHG      |3537665 | 529352 |4067017
           |  14.00 |   2.09 |  16.09
           |  86.98 |  13.02 |
           |  15.87 |  17.74 |
  ---------+--------+--------+
  CHH      |1.551E7 |1975471 |1.749E7
           |  61.37 |   7.82 |  69.19
           |  88.70 |  11.30 |
           |  69.58 |  66.21 |
  ---------+--------+--------+
  Total     2.229E7  2983797  2.527E7
              88.19    11.81   100.00

  site_type
            flag_methdiff_gt_20perc

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
  CG       |3628037 |  93279 |3721316
           |  14.35 |   0.37 |  14.72
           |  97.49 |   2.51 |
           |  14.70 |  15.84 |
  ---------+--------+--------+
  CHG      |3951517 | 115500 |4067017
           |  15.63 |   0.46 |  16.09
           |  97.16 |   2.84 |
           |  16.01 |  19.62 |
  ---------+--------+--------+
  CHH      |1.711E7 | 379957 |1.749E7
           |  67.68 |   1.50 |  69.19
           |  97.83 |   2.17 |
           |  69.30 |  64.54 |
  ---------+--------+--------+
  Total     2.469E7   588736  2.527E7
              97.67     2.33   100.00


 Table of site_type by flag_methdiff_gt_30pe

      site_type
                flag_methdiff_gt_30perc

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
      CG       |3698944 |  22372 |3721316
               |  14.63 |   0.09 |  14.72
               |  99.40 |   0.60 |
               |  14.72 |  15.52 |
      ---------+--------+--------+
      CHG      |4038430 |  28587 |4067017
               |  15.98 |   0.11 |  16.09
               |  99.30 |   0.70 |
               |  16.07 |  19.83 |
      ---------+--------+--------+
      CHH      |1.739E7 |  93201 |1.749E7
               |  68.82 |   0.37 |  69.19
               |  99.47 |   0.53 |
               |  69.21 |  64.65 |
      ---------+--------+--------+
      Total     2.513E7   144160  2.527E7
                  99.43     0.57   100.00


  Table of site_type by flag_methdiff_gt_50perc

       site_type
                 flag_methdiff_gt_50perc

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
       CG       |3720046 |   1270 |3721316
                |  14.72 |   0.01 |  14.72
                |  99.97 |   0.03 |
                |  14.72 |  14.57 |
       ---------+--------+--------+
       CHG      |4065490 |   1527 |4067017
                |  16.09 |   0.01 |  16.09
                |  99.96 |   0.04 |
                |  16.09 |  17.52 |
       ---------+--------+--------+
       CHH      |1.748E7 |   5920 |1.749E7
                |  69.16 |   0.02 |  69.19
                |  99.97 |   0.03 |
                |  69.19 |  67.91 |
       ---------+--------+--------+
       Total     2.527E7     8717  2.527E7
                   99.97     0.03   100.00
*/

proc freq data=methdiff_0Gy_1Gy;
  tables  site_type site_type*(flag_methdiff_gt_10perc flag_methdiff_gt_20perc
          flag_methdiff_gt_30perc flag_methdiff_gt_50perc);
run;



/*
                                       Cumulative    Cumulative
 site_type    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 CG            3309366       14.86       3309366        14.86
 CHG           3591816       16.13       6901182        30.98
 CHH          15373571       69.02      22274753       100.00


         Table of site_type by flag_methdiff_gt_10perc

              site_type
                        flag_methdiff_gt_10perc

              Frequency|
              Percent  |
              Row Pct  |
              Col Pct  |       0|       1|  Total
              ---------+--------+--------+
              CG       |2883115 | 426251 |3309366
                       |  12.94 |   1.91 |  14.86
                       |  87.12 |  12.88 |
                       |  14.72 |  15.82 |
              ---------+--------+--------+
              CHG      |3122165 | 469651 |3591816
                       |  14.02 |   2.11 |  16.13
                       |  86.92 |  13.08 |
                       |  15.95 |  17.43 |
              ---------+--------+--------+
              CHH      |1.358E7 |1798500 |1.537E7
                       |  60.94 |   8.07 |  69.02
                       |  88.30 |  11.70 |
                       |  69.33 |  66.75 |
              ---------+--------+--------+
              Total     1.958E7  2694402  2.227E7
                          87.90    12.10   100.00

 Table of site_type by flag_methdiff_gt_20perc

      site_type
                flag_methdiff_gt_20perc

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
      CG       |3222299 |  87067 |3309366
               |  14.47 |   0.39 |  14.86
               |  97.37 |   2.63 |
               |  14.85 |  15.02 |
      ---------+--------+--------+
      CHG      |3477931 | 113885 |3591816
               |  15.61 |   0.51 |  16.13
               |  96.83 |   3.17 |
               |  16.03 |  19.64 |
      ---------+--------+--------+
      CHH      |1.499E7 | 378910 |1.537E7
               |  67.32 |   1.70 |  69.02
               |  97.54 |   2.46 |
               |  69.12 |  65.34 |
      ---------+--------+--------+
      Total     2.169E7   579862  2.227E7
                  97.40     2.60   100.00

   Table of site_type by flag_methdiff_gt_30perc

        site_type
                  flag_methdiff_gt_30perc

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
        CG       |3288343 |  21023 |3309366
                 |  14.76 |   0.09 |  14.86
                 |  99.36 |   0.64 |
                 |  14.86 |  13.72 |
        ---------+--------+--------+
        CHG      |3560306 |  31510 |3591816
                 |  15.98 |   0.14 |  16.13
                 |  99.12 |   0.88 |
                 |  16.09 |  20.56 |
        ---------+--------+--------+
        CHH      |1.527E7 | 100723 |1.537E7
                 |  68.57 |   0.45 |  69.02
                 |  99.34 |   0.66 |
                 |  69.04 |  65.72 |
        ---------+--------+--------+
        Total     2.212E7   153256  2.227E7
                    99.31     0.69   100.00

  Table of site_type by flag_methdiff_gt_50perc

       site_type
                 flag_methdiff_gt_50perc

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
       CG       |3308037 |   1329 |3309366
                |  14.85 |   0.01 |  14.86
                |  99.96 |   0.04 |
                |  14.86 |  11.32 |
       ---------+--------+--------+
       CHG      |3589722 |   2094 |3591816
                |  16.12 |   0.01 |  16.13
                |  99.94 |   0.06 |
                |  16.12 |  17.84 |
       ---------+--------+--------+
       CHH      |1.537E7 |   8313 |1.537E7
                |  68.98 |   0.04 |  69.02
                |  99.95 |   0.05 |
                |  69.02 |  70.83 |
       ---------+--------+--------+
       Total     2.226E7    11736  2.227E7
                   99.95     0.05   100.00


*/

data flag_small_n_01;
  set methdiff_0gy_01gy;
  if methyl_C_0Gy < 5 or methyl_C_01Gy < 5 or 
     (total_C_0Gy - methyl_C_0Gy) < 5 or (total_C_01Gy - methyl_C_01Gy) < 5;
  keep site_type chr start_pos stop_pos;
run;

data flag_small_n_1;
  set methdiff_0gy_1gy;
  if methyl_C_0Gy < 5 or methyl_C_1Gy < 5 or 
     (total_C_0Gy - methyl_C_0Gy) < 5 or (total_C_1Gy - methyl_C_1Gy) < 5;
  keep site_type chr start_pos stop_pos;
run;




/* prep FET DMC results */

data fet_dmc_01;
  set wgbsA.results_by_dmc_0gy_01gy_: ;
  where P_PCHI ne . ;
  keep site_type chr start_pos stop_pos XP2_FISH P_PCHI;
run;

data fet_dmc_1;
  set wgbsA.results_by_dmc_0gy_1gy_: ;
  where P_PCHI ne . ;
  keep site_type chr start_pos stop_pos XP2_FISH P_PCHI;
run;

proc sort data=fet_dmc_01;
  by site_type chr start_pos stop_pos ;
proc sort data=fet_dmc_1;
  by site_type chr start_pos stop_pos ;
proc sort data=flag_small_n_01;
  by site_type chr start_pos stop_pos ;
proc sort data=flag_small_n_1;
  by site_type chr start_pos stop_pos ;
run;


data fet_dmc_01_2;
  merge fet_dmc_01 (in=in1) flag_small_n_01 (in=in2);
  by site_type chr start_pos stop_pos;
  if in2 then flag_use_FET=1; else flag_use_FET=0;
  if in1 then output;
run;

data fet_dmc_1_2;
  merge fet_dmc_1 (in=in1) flag_small_n_1 (in=in2);
  by site_type chr start_pos stop_pos;
  if in2 then flag_use_FET=1; else flag_use_FET=0;
  if in1 then output;
run;

proc freq data=fet_dmc_01_2;
  tables site_type*flag_use_FET;
proc freq data=fet_dmc_1_2;
  tables site_type*flag_use_FET;
run;


/*
    Table of site_type by flag_use_FET

  site_type     flag_use_FET

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
  CG       | 450707 |2674765 |3125472
           |   2.17 |  12.90 |  15.08
           |  14.42 |  85.58 |
           |  17.28 |  14.76 |
  ---------+--------+--------+
  CHG      | 684376 |2642894 |3327270
           |   3.30 |  12.75 |  16.05
           |  20.57 |  79.43 |
           |  26.23 |  14.58 |
  ---------+--------+--------+
  CHH      |1473805 | 1.28E7 |1.428E7
           |   7.11 |  61.76 |  68.87
           |  10.32 |  89.68 |
           |  56.49 |  70.65 |
  ---------+--------+--------+
  Total     2608888  1.812E7  2.073E7
              12.59    87.41   100.00



 Table of site_type by flag_use_FET

site_type     flag_use_FET

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
CG       | 417461 |2232736 |2650197
         |   2.41 |  12.90 |  15.32
         |  15.75 |  84.25 |
         |  16.19 |  15.16 |
---------+--------+--------+
CHG      | 637801 |2141530 |2779331
         |   3.69 |  12.38 |  16.06
         |  22.95 |  77.05 |
         |  24.73 |  14.54 |
---------+--------+--------+
CHH      |1523522 |1.035E7 |1.187E7
         |   8.80 |  59.82 |  68.62
         |  12.83 |  87.17 |
         |  59.08 |  70.29 |
---------+--------+--------+
Total     2578784  1.473E7   1.73E7
            14.90    85.10   100.00

*/


data fet_dmc_01_3;
  set fet_dmc_01_2;
  if flag_use_FET=1 then Pval = XP2_FISH;
  else Pval = P_PCHI;
run;

data fet_dmc_1_3;
  set fet_dmc_1_2;
  if flag_use_FET=1 then Pval = XP2_FISH;
  else Pval = P_PCHI;
run;



ods graphics / antialiasmax=26000000  obsmax=26000000;

proc sgplot data=fet_dmc_01_3;
  by site_type;
  scatter x=XP2_FISH y=P_PCHI;
run;

proc sgplot data=fet_dmc_1_3;
  by site_type;
  scatter x=XP2_FISH y=P_PCHI;
run;



proc sgplot data=fet_dmc_01_3;
  by site_type;
  scatter x=XP2_FISH y=Pval;
run;

proc sgplot data=fet_dmc_1_3;
  by site_type;
  scatter x=XP2_FISH y=Pval;
run;


proc multtest inpvalues(XP2_FISH)=fet_dmc_01_3 fdr out=fet_dmc_01_fdr(rename=(fdr_p=FET_fdr_p)) noprint;
  by site_type;
run;

proc multtest inpvalues(XP2_FISH)=fet_dmc_1_3 fdr out=fet_dmc_1_fdr(rename=(fdr_p=FET_fdr_p)) noprint;
  by site_type;
run;


proc multtest inpvalues(Pval)=fet_dmc_01_fdr fdr out=fet_dmc_01_fdr2(rename=(fdr_p=Pval_fdr_p)) noprint;
  by site_type;
run;

proc multtest inpvalues(Pval)=fet_dmc_1_fdr fdr out=fet_dmc_1_fdr2(rename=(fdr_p=Pval_fdr_p)) noprint;
  by site_type;
run;



proc sgplot data=fet_dmc_01_fdr2;
  by site_type;
  scatter x=FET_fdr_p y=Pval_fdr_p;
run;

proc sgplot data=fet_dmc_1_fdr2;
  by site_type;
  scatter x=FET_fdr_p y=Pval_fdr_p;
run;

/* Flag sig */

data fet_dmc_01_fdr3;
  set fet_dmc_01_fdr;
  where XP2_FISH ne .;
  if FET_fdr_p < 0.05 then flag_FET_FDR05=1; else flag_FET_FDR05=0;
  if FET_fdr_p < 0.01 then flag_FET_FDR01=1; else flag_FET_FDR01=0;
run;


data fet_dmc_1_fdr3;
  set fet_dmc_1_fdr;
  where XP2_FISH ne .;
  if FET_fdr_p < 0.05 then flag_FET_FDR05=1; else flag_FET_FDR05=0;
  if FET_fdr_p < 0.01 then flag_FET_FDR01=1; else flag_FET_FDR01=0;
run;

proc freq data=fet_dmc_01_fdr3;
  tables site_type * (flag_FET_FDR05 flag_FET_FDR01);
run;

proc freq data=fet_dmc_1_fdr3;
  tables site_type * (flag_FET_FDR05 flag_FET_FDR01);
run;


/*
    Table of site_type by flag_FET_FDR05

    site_type     flag_FET_FDR05

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
    CG       |3124842 |    630 |3125472
             |  15.07 |   0.00 |  15.08
             |  99.98 |   0.02 |
             |  15.08 |  50.08 |
    ---------+--------+--------+
    CHG      |3327106 |    164 |3327270
             |  16.05 |   0.00 |  16.05
             | 100.00 |   0.00 |
             |  16.05 |  13.04 |
    ---------+--------+--------+
    CHH      |1.428E7 |    464 |1.428E7
             |  68.87 |   0.00 |  68.87
             | 100.00 |   0.00 |
             |  68.87 |  36.88 |
    ---------+--------+--------+
    Total     2.073E7     1258  2.073E7
                99.99     0.01   100.00



 Table of site_type by flag_FET_FDR01

 site_type     flag_FET_FDR01

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
 CG       |3124986 |    486 |3125472
          |  15.07 |   0.00 |  15.08
          |  99.98 |   0.02 |
          |  15.08 |  51.48 |
 ---------+--------+--------+
 CHG      |3327141 |    129 |3327270
          |  16.05 |   0.00 |  16.05
          | 100.00 |   0.00 |
          |  16.05 |  13.67 |
 ---------+--------+--------+
 CHH      |1.428E7 |    329 |1.428E7
          |  68.87 |   0.00 |  68.87
          | 100.00 |   0.00 |
          |  68.87 |  34.85 |
 ---------+--------+--------+
 Total     2.073E7      944  2.073E7
            100.00     0.00   100.00

            The SAS System


  Table of site_type by flag_FET_FDR0

  site_type     flag_FET_FDR05

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
  CG       |2650083 |    114 |2650197
           |  15.31 |   0.00 |  15.32
           | 100.00 |   0.00 |
           |  15.32 |  14.77 |
  ---------+--------+--------+
  CHG      |2779328 |      3 |2779331
           |  16.06 |   0.00 |  16.06
           | 100.00 |   0.00 |
           |  16.06 |   0.39 |
  ---------+--------+--------+
  CHH      |1.187E7 |    655 |1.187E7
           |  68.62 |   0.00 |  68.62
           |  99.99 |   0.01 |
           |  68.62 |  84.84 |
  ---------+--------+--------+
  Total      1.73E7      772   1.73E7
             100.00     0.00   100.00


 Table of site_type by flag_FET_FDR01

 site_type     flag_FET_FDR01

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
 CG      it's been done |2650138 |     59 |2650197
          |  15.32 |   0.00 |  15.32
          | 100.00 |   0.00 |
          |  15.32 |  11.41 |
 ---------+--------+--------+
 CHG      |2779328 |      3 |2779331
          |  16.06 |   0.00 |  16.06
          | 100.00 |   0.00 |
          |  16.06 |   0.58 |
 ---------+--------+--------+
 CHH      |1.187E7 |    455 |1.187E7
          |  68.62 |   0.00 |  68.62
          | 100.00 |   0.00 |
          |  68.62 |  88.01 |
 ---------+--------+--------+
 Total      1.73E7      517   1.73E7
            100.00     0.00   100.00

*/

/* prep binomial DMC results */

data binomial_dmc_01;
  set wgbsA.binom_dmc_0gy_01gy_: ;
  if ProbChiSq = . then delete;
run;

data binomial_dmc_1;
  set wgbsA.binom_dmc_0gy_1gy_: ;
  if ProbChiSq = . then delete;
run;

proc sort data=binomial_dmc_01;
  by site_type;
proc multtest inpvalues(ProbChiSq)=binomial_dmc_01 fdr out=binomial_dmc_01_fdr(rename=(fdr_p=binomial_fdr_p)) noprint;
  by site_type;
run;

proc sort data=binomial_dmc_1;
  by site_type;
proc multtest inpvalues(ProbChiSq)=binomial_dmc_1 fdr out=binomial_dmc_1_fdr(rename=(fdr_p=binomial_fdr_p)) noprint;
  by site_type;
run;





data binomial_dmc_01_fdr2;
  set binomial_dmc_01_fdr;
  if binomial_fdr_p < 0.05 then flag_binomial_FDR05=1; else flag_binomial_FDR05=0;
  if binomial_fdr_p < 0.01 then flag_binomial_FDR01=1; else flag_binomial_FDR01=0;
run;


data binomial_dmc_1_fdr2;
  set binomial_dmc_1_fdr;
  if binomial_fdr_p < 0.05 then flag_binomial_FDR05=1; else flag_binomial_FDR05=0;
  if binomial_fdr_p < 0.01 then flag_binomial_FDR01=1; else flag_binomial_FDR01=0;
run;

proc freq data=binomial_dmc_01_fdr2;
  tables site_type * (flag_binomial_FDR05 flag_binomial_FDR01);
run;

proc freq data=binomial_dmc_1_fdr2;
  tables site_type * (flag_binomial_FDR05 flag_binomial_FDR01);
run;

/*
 Table of site_type by flag_binomial_FDR05

   site_type     flag_binomial_FDR05

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
   CG       |2104298 |    577 |2104875
            |  17.70 |   0.00 |  17.70
            |  99.97 |   0.03 |
            |  17.70 |  47.41 |
   ---------+--------+--------+
   CHG      |1966146 |    164 |1966310
            |  16.53 |   0.00 |  16.53
            |  99.99 |   0.01 |
            |  16.54 |  13.48 |
   ---------+--------+--------+
   CHH      |7820331 |    476 |7820807
            |  65.76 |   0.00 |  65.77
            |  99.99 |   0.01 |
            |  65.77 |  39.11 |
   ---------+--------+--------+
   Total     1.189E7     1217  1.189E7
               99.99     0.01   100.00


 Table of site_type by flag_binomial_FDR0

    site_type     flag_binomial_FDR01

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
    CG       |2104430 |    445 |2104875
             |  17.70 |   0.00 |  17.70
             |  99.98 |   0.02 |
             |  17.70 |  47.24 |
    ---------+--------+--------+
    CHG      |1966176 |    134 |1966310
             |  16.53 |   0.00 |  16.53
             |  99.99 |   0.01 |
             |  16.53 |  14.23 |
    ---------+--------+--------+
    CHH      |7820444 |    363 |7820807
             |  65.76 |   0.00 |  65.77
             | 100.00 |   0.00 |
             |  65.77 |  38.54 |
    ---------+--------+--------+
    Total     1.189E7      942  1.189E7
                99.99     0.01   100.00


   Table of site_type by flag_binomial_FDR05

      site_type     flag_binomial_FDR05

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
      CG       |1800124 |     89 |1800213
               |  18.26 |   0.00 |  18.26
               | 100.00 |   0.00 |
               |  18.26 |  22.42 |
      ---------+--------+--------+
      CHG      |1630725 |      3 |1630728
               |  16.54 |   0.00 |  16.54
               | 100.00 |   0.00 |
               |  16.54 |   0.76 |
      ---------+--------+--------+
      CHH      |6427470 |    305 |6427775
               |  65.20 |   0.00 |  65.20
               | 100.00 |   0.00 |
               |  65.20 |  76.83 |
      ---------+--------+--------+
      Total     9858319      397  9858716
                 100.00     0.00   100.00


 Table of site_type by flag_binomial_FDR01

    site_type     flag_binomial_FDR01

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
    CG       |1800160 |     53 |1800213
             |  18.26 |   0.00 |  18.26
             | 100.00 |   0.00 |
             |  18.26 |  17.79 |
    ---------+--------+--------+
    CHG      |1630725 |      3 |1630728
             |  16.54 |   0.00 |  16.54
             | 100.00 |   0.00 |
             |  16.54 |   1.01 |
    ---------+--------+--------+
    CHH      |6427533 |    242 |6427775
             |  65.20 |   0.00 |  65.20
             | 100.00 |   0.00 |
             |  65.20 |  81.21 |
    ---------+--------+--------+
    Total     9858418      298  9858716
               100.00     0.00   100.00


*/



/* prep binary logit DMC results */

data logit_dmc_01;
  set wgbsA.logit_dmc_0gy_01gy_: ;
  if ProbChiSq = . then delete;
run;

data logit_dmc_1;
  set wgbsA.logit_dmc_0gy_1gy_: ;
  if ProbChiSq = . then delete;
run;

proc sort data=logit_dmc_01;
  by site_type;
proc multtest inpvalues(ProbChiSq)=logit_dmc_01 fdr out=logit_dmc_01_fdr(rename=(fdr_p=logit_fdr_p)) noprint;
  by site_type;
run;

proc sort data=logit_dmc_1;
  by site_type;
proc multtest inpvalues(ProbChiSq)=logit_dmc_1 fdr out=logit_dmc_1_fdr(rename=(fdr_p=logit_fdr_p)) noprint;
  by site_type;
run;





data logit_dmc_01_fdr2;
  set logit_dmc_01_fdr;
  if logit_fdr_p < 0.05 then flag_logit_FDR05=1; else flag_logit_FDR05=0;
  if logit_fdr_p < 0.01 then flag_logit_FDR01=1; else flag_logit_FDR01=0;
run;


data logit_dmc_1_fdr2;
  set logit_dmc_1_fdr;
  if logit_fdr_p < 0.05 then flag_logit_FDR05=1; else flag_logit_FDR05=0;
  if logit_fdr_p < 0.01 then flag_logit_FDR01=1; else flag_logit_FDR01=0;
run;

proc freq data=logit_dmc_01_fdr2;
  tables site_type * (flag_logit_FDR05 flag_logit_FDR01);
run;

proc freq data=logit_dmc_1_fdr2;
  tables site_type * (flag_logit_FDR05 flag_logit_FDR01);
run;

/*
 Table of site_type by flag_logit_FDR05

 site_type     flag_logit_FDR05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
 CG       |2104298 |    577 |2104875
          |  17.70 |   0.00 |  17.70
          |  99.97 |   0.03 |
          |  17.70 |  47.41 |
 ---------+--------+--------+
 CHG      |1966146 |    164 |1966310
          |  16.53 |   0.00 |  16.53
          |  99.99 |   0.01 |
          |  16.54 |  13.48 |
 ---------+--------+--------+
 CHH      |7820331 |    476 |7820807
          |  65.76 |   0.00 |  65.77
          |  99.99 |   0.01 |
          |  65.77 |  39.11 |
 ---------+--------+--------+
 Total     1.189E7     1217  1.189E7
             99.99     0.01   100.00


 Table of site_type by flag_logit_FDR01

  site_type     flag_logit_FDR01

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
  CG       |2104430 |    445 |2104875
           |  17.70 |   0.00 |  17.70
           |  99.98 |   0.02 |
           |  17.70 |  47.24 |
  ---------+--------+--------+
  CHG      |1966176 |    134 |1966310
           |  16.53 |   0.00 |  16.53
           |  99.99 |   0.01 |
           |  16.53 |  14.23 |
  ---------+--------+--------+
  CHH      |7820444 |    363 |7820807
           |  65.76 |   0.00 |  65.77
           | 100.00 |   0.00 |
           |  65.77 |  38.54 |
  ---------+--------+--------+
  Total     1.189E7      942  1.189E7
              99.99     0.01   100.00


Table of site_type by flag_logit_FDR05

 site_type     flag_logit_FDR05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
 CG       |1800124 |     89 |1800213
          |  18.26 |   0.00 |  18.26
          | 100.00 |   0.00 |
          |  18.26 |  22.42 |
 ---------+--------+--------+
 CHG      |1630725 |      3 |1630728
          |  16.54 |   0.00 |  16.54
          | 100.00 |   0.00 |
          |  16.54 |   0.76 |
 ---------+--------+--------+
 CHH      |6427470 |    305 |6427775
          |  65.20 |   0.00 |  65.20
          | 100.00 |   0.00 |
          |  65.20 |  76.83 |
 ---------+--------+--------+
 Total     9858319      397  9858716
            100.00     0.00   100.00


 Table of site_type by flag_logit_FDR01

  site_type     flag_logit_FDR01

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
  CG       |1800160 |     53 |1800213
           |  18.26 |   0.00 |  18.26
           | 100.00 |   0.00 |
           |  18.26 |  17.79 |
  ---------+--------+--------+
  CHG      |1630725 |      3 |1630728
           |  16.54 |   0.00 |  16.54
           | 100.00 |   0.00 |
           |  16.54 |   1.01 |
  ---------+--------+--------+
  CHH      |6427533 |    242 |6427775
           |  65.20 |   0.00 |  65.20
           | 100.00 |   0.00 |
           |  65.20 |  81.21 |
  ---------+--------+--------+
  Total     9858418      298  9858716
             100.00     0.00   100.00




*/



data credible_01_2merge;
  set methdiff_0Gy_01Gy;
  keep chr start_pos stop_pos site_type perc_methyl_0gy perc_methyl_01gy methdiff_01Gy_0Gy flag_methdiff_: ;
run;

data fet_01_2merge;
  set fet_dmc_01_fdr3;
  keep chr start_pos stop_pos site_type XP2_FISH FET_fdr_p flag_FET_fdr05 flag_FET_fdr01;
  rename XP2_FISH=FET_P;
run;

data binomial_01_2merge;
  set binomial_dmc_01_fdr2;
  keep chr start_pos stop_pos site_type ProbChiSq binomial_fdr_p flag_binomial_fdr05 flag_binomial_fdr01;
  rename ProbChiSq=binomial_P;
run;


data logit_01_2merge;
  set logit_dmc_01_fdr2;
  keep chr start_pos stop_pos site_type ProbChiSq logit_fdr_p flag_logit_fdr05 flag_logit_fdr01;
  rename ProbChiSq=logit_P;
run;

proc sort data=credible_01_2merge;
  by site_type chr start_pos stop_pos;
proc sort data=fet_01_2merge;
  by site_type chr start_pos stop_pos;
proc sort data=binomial_01_2merge;
  by site_type chr start_pos stop_pos;
proc sort data=logit_01_2merge;
  by site_type chr start_pos stop_pos;
run;

data all_dmc_tests_01;
  merge credible_01_2merge (in=in1) fet_01_2merge (in=in2) binomial_01_2merge (in=in3) logit_01_2merge (in=in4);
  by site_type chr start_pos stop_pos; 
  if in1 then flag_in_credible=1; else flag_in_credible=0;
  if in2 then flag_in_fet=1; else flag_in_fet=0;
  if in3 then flag_in_binomial=1; else flag_in_binomial=0;
  if in4 then flag_in_logit=1; else flag_in_logit=0;
run;







data credible_1_2merge;
  set methdiff_0Gy_1Gy;
  keep chr start_pos stop_pos site_type perc_methyl_0gy perc_methyl_1gy methdiff_1Gy_0Gy flag_methdiff_: ;
run;

data fet_1_2merge;
  set fet_dmc_1_fdr3;
  keep chr start_pos stop_pos site_type XP2_FISH FET_fdr_p flag_FET_fdr05 flag_FET_fdr01;
  rename XP2_FISH=FET_P;
run;

data binomial_1_2merge;
  set binomial_dmc_1_fdr2;
  keep chr start_pos stop_pos site_type ProbChiSq binomial_fdr_p flag_binomial_fdr05 flag_binomial_fdr01;
  rename ProbChiSq=binomial_P;
run;


data logit_1_2merge;
  set logit_dmc_1_fdr2;
  keep chr start_pos stop_pos site_type ProbChiSq logit_fdr_p flag_logit_fdr05 flag_logit_fdr01;
  rename ProbChiSq=logit_P;
run;

proc sort data=credible_1_2merge;
  by site_type chr start_pos stop_pos;
proc sort data=fet_1_2merge;
  by site_type chr start_pos stop_pos;
proc sort data=binomial_1_2merge;
  by site_type chr start_pos stop_pos;
proc sort data=logit_1_2merge;
  by site_type chr start_pos stop_pos;
run;

data all_dmc_tests_1;
  merge credible_1_2merge (in=in1) fet_1_2merge (in=in2) binomial_1_2merge (in=in3) logit_1_2merge (in=in4);
  by site_type chr start_pos stop_pos; 
  if in1 then flag_in_credible=1; else flag_in_credible=0;
  if in2 then flag_in_fet=1; else flag_in_fet=0;
  if in3 then flag_in_binomial=1; else flag_in_binomial=0;
  if in4 then flag_in_logit=1; else flag_in_logit=0;
run;






proc freq data=all_dmc_tests_01 noprint;
  tables  site_type * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit / out=dmc01_meth0_fdr05 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit / out=dmc01_meth10_fdr05 ;
  tables  site_type * flag_methdiff_gt_20perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc01_meth20_fdr05 ;
  tables  site_type * flag_methdiff_gt_30perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc01_meth30_fdr05 ;
  tables  site_type * flag_methdiff_gt_50perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc01_meth50_fdr05 ;

run;



proc freq data=all_dmc_tests_1 noprint;
  tables  site_type * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit / out=dmc1_meth0_fdr05 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc1_meth10_fdr05 ;
  tables  site_type * flag_methdiff_gt_20perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05 * flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc1_meth20_fdr05 ;
  tables  site_type * flag_methdiff_gt_30perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc1_meth30_fdr05 ;
  tables  site_type * flag_methdiff_gt_50perc * flag_FET_fdr05*flag_binomial_fdr05 *flag_logit_fdr05* flag_in_credible * flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc1_meth50_fdr05 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr01 * flag_in_credible 
* flag_in_fet * flag_in_binomial * flag_in_logit/ out=dmc1_meth50_fdr052 ;
run;

proc print data=dmc01_meth0_fdr05;
proc print data=dmc01_meth10_fdr05;
proc print data=dmc01_meth20_fdr05;
proc print data=dmc01_meth30_fdr05;
proc print data=dmc01_meth50_fdr05;

proc print data=dmc1_meth0_fdr05;
proc print data=dmc1_meth10_fdr05;
proc print data=dmc1_meth20_fdr05;
proc print data=dmc1_meth30_fdr05;
proc print data=dmc1_meth50_fdr052;
run;


/*






*/



data all_dmc_tests_01_2;
  set all_dmc_tests_01;
  if flag_FET_fdr05 = . then flag_FET_fdr05 = 0;
  if flag_FET_fdr01 = . then flag_FET_fdr01 = 0;
  if flag_binomial_fdr05 = . then flag_binomial_fdr05 = 0;
  if flag_binomial_fdr01 = . then flag_binomial_fdr01 = 0;
run;

data all_dmc_tests_1_2;
  set all_dmc_tests_1;
  if flag_FET_fdr05 = . then flag_FET_fdr05 = 0;
  if flag_FET_fdr01 = . then flag_FET_fdr01 = 0;
  if flag_binomial_fdr05 = . then flag_binomial_fdr05 = 0;
  if flag_binomial_fdr01 = . then flag_binomial_fdr01 = 0;
run;






proc freq data=all_dmc_tests_01_2 noprint;
  tables  site_type  * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05  / out=dmc01_meth0_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05  / out=dmc01_meth10_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_20perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05 / out=dmc01_meth20_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_30perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05 / out=dmc01_meth30_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_50perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05 / out=dmc01_meth50_fdr05_2 ;
  tables  site_type  * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc01_meth0_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc01_meth10_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_20perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc01_meth20_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_30perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc01_meth30_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_50perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01 / out=dmc01_meth50_fdr01_2 ;
run;



proc freq data=all_dmc_tests_1_2 noprint;
  tables  site_type  * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05  / out=dmc1_meth0_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05   / out=dmc1_meth10_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_20perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05   / out=dmc1_meth20_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_30perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05   / out=dmc1_meth30_fdr05_2 ;
  tables  site_type * flag_methdiff_gt_50perc * flag_FET_fdr05*flag_binomial_fdr05*flag_logit_fdr05   / out=dmc1_meth50_fdr05_2 ;
  tables  site_type  * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc1_meth0_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_10perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc1_meth10_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_20perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc1_meth20_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_30perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc1_meth30_fdr01_2 ;
  tables  site_type * flag_methdiff_gt_50perc * flag_FET_fdr01*flag_binomial_fdr01*flag_logit_fdr01  / out=dmc1_meth50_fdr01_2 ;
run;


proc print data=dmc01_meth0_fdr05_2;
proc print data=dmc01_meth10_fdr05_2;
proc print data=dmc01_meth20_fdr05_2;
proc print data=dmc01_meth30_fdr05_2;
proc print data=dmc01_meth50_fdr05_2;
proc print data=dmc01_meth0_fdr01_2;
proc print data=dmc01_meth10_fdr01_2;
proc print data=dmc01_meth20_fdr01_2;
proc print data=dmc01_meth30_fdr01_2;
proc print data=dmc01_meth50_fdr01_2;

proc print data=dmc1_meth0_fdr05_2;
proc print data=dmc1_meth10_fdr05_2;
proc print data=dmc1_meth20_fdr05_2;
proc print data=dmc1_meth30_fdr05_2;
proc print data=dmc1_meth50_fdr05_2;
proc print data=dmc1_meth0_fdr01_2;
proc print data=dmc1_meth10_fdr01_2;
proc print data=dmc1_meth20_fdr01_2;
proc print data=dmc1_meth30_fdr01_2;
proc print data=dmc1_meth50_fdr01_2;
run;


/*  DACs */


%macro prepData(trt,units);

data counts;
set wgbslocA.methylation_data_gc;
if chr = "Mt" then delete;
if chr = "Pt" then delete;
where units="&units." and treatment="&trt.";
run;

proc sort data=counts;
by chr start_pos stop_pos site_type treatment units;
proc means data=counts noprint;
by chr start_pos stop_pos site_type treatment units;
var total_C methyl_C perc_methyl;
output out=counts2 sum(total_C)=total_C_&trt._&units. sum(methyl_C)=methyl_C_&trt._&units. mean(perc_methyl)=perc_methyl_&trt._&units.;
run;

data cov;
set wgbslocA.flag_coverage_10x_gc;
where flag_&trt._&units._coverage_ge_10x=1;
keep site_type chr start_pos stop_pos;
run;

proc sort data=cov;
by site_type chr start_pos stop_pos;
proc sort data=counts2;
by site_type chr start_pos stop_pos;
run;

data counts_&trt._&units.;
merge counts2 (in=in1) cov (in=in2);
by site_type chr start_pos stop_pos;
  if in1 and in2;
  run;

%mend;



%prepData(0Gy,0U);
%prepData(01Gy,0U);
%prepData(1Gy,0U);
%prepData(0Gy,100U);
%prepData(01Gy,100U);
%prepData(1Gy,100U);


proc sort data=counts_0gy_0U;
   by site_type chr start_pos stop_pos;
proc sort data=counts_0gy_100U;
   by site_type chr start_pos stop_pos;
proc sort data=counts_01gy_0U;
   by site_type chr start_pos stop_pos;
proc sort data=counts_01gy_100U;
   by site_type chr start_pos stop_pos;
proc sort data=counts_1gy_0U;
   by site_type chr start_pos stop_pos;
proc sort data=counts_1gy_100U;
   by site_type chr start_pos stop_pos;
run;


data counts_0gy_01gy;
  merge counts_0gy_0U (in=in1) counts_0gy_100U (in=in2)
        counts_01gy_0U (in=in3) counts_01gy_100U (in=in4);
   by site_type chr start_pos stop_pos;
   if in1 and in2 and in3 and in4;
run;

data counts_0gy_1gy;
  merge counts_0gy_0U (in=in1) counts_0gy_100U (in=in2)
        counts_1gy_0U (in=in3) counts_1gy_100U (in=in4);
   by site_type chr start_pos stop_pos;
   if in1 and in2 and in3 and in4;
run;

data counts_0gy_01gy_methdiff;
  set counts_0gy_01gy;
   methdiff_0Gy = perc_methyl_0gy_100U - perc_methyl_0gy_0U ;
   methdiff_01Gy = perc_methyl_01gy_100U - perc_methyl_01gy_0U ;
   methdiff_01Gy_0Gy = methdiff_01Gy - methdiff_0Gy ;
   if abs(methdiff_01Gy_0Gy) > 0.1 then flag_methdiff_gt_10perc = 1;
   else flag_methdiff_gt_10perc = 0;
   if abs(methdiff_01Gy_0Gy) > 0.2 then flag_methdiff_gt_20perc = 1;
   else flag_methdiff_gt_20perc = 0;
   if abs(methdiff_01Gy_0Gy) > 0.3 then flag_methdiff_gt_30perc = 1;
   else flag_methdiff_gt_30perc = 0;
   if abs(methdiff_01Gy_0Gy) > 0.5 then flag_methdiff_gt_50perc = 1;
   else flag_methdiff_gt_50perc = 0;
run;

data counts_0gy_1gy_methdiff;
  set counts_0gy_1gy;
   methdiff_0Gy = perc_methyl_0gy_100U - perc_methyl_0gy_0U ;
   methdiff_1Gy = perc_methyl_1gy_100U - perc_methyl_1gy_0U ;
   methdiff_1Gy_0Gy = methdiff_1Gy - methdiff_0Gy ;
   if abs(methdiff_1Gy_0Gy) > 0.1 then flag_methdiff_gt_10perc = 1;
   else flag_methdiff_gt_10perc = 0;
   if abs(methdiff_1Gy_0Gy) > 0.2 then flag_methdiff_gt_20perc = 1;
   else flag_methdiff_gt_20perc = 0;
   if abs(methdiff_1Gy_0Gy) > 0.3 then flag_methdiff_gt_30perc = 1;
   else flag_methdiff_gt_30perc = 0;
   if abs(methdiff_1Gy_0Gy) > 0.5 then flag_methdiff_gt_50perc = 1;
   else flag_methdiff_gt_50perc = 0;
run;


/* Open/closed chromatin */

data FETbyDose_0gy;
set wgbsA.FET_DAC_0gy_100U_0gy_0U_: ;
where XP2_FISH ne .;
keep site_type chr start_pos stop_pos XP2_FISH;
rename XP2_FISH = FET_P_0gy_100U_0U;
run;

data FETbyDose_01gy;
set wgbsA.FET_DAC_01gy_100U_01gy_0U_: ;
where XP2_FISH ne .;
keep site_type chr start_pos stop_pos XP2_FISH;
rename XP2_FISH = FET_P_01gy_100U_0U;
run;

data FETbyDose_1gy;
set wgbsA.FET_DAC_1gy_100U_1gy_0U_: ;
where XP2_FISH ne .;
keep site_type chr start_pos stop_pos XP2_FISH;
rename XP2_FISH = FET_P_1gy_100U_0U;
run;

proc multtest inpvalues(FET_P_0gy_100U_0U)=FETbyDose_0gy
              fdr noprint
              out=FETbyDose_0gy_fdr(rename=(fdr_p=FET_FDR_P_0gy_100U_0U));
run;

proc multtest inpvalues(FET_P_01gy_100U_0U)=FETbyDose_01gy
              fdr noprint
              out=FETbyDose_01gy_fdr(rename=(fdr_p=FET_FDR_P_01gy_100U_0U));
run;

proc multtest inpvalues(FET_P_1gy_100U_0U)=FETbyDose_1gy
              fdr noprint
              out=FETbyDose_1gy_fdr(rename=(fdr_p=FET_FDR_P_1gy_100U_0U));
run;

/* CMH test */

data cmh1_01gy_0gy_fet_0 cmh1_01gy_0gy_fet_01 cmh1_01gy_0gy;
   set wgbsA.cmh_dac_01gy_0gy_: ;
   if condition = "01Gy" and XP2_FISH ne . then output cmh1_01gy_0gy_fet_01;
   if condition = "0Gy" and XP2_FISH ne . then output cmh1_01gy_0gy_fet_0;
   if P_CMHGA ne . then output cmh1_01gy_0gy;
run;

data cmh2_01gy_0gy_fet_0 cmh2_01gy_0gy_fet_100 cmh2_01gy_0gy;
   set wgbsA.cmh2_dac_01gy_0gy_: ;
   if unit = "0U" and XP2_FISH ne . then output cmh2_01gy_0gy_fet_0;
   if unit = "100U" and XP2_FISH ne . then output cmh2_01gy_0gy_fet_100;
   if P_CMHGA ne . then output cmh2_01gy_0gy;
run;

data cmh1_01gy_0gy_fet_0_2;
  set cmh1_01gy_0gy_fet_0;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_0gy_100U_0U _RROR_=OR_0gy_100U_0U;
run;

data cmh1_01gy_0gy_fet_01_2;
  set cmh1_01gy_0gy_fet_01;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_01gy_100U_0U _RROR_=OR_01gy_100U_0U;
run;

data cmh1_01gy_0gy_2;
  set cmh1_01gy_0gy;
  keep site_type chr start_pos stop_pos P_CMHGA P_BDCHI _MHOR_;
  rename  P_CMHGA = CMH_P_01Gy_0gy P_BDCHI=BreslowDay_P_01Gy_0gy  _MHOR_=OR_CMH_01Gy_0Gy; 
run;


data cmh2_01gy_0gy_fet_0_2;
  set cmh2_01gy_0gy_fet_0;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_0U_01Gy_0Gy _RROR_=OR_0gy_100U_0U;
run;

data cmh2_01gy_0gy_fet_100_2;
  set cmh2_01gy_0gy_fet_100;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_100U_01Gy_0Gy _RROR_=OR_01gy_100U_0U;
run;

data cmh2_01gy_0gy_2;
  set cmh2_01gy_0gy;
  keep site_type chr start_pos stop_pos P_CMHGA P_BDCHI _MHOR_;
  rename  P_CMHGA = CMH_P_01Gy_0gy_v2  P_BDCHI=BreslowDay_P_01Gy_0gy_v2  _MHOR_=OR_CMH_01Gy_0Gy_v2;
run;


proc multtest inpvalues(FET2_P_0gy_100U_0U)=cmh1_01gy_0gy_fet_0_2 fdr noprint out=cmh1_01gy_0gy_fet_0_fdr(rename=(fdr_p=FET2_FDR_P_0gy_100U_0U)); run;
proc multtest inpvalues(FET2_P_01gy_100U_0U)=cmh1_01gy_0gy_fet_01_2 fdr noprint out=cmh1_01gy_0gy_fet_01_fdr(rename=(fdr_p=FET2_FDR_P_01gy_100U_0U)); run;
proc multtest inpvalues(CMH_P_01Gy_0gy)=cmh1_01gy_0gy_2 fdr noprint out=cmh1_01gy_0gy_fdr(rename=(fdr_p=CMH_FDR_P_01Gy_0gy)); run;

proc multtest inpvalues(FET2_P_0U_01Gy_0Gy)=cmh2_01gy_0gy_fet_0_2 fdr noprint out=cmh2_01gy_0gy_fet_0_fdr(rename=(fdr_p=FET2_FDR_P_0U_01Gy_0Gy)); run;
proc multtest inpvalues(FET2_P_100U_01Gy_0Gy)=cmh2_01gy_0gy_fet_100_2 fdr noprint out=cmh2_01gy_0gy_fet_100_fdr(rename=(fdr_p=FET2_FDR_P_100U_01Gy_0Gy)); run;
proc multtest inpvalues(CMH_P_01Gy_0gy_v2)=cmh2_01gy_0gy_2 fdr noprint out=cmh2_01gy_0gy_fdr(rename=(fdr_p=CMH_FDR_P_01Gy_0gy_v2)); run;


proc multtest inpvalues(BreslowDay_P_01Gy_0gy)=cmh1_01gy_0gy_fdr fdr noprint out=cmh1_01gy_0gy_fdr2(rename=(fdr_p=BD_FDR_P_01Gy_0gy)); run;
proc multtest inpvalues(BreslowDay_P_01Gy_0gy_v2)=cmh2_01gy_0gy_fdr fdr noprint out=cmh2_01gy_0gy_fdr2(rename=(fdr_p=BD_FDR_P_01Gy_0gy_v2)); run;


proc sort data=cmh1_01gy_0gy_fet_0_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh1_01gy_0gy_fet_01_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh1_01gy_0gy_fdr2;
  by site_type chr start_pos stop_pos;
proc sort data=cmh2_01gy_0gy_fet_0_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh2_01gy_0gy_fet_100_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh2_01gy_0gy_fdr2;
  by site_type chr start_pos stop_pos;
run;

data all_cmh_01gy_0gy;
   merge cmh1_01gy_0gy_fet_0_fdr (in=in1) cmh1_01gy_0gy_fet_01_fdr (in=in2)
         cmh1_01gy_0gy_fdr2 (in=in3) cmh2_01gy_0gy_fet_0_fdr (in=in4)
         cmh2_01gy_0gy_fet_100_fdr (in=in5) cmh2_01gy_0gy_fdr2 (in=in6);
   by site_type chr start_pos stop_pos;
run;






data cmh1_1gy_0gy_fet_0 cmh1_1gy_0gy_fet_1 cmh1_1gy_0gy;
   set wgbsA.cmh_dac_1gy_0gy_: ;
   if condition = "1Gy" and XP2_FISH ne . then output cmh1_1gy_0gy_fet_1;
   if condition = "0Gy" and XP2_FISH ne . then output cmh1_1gy_0gy_fet_0;
   if P_CMHGA ne . then output cmh1_1gy_0gy;
run;

data cmh2_1gy_0gy_fet_0 cmh2_1gy_0gy_fet_100 cmh2_1gy_0gy;
   set wgbsA.cmh2_dac_1gy_0gy_: ;
   if unit = "0U" and XP2_FISH ne . then output cmh2_1gy_0gy_fet_0;
   if unit = "100U" and XP2_FISH ne . then output cmh2_1gy_0gy_fet_100;
   if P_CMHGA ne . then output cmh2_1gy_0gy;
run;


data cmh1_1gy_0gy_fet_0_2;
  set cmh1_1gy_0gy_fet_0;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_0gy_100U_0U _RROR_=OR_0gy_100U_0U;
run;

data cmh1_1gy_0gy_fet_1_2;
  set cmh1_1gy_0gy_fet_1;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_1gy_100U_0U _RROR_=OR_1gy_100U_0U;
run;

data cmh1_1gy_0gy_2;
  set cmh1_1gy_0gy;
  keep site_type chr start_pos stop_pos P_CMHGA P_BDCHI _MHOR_;
  rename  P_CMHGA = CMH_P_1Gy_0gy P_BDCHI=BreslowDay_P_1Gy_0gy  _MHOR_=OR_CMH_1Gy_0Gy; 
run;


data cmh2_1gy_0gy_fet_0_2;
  set cmh2_1gy_0gy_fet_0;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_0U_1Gy_0Gy _RROR_=OR_0gy_100U_0U;
run;

data cmh2_1gy_0gy_fet_100_2;
  set cmh2_1gy_0gy_fet_100;
  keep site_type chr start_pos stop_pos XP2_FISH _RROR_;
  rename  XP2_FISH = FET2_P_100U_1Gy_0Gy _RROR_=OR_1gy_100U_0U;
run;

data cmh2_1gy_0gy_2;
  set cmh2_1gy_0gy;
  keep site_type chr start_pos stop_pos P_CMHGA P_BDCHI _MHOR_;
  rename  P_CMHGA = CMH_P_1Gy_0gy_v2  P_BDCHI=BreslowDay_P_1Gy_0gy_v2  _MHOR_=OR_CMH_1Gy_0Gy_v2;
run;



proc multtest inpvalues(FET2_P_0gy_100U_0U)=cmh1_1gy_0gy_fet_0_2 fdr noprint out=cmh1_1gy_0gy_fet_0_fdr(rename=(fdr_p=FET2_FDR_P_0gy_100U_0U)); run;
proc multtest inpvalues(FET2_P_1gy_100U_0U)=cmh1_1gy_0gy_fet_1_2 fdr noprint out=cmh1_1gy_0gy_fet_1_fdr(rename=(fdr_p=FET2_FDR_P_1gy_100U_0U)); run;
proc multtest inpvalues(CMH_P_1Gy_0gy)=cmh1_1gy_0gy_2 fdr noprint out=cmh1_1gy_0gy_fdr(rename=(fdr_p=CMH_FDR_P_1Gy_0gy)); run;

proc multtest inpvalues(FET2_P_0U_1Gy_0Gy)=cmh2_1gy_0gy_fet_0_2 fdr noprint out=cmh2_1gy_0gy_fet_0_fdr(rename=(fdr_p=FET2_FDR_P_0U_1Gy_0Gy)); run;
proc multtest inpvalues(FET2_P_100U_1Gy_0Gy)=cmh2_1gy_0gy_fet_100_2 fdr noprint out=cmh2_1gy_0gy_fet_100_fdr(rename=(fdr_p=FET2_FDR_P_100U_1Gy_0Gy)); run;
proc multtest inpvalues(CMH_P_1Gy_0gy_v2)=cmh2_1gy_0gy_2 fdr noprint out=cmh2_1gy_0gy_fdr(rename=(fdr_p=CMH_FDR_P_1Gy_0gy_v2)); run;


proc multtest inpvalues(BreslowDay_P_1Gy_0gy)=cmh1_1gy_0gy_fdr fdr noprint out=cmh1_1gy_0gy_fdr2(rename=(fdr_p=BD_FDR_P_1Gy_0gy)); run;
proc multtest inpvalues(BreslowDay_P_1Gy_0gy_v2)=cmh2_1gy_0gy_fdr fdr noprint out=cmh2_1gy_0gy_fdr2(rename=(fdr_p=BD_FDR_P_1Gy_0gy_v2)); run;


proc sort data=cmh1_1gy_0gy_fet_0_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh1_1gy_0gy_fet_1_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh1_1gy_0gy_fdr2;
  by site_type chr start_pos stop_pos;
proc sort data=cmh2_1gy_0gy_fet_0_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh2_1gy_0gy_fet_100_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=cmh2_1gy_0gy_fdr2;
  by site_type chr start_pos stop_pos;
run;

data all_cmh_1gy_0gy;
   merge cmh1_1gy_0gy_fet_0_fdr (in=in1) cmh1_1gy_0gy_fet_1_fdr (in=in2)
         cmh1_1gy_0gy_fdr2 (in=in3) cmh2_1gy_0gy_fet_0_fdr (in=in4)
         cmh2_1gy_0gy_fet_100_fdr (in=in5) cmh2_1gy_0gy_fdr2 (in=in6);
   by site_type chr start_pos stop_pos;
run;


/* binomial DACs */



data binomial_dac_01;
  set wgbsA.binom_dac_0gy_01gy_: ;
  where effect = "treatment*units";
  if ProbChiSq = . then delete;
run;

data binomial_dac_1;
  set wgbsA.binom_dac_0gy_1gy_: ;
  where effect = "treatment*units";
  if ProbChiSq = . then delete;
run;

proc sort data=binomial_dac_01;
  by site_type;
proc multtest inpvalues(ProbChiSq)=binomial_dac_01 fdr out=binomial_dac_01_fdr(rename=(fdr_p=binomial_fdr_p ProbChiSq=binomial_p)) noprint;
  by site_type;
run;

proc sort data=binomial_dac_1;
  by site_type;
proc multtest inpvalues(ProbChiSq)=binomial_dac_1 fdr out=binomial_dac_1_fdr(rename=(fdr_p=binomial_fdr_p ProbChiSq=binomial_p)) noprint;
  by site_type;
run;



/* logit DACs */



data logit_dac_01;
  set wgbsA.logit_dac_0gy_01gy_: ;
  where effect = "treatment*units";
  if ProbChiSq = . then delete;
run;

data logit_dac_1;
  set wgbsA.logit_dac_0gy_1gy_: ;
  where effect = "treatment*units";
  if ProbChiSq = . then delete;
run;

proc sort data=logit_dac_01;
  by site_type;
proc multtest inpvalues(ProbChiSq)=logit_dac_01 fdr out=logit_dac_01_fdr(rename=(fdr_p=logit_fdr_p ProbChiSq=logit_p)) noprint;
  by site_type;
run;

proc sort data=logit_dac_1;
  by site_type;
proc multtest inpvalues(ProbChiSq)=logit_dac_1 fdr out=logit_dac_1_fdr(rename=(fdr_p=logit_fdr_p ProbChiSq=logit_p)) noprint;
  by site_type;
run;


proc sort data=FETbyDose_0gy_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=FETbyDose_01gy_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=FETbyDose_1gy_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=all_cmh_01gy_0gy;
  by site_type chr start_pos stop_pos;
proc sort data=all_cmh_1gy_0gy;
  by site_type chr start_pos stop_pos;
proc sort data=binomial_dac_01_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=binomial_dac_1_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=logit_dac_01_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=logit_dac_1_fdr;
  by site_type chr start_pos stop_pos;
proc sort data=counts_0gy_01gy_methdiff;
  by site_type chr start_pos stop_pos;
proc sort data=counts_0gy_1gy_methdiff;
  by site_type chr start_pos stop_pos;
run;



data all_DAC_tests_01Gy_0gy;
  merge FETbyDose_0gy_fdr FETbyDose_01gy_fdr all_cmh_01gy_0gy binomial_dac_01_fdr counts_0gy_01gy_methdiff logit_dac_01_fdr;
  by site_type chr start_pos stop_pos;
run;


data all_DAC_tests_1Gy_0gy;
  merge FETbyDose_0gy_fdr FETbyDose_1gy_fdr all_cmh_1gy_0gy binomial_dac_1_fdr counts_0gy_1gy_methdiff logit_dac_1_fdr; 
  by site_type chr start_pos stop_pos;
run;



data dac_01_check;
  set all_DAC_tests_01Gy_0gY;
  where FET_FDR_P_0gy_100U_0U ne . and  FET_FDR_P_01gy_100U_0U ne . and
        FET2_FDR_P_0gy_100U_0U ne . and  FET2_FDR_P_01gy_100U_0U ne . and
        CMH_FDR_P_01gy_0gy ne . and BD_FDR_P_01Gy_0gy ne . and FET2_FDR_P_0U_01gy_0gy ne .
        and FET2_FDR_P_100U_01gy_0gy ne . and CMH_FDR_P_01gy_0gy_v2 ne . and BD_FDR_P_01Gy_0gy_v2 ne .;
run;

*3535414 obs;

proc corr data=dac_01_check pearson out=dac_01_check_pearson noprint;
  by site_type;
  var FET_FDR_P_0gy_100U_0U FET_FDR_P_01gy_100U_0U 
        FET2_FDR_P_0gy_100U_0U FET2_FDR_P_01gy_100U_0U
        CMH_FDR_P_01gy_0gy BD_FDR_P_01Gy_0gy FET2_FDR_P_0U_01gy_0gy 
        FET2_FDR_P_100U_01gy_0gy  CMH_FDR_P_01gy_0gy_v2  BD_FDR_P_01Gy_0gy_v2 binomial_fdr_p logit_fdr_p;
run;

proc print data=dac_01_check_pearson;
run;

/*

                                                                          FET_FDR_                    FET2_FDR_
          site_                                           FET_FDR_P_       P_01gy_    FET2_FDR_P_       P_01gy_     CMH_FDR_P_
   Obs    type     _TYPE_    _NAME_                      0gy_100U_0U       100U_0U    0gy_100U_0U       100U_0U       01Gy_0gy

     1     GC       MEAN                                        0.30          0.09           0.30          0.09           0.07
     2     GC       STD                                         0.33          0.22           0.33          0.22           0.18
     3     GC       N                                     3535414.00    3535414.00     3535414.00    3535414.00     3535414.00
     4     GC       CORR     FET_FDR_P_0gy_100U_0U              1.00          0.37           1.00          0.37           0.47
     5     GC       CORR     FET_FDR_P_01gy_100U_0U             0.37          1.00           0.37          1.00           0.78
     6     GC       CORR     FET2_FDR_P_0gy_100U_0U             1.00          0.37           1.00          0.37           0.47
     7     GC       CORR     FET2_FDR_P_01gy_100U_0U            0.37          1.00           0.37          1.00           0.78
     8     GC       CORR     CMH_FDR_P_01Gy_0gy                 0.47          0.78           0.47          0.78           1.00
     9     GC       CORR     BD_FDR_P_01Gy_0gy                 -0.24          0.15          -0.24          0.15           0.03
    10     GC       CORR     FET2_FDR_P_0U_01Gy_0Gy             0.00          0.00           0.00          0.00          -0.00
    11     GC       CORR     FET2_FDR_P_100U_01Gy_0Gy          -0.00          0.32          -0.00          0.32           0.22
    12     GC       CORR     CMH_FDR_P_01Gy_0gy_v2              0.10          0.23           0.10          0.23           0.20
    13     GC       CORR     BD_FDR_P_01Gy_0gy_v2              -0.25          0.15          -0.25          0.15           0.03
    14     GC       CORR     binomial_fdr_p                    -0.11          0.11          -0.11          0.11           0.02
    15     GC       CORR     logit_fdr_p                       -0.11          0.11          -0.11          0.11           0.02

                                          FET2_FDR_    CMH_FDR_P_
            BD_FDR_P_    FET2_FDR_P_        P_100U_     01Gy_0gy_     BD_FDR_P_       binomial_
   Obs       01Gy_0gy    0U_01Gy_0Gy       01Gy_0Gy        v2        01Gy_0gy_v2          fdr_p    logit_fdr_p

     1           0.76           1.00           0.51          0.45           0.75           0.97           0.97
     2           0.12           0.01           0.32          0.29           0.14           0.05           0.05
     3     3535414.00     3535414.00     3535414.00    3535414.00     3535414.00     2073271.00     2073271.00
     4          -0.24           0.00          -0.00          0.10          -0.25          -0.11          -0.11
     5           0.15           0.00           0.32          0.23           0.15           0.11           0.11
     6          -0.24           0.00          -0.00          0.10          -0.25          -0.11          -0.11
     7           0.15           0.00           0.32          0.23           0.15           0.11           0.11
     8           0.03          -0.00           0.22          0.20           0.03           0.02           0.02
     9           1.00           0.06           0.27          0.01           0.99           0.69           0.69
    10           0.06           1.00           0.01          0.01           0.05           0.20           0.20
    11           0.27           0.01           1.00          0.82           0.29           0.25           0.25
    12           0.01           0.01           0.82          1.00           0.02           0.07           0.07
    13           0.99           0.05           0.29          0.02           1.00           0.67           0.67
    14           0.69           0.20           0.25          0.07           0.67           1.00           1.00
    15           0.69           0.20           0.25          0.07           0.67           1.00           1.00


*/
ods graphics / maxobs=3535414 antialiasmax=3535414;

proc sgplot data=dac_01_check;
   scatter x=FET_FDR_P_0gy_100U_0U y=FET2_FDR_P_0gy_100U_0U;
  label FET_FDR_P_0gy_100U_0U="FET 100U vs 0U, 0Gy (standalone)"
        FET2_FDR_P_0gy_100U_0U="FET 100U vs 0U, 0Gy (CMH)"
        ;
run;

proc sgplot data=dac_01_check;
   scatter x=FET_FDR_P_01gy_100U_0U y=FET2_FDR_P_01gy_100U_0U;
  label FET_FDR_P_01gy_100U_0U="FET 100U vs 0U, 10cGy (standalone)"
        FET2_FDR_P_01gy_100U_0U="FET 100U vs 0U, 10cGy (CMH)"
        ;
run;


proc sgplot data=dac_01_check;
   scatter x=FET_FDR_P_0gy_100U_0U y=FET2_FDR_P_0U_01Gy_0Gy;
  label FET_FDR_P_0gy_100U_0U="FET 100U vs 0U, 10cGy (standalone)"
        FET2_FDR_P_0U_01Gy_0Gy="FET 0U, 10cGy vs 0cGy (CMH)"
        ;
run;

proc sgplot data=dac_01_check;
   scatter x=FET_FDR_P_01gy_100U_0U y=FET2_FDR_P_0U_01Gy_0Gy;
  label FET_FDR_P_0gy_100U_0U="FET 100U vs 0U, 10cGy (standalone)"
        FET2_FDR_P_0U_01Gy_0Gy="FET 100U, 10cGy vs 0cGy (CMH)"
        ;
run;


proc sgplot data=dac_01_check;
   scatter x=CMH_FDR_P_01gy_0gy y=CMH_FDR_P_01gy_0gy_v2;
  label CMH_FDR_P_01gy_0gy="CMH FDR P test (stratified on IR dose)"
        CMH_FDR_P_01gy_0gy_v2="CMH FDR P test (stratified on M.CviPI)"
        ;
run;


data dac_01_check;
  set all_DAC_tests_01Gy_0gY;
  where FET_FDR_P_0gy_100U_0U ne . and  FET_FDR_P_01gy_100U_0U ne . and
        FET2_FDR_P_0gy_100U_0U ne . and  FET2_FDR_P_01gy_100U_0U ne . and
        CMH_FDR_P_01gy_0gy ne . and BD_FDR_P_01Gy_0gy ne . and FET2_FDR_P_0U_01gy_0gy ne .
        and FET2_FDR_P_100U_01gy_0gy ne . and CMH_FDR_P_01gy_0gy_v2 ne . and BD_FDR_P_01Gy_0gy_v2 ne .;
run;

*3535414 obs;

proc corr data=dac_01_check pearson out=dac_01_check_pearson noprint;
  by site_type;
  var FET_FDR_P_0gy_100U_0U FET_FDR_P_01gy_100U_0U 
        FET2_FDR_P_0gy_100U_0U FET2_FDR_P_01gy_100U_0U
        CMH_FDR_P_01gy_0gy BD_FDR_P_01Gy_0gy FET2_FDR_P_0U_01gy_0gy 
        FET2_FDR_P_100U_01gy_0gy  CMH_FDR_P_01gy_0gy_v2  BD_FDR_P_01Gy_0gy_v2;
run;

proc print data=dac_01_check_pearson;
run;




/* Import results */

%let path=/TB14/TB14/sandbox/dtra_sandbox;

proc import datafile="&path./at_rad_dss_0cGy_01cGy_CG.results.txt" out=dss_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_01cGy_CHG.results.txt" out=dss_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_01cGy_CHH.results.txt" out=dss_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylkit_0gy_01gy_CG.results.overdispersion.txt" out=mk_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylkit_0cGy_01cGy_CHG.results.overdispersion.txt" out=mk_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylkit_0cGy_01cGy_CHH.results.overdispersion.txt" out=mk_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_01cGy_CG.binomial.txt" out=ms_bin_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_01cGy_CG.betabinomial.txt" out=ms_bb_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_01cGy_CHG.binomial.txt" out=ms_bin_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_01cGy_CHG.betabinomial.txt" out=ms_bb_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_01cGy_CHH.binomial.txt" out=ms_bin_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_01cGy_CHH.betabinomial.txt" out=ms_bb_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_01gy_CG.output.txt" out=met_cg dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_01gy_CHG.output.txt" out=met_chg dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_01gy_CHH.output.txt" out=met_chh dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./at_rad_methylkit_0cGy_01cGy_GC.results2.overdispersion.txt" out=mk_gc dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_01cGy_GC.treatment.txt" out=dss_gc_trt dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_01cGy_GC.trt_by_units.txt" out=dss_gc_int dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_01gy_GC.output.txt" out=met_gc dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_01gy_GC_no_negative.output.txt" out=met_gc_noneg dbms=tab replace; guessingrows=10000; getnames=no; run;








/* Stack results */

data methylkit_dmc_all;
   set mk_cg (in=in1) mk_chg (in=in2) mk_chh (in=in3);
   length site_type $3.;
   length chr2 $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   chr2 = compress(chr);
   if qvalue < 0.05 and qvalue > 0 then flag_methylkit_q05=1; else flag_methylkit_q05=0;
   keep site_type chr2 end pvalue qvalue meth_diff flag_methylkit_q05 ;
   rename chr2=chr end=stop_pos pvalue=methylkit_p qvalue=methylkit_q meth_diff=methylkit_meth_diff;
run;


data methylsig_dmc_all_bin;
   set ms_bin_cg (in=in1) ms_bin_chg (in=in2) ms_bin_chh (in=in3);
   length site_type $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   if fdr="NA" then delete;
   fdr2=fdr*1;
   if fdr2 < 0.05 and fdr2 > 0 then flag_methylsig_bin_fdr05=1;
   else if index(fdr, "e") > 0 then flag_methylsig_bin_fdr05=1; 
   else flag_methylsig_bin_fdr05=0;
   keep site_type seqnames end meth_diff pvalue  fdr2 flag_methylsig_bin_fdr05 ;
   rename seqnames=chr end=stop_pos meth_diff=methylsig_bin_meth_diff pvalue=methylsig_bin_p fdr2=methylsig_bin_fdr_p;run;

data check;
  set methylsig_dmc_all_bin;
  where flag_methylsig_bin_fdr05=1;
run;


data methylsig_dmc_all_bb;
   set ms_bb_cg (in=in1) ms_bb_chg (in=in2) ms_bb_chh (in=in3);
   length site_type $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   if fdr="NA" then delete;
   fdr2=fdr*1;
   if fdr2 < 0.05 and fdr2 > 0 then flag_methylsig_bb_fdr05=1; 
   else if index(fdr, "e") > 0 then flag_methylsig_bb_fdr05=1; 
   else flag_methylsig_bb_fdr05=0;
   keep site_type seqnames end meth_diff pvalue  fdr2 flag_methylsig_bb_fdr05 ;
   rename seqnames=chr end=stop_pos meth_diff=methylsig_bb_meth_diff pvalue=methylsig_bb_p fdr2=methylsig_bb_fdr_p;
run;


data metilene_dmc_all;
   set met_cg (in=in1) met_chg (in=in2) met_chh (in=in3);
   length site_type $3.;
   length chr2 $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   chr2 = compress(VAR1);
   if VAR4 < 0.05 then flag_metilene_q05=1; else flag_metilene_q05=0;
   keep site_type chr2 VAR3 VAR4 VAR5 VAR7 flag_metilene_q05;
   rename chr2=chr VAR3=stop_pos VAR4=metilene_q VAR5=metilene_meth_diff VAR7=metilene_MWU_P;
run;


data dss_dmc_all;
   set dss_cg (in=in1) dss_chg (in=in2) dss_chh (in=in3);
   length site_type $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH"; 
   if fdr < 0.05 then flag_dss_fdr05=1; else flag_dss_fdr05=0;
   keep site_type chr pos diff pval fdr flag_dss_fdr05;
   rename pos=stop_pos diff=dss_meth_diff pval=dss_P fdr=dss_FDR;
run;




proc sort data=methylkit_dmc_all;
  by site_type chr stop_pos;
proc sort data=methylsig_dmc_all_bin;
  by site_type chr stop_pos;
proc sort data=methylsig_dmc_all_bb;
  by site_type chr stop_pos;
proc sort data=metilene_dmc_all;
  by site_type chr stop_pos;
proc sort data=dss_dmc_all;
  by site_type chr stop_pos;
proc sort data=all_dmc_tests_01_2;
  by site_type chr stop_pos;
run;

data all_all_dmc_tests_01gy;
  merge all_dmc_tests_01_2 methylkit_dmc_all methylsig_dmc_all_bin methylsig_dmc_all_bb
   metilene_dmc_all dss_dmc_all;
   by site_type chr stop_pos;
run;


data methylkit_dac_all;
   set mk_gc ;
  length chr2 $3.;
   chr2=compress(chr);
   if qvalue < 0.05 and qvalue > 0 then flag_methylkit_q05=1; else flag_methylkit_q05=0;
   keep chr2 end pvalue qvalue meth_diff flag_methylkit_q05 ;
   rename  chr2=chr end=stop_pos pvalue=methylkit_p qvalue=methylkit_q meth_diff=methylkit_meth_diff;
run;

data metilene_dac_all;
   set met_gc;
  length chr2 $3.;
   chr2 = compress(VAR1);
   if VAR4 < 0.05 then flag_metilene_q05=1; else flag_metilene_q05=0;
   keep chr2 VAR3 VAR4 VAR5 VAR7 flag_metilene_q05;
   rename chr2=chr VAR3=stop_pos VAR4=metilene_q VAR5=metilene_meth_diff VAR7=metilene_MWU_P;
run;


data metilene_dac_all2;
   set met_gc_noneg;
  length chr2 $3.;
   chr2 = compress(VAR1);
   if VAR4 < 0.05 then flag_metilene_noneg_q05=1; else flag_metilene_noneg_q05=0;
   keep chr2 VAR3 VAR4 VAR5 VAR7 flag_metilene_noneg_q05;
   rename chr2=chr VAR3=stop_pos VAR4=metilene_noneg_q VAR5=metilene_noneg_meth_diff VAR7=metilene_noneg_MWU_P;
run;

data dss_dac_all;
   set dss_gc_trt;
   if pvals = "NA" then delete;
   pval = pvals * 1;
   fdr2 = fdrs * 1;
   if fdr2 < 0.05 then flag_dss_trt_fdr05=1; else flag_dss_trt_fdr05=0;
   keep chr pos  pval fdrs flag_dss_trt_fdr05;
   rename pos=stop_pos pval=dss_trt_P fdrs=dss_trt_FDR_P;
run;

data dss_dac_all2;
   set dss_gc_int;
   if pvals = "NA" then delete;
   pval = pvals * 1;
   fdr2 = fdrs * 1;
   if fdr2 < 0.05 then flag_dss_int_fdr05=1; else flag_dss_int_fdr05=0;
   keep chr pos pval fdrs flag_dss_int_fdr05;
   rename pos=stop_pos pval=dss_int_P fdrs=dss_int_FDR_P;
run;




proc sort data=methylkit_dac_all;
  by  chr stop_pos;
proc sort data=metilene_dac_all;
  by  chr stop_pos;
proc sort data=metilene_dac_all2;
  by  chr stop_pos;
proc sort data=dss_dac_all;
  by  chr stop_pos;
proc sort data=dss_dac_all2;
  by  chr stop_pos;
proc sort data=all_DAC_tests_01Gy_0gy;
  by  chr stop_pos;
run;

data all_all_dac_tests_01gy;
  merge all_DAC_tests_01Gy_0gy methylkit_dac_all metilene_dac_all metilene_dac_all2
   dss_dac_all dss_dac_all2;
   by  chr stop_pos;
   if (FET_FDR_P_01Gy_100U_0U ne . and FET_FDR_P_01Gy_100U_0U  < 0.05) or
    (FET_FDR_P_0Gy_100U_0U ne . and FET_FDR_P_0Gy_100U_0U  < 0.05) 
   then flag_FET_either_fdr05=1;
   else  flag_FET_either_fdr05=0;

   if (FET_FDR_P_01Gy_100U_0U ne . and FET_FDR_P_01Gy_100U_0U  < 0.01) or
    (FET_FDR_P_0Gy_100U_0U ne . and FET_FDR_P_0Gy_100U_0U  < 0.01) 
   then flag_FET_either_fdr01=1;
   else  flag_FET_either_fdr01=0;

   if CMH_FDR_P_01Gy_0gy ne . and CMH_FDR_P_01Gy_0gy < 0.05 then flag_CMH_FDR_P05=1;
   else  flag_CMH_FDR_P05=0;
   if CMH_FDR_P_01Gy_0gy ne . and CMH_FDR_P_01Gy_0gy < 0.01 then flag_CMH_FDR_P01=1;
   else  flag_CMH_FDR_P01=0;
run;



data wgbsA.all_dmc_tests_01gy_compare;
   set all_all_dmc_tests_01gy;
run;

data wgbsA.all_dac_tests_01gy_compare;
   set all_all_dac_tests_01gy;
run;












proc import datafile="&path./at_rad_dss_0cGy_1cGy_CG.results.txt" out=dss_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_1cGy_CHG.results.txt" out=dss_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_1cGy_CHH.results.txt" out=dss_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylkit_0gy_1gy_CG.results.overdispersion.txt" out=mk_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylkit_0cGy_1cGy_CHG.results.overdispersion.txt" out=mk_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylkit_0cGy_1cGy_CHH.results.overdispersion.txt" out=mk_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_1cGy_CG.binomial.txt" out=ms_bin_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_1cGy_CG.betabinomial.txt" out=ms_bb_cg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_1cGy_CHG.binomial.txt" out=ms_bin_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_1cGy_CHG.betabinomial.txt" out=ms_bb_chg dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_1cGy_CHH.binomial.txt" out=ms_bin_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_methylsig_0cGy_1cGy_CHH.betabinomial.txt" out=ms_bb_chh dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_1gy_CG.output.txt" out=met_cg dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_1gy_CHG.output.txt" out=met_chg dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_1gy_CHH.output.txt" out=met_chh dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./at_rad_methylkit_0cGy_1cGy_GC.results2.overdispersion.txt" out=mk_gc dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_1cGy_GC.treatment.txt" out=dss_gc_trt dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./at_rad_dss_0cGy_1cGy_GC.trt_by_units.txt" out=dss_gc_int dbms=tab replace; guessingrows=10000; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_1gy_GC.output.txt" out=met_gc dbms=tab replace; guessingrows=10000; getnames=no; run;
proc import datafile="&path./metilene_output/at_rad_metilene_0gy_1gy_GC_no_negative.output.txt" out=met_gc_noneg dbms=tab replace; guessingrows=10000; getnames=no; run;







/* Stack results */

data methylkit_dmc_all;
   set mk_cg (in=in1) mk_chg (in=in2) mk_chh (in=in3);
   length site_type $3.;
   length chr2 $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   chr2 = compress(chr);
   if qvalue < 0.05 and qvalue > 0 then flag_methylkit_q05=1; else flag_methylkit_q05=0;
   keep site_type chr2 end pvalue qvalue meth_diff flag_methylkit_q05 ;
   rename chr2=chr end=stop_pos pvalue=methylkit_p qvalue=methylkit_q meth_diff=methylkit_meth_diff;
run;


data methylsig_dmc_all_bin;
   set ms_bin_cg (in=in1) ms_bin_chg (in=in2) ms_bin_chh (in=in3);
   length site_type $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   if fdr="NA" then delete;
   fdr2=fdr*1;
   if fdr2 < 0.05 and fdr2 > 0 then flag_methylsig_bin_fdr05=1;
   else if index(fdr, "e") > 0 then flag_methylsig_bin_fdr05=1; 
   else flag_methylsig_bin_fdr05=0;
   keep site_type seqnames end meth_diff pvalue  fdr2 flag_methylsig_bin_fdr05 ;
   rename seqnames=chr end=stop_pos meth_diff=methylsig_bin_meth_diff pvalue=methylsig_bin_p fdr2=methylsig_bin_fdr_p;run;

data check;
  set methylsig_dmc_all_bin;
  where flag_methylsig_bin_fdr05=1;
run;


data methylsig_dmc_all_bb;
   set ms_bb_cg (in=in1) ms_bb_chg (in=in2) ms_bb_chh (in=in3);
   length site_type $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   if fdr="NA" then delete;
   fdr2=fdr*1;
   if fdr2 < 0.05 and fdr2 > 0 then flag_methylsig_bb_fdr05=1; 
   else if index(fdr, "e") > 0 then flag_methylsig_bb_fdr05=1; 
   else flag_methylsig_bb_fdr05=0;
   keep site_type seqnames end meth_diff pvalue  fdr2 flag_methylsig_bb_fdr05 ;
   rename seqnames=chr end=stop_pos meth_diff=methylsig_bb_meth_diff pvalue=methylsig_bb_p fdr2=methylsig_bb_fdr_p;
run;


data metilene_dmc_all;
   set met_cg (in=in1) met_chg (in=in2) met_chh (in=in3);
   length site_type $3.;
   length chr2 $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH";
   chr2 = compress(VAR1);
   if VAR4 < 0.05 then flag_metilene_q05=1; else flag_metilene_q05=0;
   keep site_type chr2 VAR3 VAR4 VAR5 VAR7 flag_metilene_q05;
   rename chr2=chr VAR3=stop_pos VAR4=metilene_q VAR5=metilene_meth_diff VAR7=metilene_MWU_P;
run;


data dss_dmc_all;
   set dss_cg (in=in1) dss_chg (in=in2) dss_chh (in=in3);
   length site_type $3.;
   if in1 then site_type="CG";
   if in2 then site_type="CHG";
   if in3 then site_type="CHH"; 
   if fdr < 0.05 then flag_dss_fdr05=1; else flag_dss_fdr05=0;
   keep site_type chr pos diff pval fdr flag_dss_fdr05;
   rename pos=stop_pos diff=dss_meth_diff pval=dss_P fdr=dss_FDR;
run;




proc sort data=methylkit_dmc_all;
  by site_type chr stop_pos;
proc sort data=methylsig_dmc_all_bin;
  by site_type chr stop_pos;
proc sort data=methylsig_dmc_all_bb;
  by site_type chr stop_pos;
proc sort data=metilene_dmc_all;
  by site_type chr stop_pos;
proc sort data=dss_dmc_all;
  by site_type chr stop_pos;
proc sort data=all_dmc_tests_1_2;
  by site_type chr stop_pos;
run;

data all_all_dmc_tests_1gy;
  merge all_dmc_tests_1_2 methylkit_dmc_all methylsig_dmc_all_bin methylsig_dmc_all_bb
   metilene_dmc_all dss_dmc_all;
   by site_type chr stop_pos;
run;




data wgbsA.all_dmc_tests_1gy_compare;
   set all_all_dmc_tests_1gy;
run;


data methylkit_dac_all;
   set mk_gc ;
  length chr2 $3.;
   chr2=compress(chr);
   if qvalue < 0.05 and qvalue > 0 then flag_methylkit_q05=1; else flag_methylkit_q05=0;
   keep chr2 end pvalue qvalue meth_diff flag_methylkit_q05 ;
   rename  chr2=chr end=stop_pos pvalue=methylkit_p qvalue=methylkit_q meth_diff=methylkit_meth_diff;
run;

data metilene_dac_all;
   set met_gc;
  length chr2 $3.;
   chr2 = compress(VAR1);
   if VAR4 < 0.05 then flag_metilene_q05=1; else flag_metilene_q05=0;
   keep chr2 VAR3 VAR4 VAR5 VAR7 flag_metilene_q05;
   rename chr2=chr VAR3=stop_pos VAR4=metilene_q VAR5=metilene_meth_diff VAR7=metilene_MWU_P;
run;


data metilene_dac_all2;
   set met_gc_noneg;
  length chr2 $3.;
   chr2 = compress(VAR1);
   if VAR4 < 0.05 then flag_metilene_noneg_q05=1; else flag_metilene_noneg_q05=0;
   keep chr2 VAR3 VAR4 VAR5 VAR7 flag_metilene_noneg_q05;
   rename chr2=chr VAR3=stop_pos VAR4=metilene_noneg_q VAR5=metilene_noneg_meth_diff VAR7=metilene_noneg_MWU_P;
run;

data dss_dac_all;
   set dss_gc_trt;
   if pvals = "NA" then delete;
   pval = pvals * 1;
   fdr2 = fdrs * 1;
   if fdr2 < 0.05 then flag_dss_trt_fdr05=1; else flag_dss_trt_fdr05=0;
   keep chr pos  pval fdrs flag_dss_trt_fdr05;
   rename pos=stop_pos pval=dss_trt_P fdrs=dss_trt_FDR_P;
run;

data dss_dac_all2;
   set dss_gc_int;
   if pvals = "NA" then delete;
   pval = pvals * 1;
   fdr2 = fdrs * 1;
   if fdr2 < 0.05 then flag_dss_int_fdr05=1; else flag_dss_int_fdr05=0;
   keep chr pos pval fdrs flag_dss_int_fdr05;
   rename pos=stop_pos pval=dss_int_P fdrs=dss_int_FDR_P;
run;




proc sort data=methylkit_dac_all;
  by  chr stop_pos;
proc sort data=metilene_dac_all;
  by  chr stop_pos;
proc sort data=metilene_dac_all2;
  by  chr stop_pos;
proc sort data=dss_dac_all;
  by  chr stop_pos;
proc sort data=dss_dac_all2;
  by  chr stop_pos;
proc sort data=all_DAC_tests_1Gy_0gy;
  by  chr stop_pos;
run;

data all_all_dac_tests_1gy;
  merge all_DAC_tests_1Gy_0gy methylkit_dac_all metilene_dac_all metilene_dac_all2
   dss_dac_all dss_dac_all2;
   by  chr stop_pos;
   if (FET_FDR_P_1Gy_100U_0U ne . and FET_FDR_P_1Gy_100U_0U  < 0.05) or
    (FET_FDR_P_0Gy_100U_0U ne . and FET_FDR_P_0Gy_100U_0U  < 0.05) 
   then flag_FET_either_fdr05=1;
   else  flag_FET_either_fdr05=0;

   if (FET_FDR_P_1Gy_100U_0U ne . and FET_FDR_P_1Gy_100U_0U  < 0.01) or
    (FET_FDR_P_0Gy_100U_0U ne . and FET_FDR_P_0Gy_100U_0U  < 0.01) 
   then flag_FET_either_fdr01=1;
   else  flag_FET_either_fdr01=0;

   if CMH_FDR_P_1Gy_0gy ne . and CMH_FDR_P_1Gy_0gy < 0.05 then flag_CMH_FDR_P05=1;
   else  flag_CMH_FDR_P05=0;
   if CMH_FDR_P_1Gy_0gy ne . and CMH_FDR_P_1Gy_0gy < 0.01 then flag_CMH_FDR_P01=1;
   else  flag_CMH_FDR_P01=0;
run;



data wgbsA.all_dac_tests_1gy_compare;
   set all_all_dac_tests_1gy;
run;




