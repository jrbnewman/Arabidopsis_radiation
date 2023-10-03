

%macro methDataGen(siteType);


/* DMR counts:
   Hyper all, exon, intron, promoter, downstream; 0.1G, 1G
   Hypo all, exon, intron, promoter, downstream; 0.1G, 1G

 */
%let siteType=CHH;

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

