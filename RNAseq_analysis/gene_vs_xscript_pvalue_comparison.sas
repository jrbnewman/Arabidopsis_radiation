/* check that DE genes from transcript-level data are the same as the previous set of genes */

ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* DE genes (from exons) */

data de_genes;
  set rs.fdr_by_gene_v2;
  if p_treatment < 0.05 and p_treatment > 0 then flag_gene_trt_p05=1;
  else flag_gene_trt_p05=0;

  if p_time < 0.05 and p_time > 0 then flag_gene_time_p05=1;
  else flag_gene_time_p05=0;

  if p_trt_by_time < 0.05 and p_trt_by_time > 0 then flag_gene_trtime_p05=1;
  else flag_gene_trtime_p05=0;

  if p_01gy_v_Mock_1h < 0.05 and p_01gy_v_Mock_1h > 0 then flag_gene_01vM_1h_p05=1;
  else flag_gene_01vM_1h_p05=0;
  if p_01gy_v_Mock_3h < 0.05 and p_01gy_v_Mock_3h > 0 then flag_gene_01vM_3h_p05=1;
  else flag_gene_01vM_3h_p05=0;
  if p_01gy_v_Mock_24h < 0.05 and p_01gy_v_Mock_24h > 0 then flag_gene_01vM_24h_p05=1;
  else flag_gene_01vM_24h_p05=0;
  if p_01gy_v_Mock_72h < 0.05 and p_01gy_v_Mock_72h > 0 then flag_gene_01vM_72h_p05=1;
  else flag_gene_01vM_72h_p05=0;

  if p_01gy_v_Mock_1h < 0.05 and p_01gy_v_Mock_1h > 0 then flag_gene_1vM_1h_p05=1;
  else flag_gene_1vM_1h_p05=0;
  if p_01gy_v_Mock_3h < 0.05 and p_01gy_v_Mock_3h > 0 then flag_gene_1vM_3h_p05=1;
  else flag_gene_1vM_3h_p05=0;
  if p_01gy_v_Mock_24h < 0.05 and p_01gy_v_Mock_24h > 0 then flag_gene_1vM_24h_p05=1;
  else flag_gene_1vM_24h_p05=0;
  if p_01gy_v_Mock_72h < 0.05 and p_01gy_v_Mock_72h > 0 then flag_gene_1vM_72h_p05=1;
  else flag_gene_1vM_72h_p05=0;

  keep gene_id flag_gene_: ;
run;


/* DE transcripts */

data de_xscript;
  set rs.fdr_by_transcript;
  if p_treatment < 0.05 and p_treatment > 0 then flag_xs_trt_p05=1;
  else flag_xs_trt_p05=0;

  if p_time < 0.05 and p_time > 0 then flag_xs_time_p05=1;
  else flag_xs_time_p05=0;

  if p_trt_by_time < 0.05 and p_trt_by_time > 0 then flag_xs_trtime_p05=1;
  else flag_xs_trtime_p05=0;

  if p_01gy_v_Mock_1h < 0.05 and p_01gy_v_Mock_1h > 0 then flag_xs_01vM_1h_p05=1;
  else flag_xs_01vM_1h_p05=0;
  if p_01gy_v_Mock_3h < 0.05 and p_01gy_v_Mock_3h > 0 then flag_xs_01vM_3h_p05=1;
  else flag_xs_01vM_3h_p05=0;
  if p_01gy_v_Mock_24h < 0.05 and p_01gy_v_Mock_24h > 0 then flag_xs_01vM_24h_p05=1;
  else flag_xs_01vM_24h_p05=0;
  if p_01gy_v_Mock_72h < 0.05 and p_01gy_v_Mock_72h > 0 then flag_xs_01vM_72h_p05=1;
  else flag_xs_01vM_72h_p05=0;

  if p_01gy_v_Mock_1h < 0.05 and p_01gy_v_Mock_1h > 0 then flag_xs_1vM_1h_p05=1;
  else flag_xs_1vM_1h_p05=0;
  if p_01gy_v_Mock_3h < 0.05 and p_01gy_v_Mock_3h > 0 then flag_xs_1vM_3h_p05=1;
  else flag_xs_1vM_3h_p05=0;
  if p_01gy_v_Mock_24h < 0.05 and p_01gy_v_Mock_24h > 0 then flag_xs_1vM_24h_p05=1;
  else flag_xs_1vM_24h_p05=0;
  if p_01gy_v_Mock_72h < 0.05 and p_01gy_v_Mock_72h > 0 then flag_xs_1vM_72h_p05=1;
  else flag_xs_1vM_72h_p05=0;

  keep transcript_id flag_xs_: ;
run;

data xs2gene;
  set tair.tair20_exons;
  length transcript_id2 $50.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output; end;
  keep transcript_id2 gene_id;
  rename transcript_id2=transcript_id;
run;

proc sort data=xs2gene nodup;
  by transcript_id gene_id;
proc sort data=de_xscript;
  by transcript_id;
run;

data de_xscript2;
  merge de_xscript (in=in1) xs2gene (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=de_xscript2;
  by gene_id;
proc means data=de_xscript2 noprint;
  by gene_id;
  var flag_xs_: ;
  output out=de_xs_by_gene max=;
run;

proc sort data=de_xs_by_gene;
  by gene_id;
proc sort data=de_genes;
  by gene_id;
run;

data de_gene_v_xs;
  merge de_genes (in=in1) de_xs_by_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;
/*
 There were 21738 observations read from the data set WORK.DE_GENES.
 There were 22790 observations read from the data set WORK.DE_XS_BY_GENE.
 The data set WORK.DE_GENE_V_XS has 20651 observations and 25 variables.
*/

proc freq data=de_gene_v_xs;
   tables flag_gene_trt_p05*flag_xs_trt_p05
          flag_gene_time_p05*flag_xs_time_p05
          flag_gene_trtime_p05*flag_xs_trtime_p05
          flag_gene_01vM_1h_p05*flag_xs_01vM_1h_p05
          flag_gene_01vM_3h_p05*flag_xs_01vM_3h_p05
          flag_gene_01vM_24h_p05*flag_xs_01vM_24h_p05
          flag_gene_01vM_72h_p05*flag_xs_01vM_72h_p05
          flag_gene_1vM_1h_p05*flag_xs_1vM_1h_p05
          flag_gene_1vM_3h_p05*flag_xs_1vM_3h_p05
          flag_gene_1vM_24h_p05*flag_xs_1vM_24h_p05
          flag_gene_1vM_72h_p05*flag_xs_1vM_72h_p05;
run;

/*
 flag_gene_trt_p05
           flag_xs_trt_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  11923 |   1780 |  13703
          |  57.74 |   8.62 |  66.36
          |  87.01 |  12.99 |
          |  88.32 |  24.89 |
 ---------+--------+--------+
        1 |   1577 |   5371 |   6948
          |   7.64 |  26.01 |  33.64
          |  22.70 |  77.30 |
          |  11.68 |  75.11 |
 ---------+--------+--------+
 Total       13500     7151    20651
             65.37    34.63   100.00

            The SAS System            11:33 Monda

 flag_gene_time_p05
           flag_xs_time_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   2931 |    750 |   3681
          |  14.19 |   3.63 |  17.82
          |  79.63 |  20.37 |
          |  71.56 |   4.53 |
 ---------+--------+--------+
        1 |   1165 |  15805 |  16970
          |   5.64 |  76.53 |  82.18
          |   6.87 |  93.13 |
          |  28.44 |  95.47 |
 ---------+--------+--------+
 Total        4096    16555    20651
             19.83    80.17   100.00

            The SAS System            1

 flag_gene_trtime_p05
           flag_xs_trtime_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  12752 |   1817 |  14569
          |  61.75 |   8.80 |  70.55
          |  87.53 |  12.47 |
          |  87.95 |  29.54 |
 ---------+--------+--------+
        1 |   1747 |   4335 |   6082
          |   8.46 |  20.99 |  29.45
          |  28.72 |  71.28 |
          |  12.05 |  70.46 |
 ---------+--------+--------+
 Total       14499     6152    20651
             70.21    29.79   100.00

            The SAS System            11:33

flag_gene_01vM_1h_p05
          flag_xs_01vM_1h_p05

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  12700 |   1420 |  14120
         |  61.50 |   6.88 |  68.37
         |  89.94 |  10.06 |
         |  89.84 |  21.80 |
---------+--------+--------+
       1 |   1437 |   5094 |   6531
         |   6.96 |  24.67 |  31.63
         |  22.00 |  78.00 |
         |  10.16 |  78.20 |
---------+--------+--------+
Total       14137     6514    20651
            68.46    31.54   100.00

flag_gene_01vM_3h_p05
          flag_xs_01vM_3h_p05

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  16462 |   1279 |  17741
         |  79.72 |   6.19 |  85.91
         |  92.79 |   7.21 |
         |  94.25 |  40.17 |
---------+--------+--------+
       1 |   1005 |   1905 |   2910
         |   4.87 |   9.22 |  14.09
         |  34.54 |  65.46 |
         |   5.75 |  59.83 |
---------+--------+--------+
Total       17467     3184    20651
            84.58    15.42   100.00

 flag_gene_01vM_24h_p05
           flag_xs_01vM_24h_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  16312 |   1403 |  17715
          |  78.99 |   6.79 |  85.78
          |  92.08 |   7.92 |
          |  94.01 |  42.53 |
 ---------+--------+--------+
        1 |   1040 |   1896 |   2936
          |   5.04 |   9.18 |  14.22
          |  35.42 |  64.58 |
          |   5.99 |  57.47 |
 ---------+--------+--------+
 Total       17352     3299    20651
             84.02    15.98   100.00

 flag_gene_01vM_72h_p05
           flag_xs_01vM_72h_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  15134 |   1316 |  16450
          |  73.28 |   6.37 |  79.66
          |  92.00 |   8.00 |
          |  92.81 |  30.29 |
 ---------+--------+--------+
        1 |   1172 |   3029 |   4201
          |   5.68 |  14.67 |  20.34
          |  27.90 |  72.10 |
          |   7.19 |  69.71 |
 ---------+--------+--------+
 Total       16306     4345    20651
             78.96    21.04   100.00

 flag_gene_1vM_1h_p05
           flag_xs_1vM_1h_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  12700 |   1420 |  14120
          |  61.50 |   6.88 |  68.37
          |  89.94 |  10.06 |
          |  89.84 |  21.80 |
 ---------+--------+--------+
        1 |   1437 |   5094 |   6531
          |   6.96 |  24.67 |  31.63
          |  22.00 |  78.00 |
          |  10.16 |  78.20 |
 ---------+--------+--------+
 Total       14137     6514    20651
             68.46    31.54   100.00


 flag_gene_1vM_3h_p05
           flag_xs_1vM_3h_p05

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  16462 |   1279 |  17741
          |  79.72 |   6.19 |  85.91
          |  92.79 |   7.21 |
          |  94.25 |  40.17 |
 ---------+--------+--------+
        1 |   1005 |   1905 |   2910
          |   4.87 |   9.22 |  14.09
          |  34.54 |  65.46 |
          |   5.75 |  59.83 |
 ---------+--------+--------+
 Total       17467     3184    20651
             84.58    15.42   100.00

flag_gene_1vM_24h_p05
          flag_xs_1vM_24h_p05

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  16312 |   1403 |  17715
         |  78.99 |   6.79 |  85.78
         |  92.08 |   7.92 |
         |  94.01 |  42.53 |
---------+--------+--------+
       1 |   1040 |   1896 |   2936
         |   5.04 |   9.18 |  14.22
         |  35.42 |  64.58 |
         |   5.99 |  57.47 |
---------+--------+--------+
Total       17352     3299    20651
            84.02    15.98   100.00


flag_gene_1vM_72h_p05
          flag_xs_1vM_72h_p05

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  15134 |   1316 |  16450
         |  73.28 |   6.37 |  79.66
         |  92.00 |   8.00 |
         |  92.81 |  30.29 |
---------+--------+--------+
       1 |   1172 |   3029 |   4201
         |   5.68 |  14.67 |  20.34
         |  27.90 |  72.10 |
         |   7.19 |  69.71 |
---------+--------+--------+
Total       16306     4345    20651
            78.96    21.04   100.00



*/
