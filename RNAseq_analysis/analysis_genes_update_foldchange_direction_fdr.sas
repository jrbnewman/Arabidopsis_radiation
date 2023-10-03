ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* update corresponding up/down indicators, given FDR-corrected P values */

data sign;
  set rs.arab_sign_by_contrast_gene;
  keep gene_id sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
       sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h;
run;

data fdr;
  set rs.fdr_by_gene;
run;

proc sort data=sign;
  by gene_id;
proc sort data=fdr;
  by gene_id;
run;

data fdr_w_sign;
  merge fdr (in=in1) sign (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Redo FC flags */

data fdr_w_sign2;
  set fdr_w_sign;
  length sign_01gy_v_Mock_1h_fdr05 $1.;
  length sign_01gy_v_Mock_3h_fdr05 $1.;
  length sign_01gy_v_Mock_24h_fdr05 $1.;
  length sign_01gy_v_Mock_72h_fdr05 $1.;
  length sign_1gy_v_Mock_1h_fdr05 $1.;
  length sign_1gy_v_Mock_3h_fdr05 $1.;
  length sign_1gy_v_Mock_24h_fdr05 $1.;
  length sign_1gy_v_Mock_72h_fdr05 $1.;

  length sign_01gy_v_Mock_1h_fdr10 $1.;
  length sign_01gy_v_Mock_3h_fdr10 $1.;
  length sign_01gy_v_Mock_24h_fdr10 $1.;
  length sign_01gy_v_Mock_72h_fdr10 $1.;
  length sign_1gy_v_Mock_1h_fdr10 $1.;
  length sign_1gy_v_Mock_3h_fdr10 $1.;
  length sign_1gy_v_Mock_24h_fdr10 $1.;
  length sign_1gy_v_Mock_72h_fdr10 $1.;


  length sign_01gy_v_Mock_1h_fdr20 $1.;
  length sign_01gy_v_Mock_3h_fdr20 $1.;
  length sign_01gy_v_Mock_24h_fdr20 $1.;
  length sign_01gy_v_Mock_72h_fdr20 $1.;
  length sign_1gy_v_Mock_1h_fdr20 $1.;
  length sign_1gy_v_Mock_3h_fdr20 $1.;
  length sign_1gy_v_Mock_24h_fdr20 $1.;
  length sign_1gy_v_Mock_72h_fdr20 $1.;

  /* FDR 5% */
  if flag_01gy_v_Mock_1h_fdr05=1 then sign_01gy_v_Mock_1h_fdr05=sign_01gy_v_Mock_1h;
  else sign_01gy_v_Mock_1h_fdr05="N";

  if flag_01gy_v_Mock_3h_fdr05=1 then sign_01gy_v_Mock_3h_fdr05=sign_01gy_v_Mock_3h;
  else sign_01gy_v_Mock_3h_fdr05="N";

  if flag_01gy_v_Mock_24h_fdr05=1 then sign_01gy_v_Mock_24h_fdr05=sign_01gy_v_Mock_24h;
  else sign_01gy_v_Mock_24h_fdr05="N";

  if flag_01gy_v_Mock_72h_fdr05=1 then sign_01gy_v_Mock_72h_fdr05=sign_01gy_v_Mock_72h;
  else sign_01gy_v_Mock_72h_fdr05="N";


  if flag_1gy_v_Mock_1h_fdr05=1 then sign_1gy_v_Mock_1h_fdr05=sign_1gy_v_Mock_1h;
  else sign_1gy_v_Mock_1h_fdr05="N";

  if flag_1gy_v_Mock_3h_fdr05=1 then sign_1gy_v_Mock_3h_fdr05=sign_1gy_v_Mock_3h;
  else sign_1gy_v_Mock_3h_fdr05="N";

  if flag_1gy_v_Mock_24h_fdr05=1 then sign_1gy_v_Mock_24h_fdr05=sign_1gy_v_Mock_24h;
  else sign_1gy_v_Mock_24h_fdr05="N";

  if flag_1gy_v_Mock_72h_fdr05=1 then sign_1gy_v_Mock_72h_fdr05=sign_1gy_v_Mock_72h;
  else sign_1gy_v_Mock_72h_fdr05="N";

  /* FDR 10% */
  if flag_01gy_v_Mock_1h_fdr10=1 then sign_01gy_v_Mock_1h_fdr10=sign_01gy_v_Mock_1h;
  else sign_01gy_v_Mock_1h_fdr10="N";

  if flag_01gy_v_Mock_3h_fdr10=1 then sign_01gy_v_Mock_3h_fdr10=sign_01gy_v_Mock_3h;
  else sign_01gy_v_Mock_3h_fdr10="N";

  if flag_01gy_v_Mock_24h_fdr10=1 then sign_01gy_v_Mock_24h_fdr10=sign_01gy_v_Mock_24h;
  else sign_01gy_v_Mock_24h_fdr10="N";

  if flag_01gy_v_Mock_72h_fdr10=1 then sign_01gy_v_Mock_72h_fdr10=sign_01gy_v_Mock_72h;
  else sign_01gy_v_Mock_72h_fdr10="N";


  if flag_1gy_v_Mock_1h_fdr10=1 then sign_1gy_v_Mock_1h_fdr10=sign_1gy_v_Mock_1h;
  else sign_1gy_v_Mock_1h_fdr10="N";

  if flag_1gy_v_Mock_3h_fdr10=1 then sign_1gy_v_Mock_3h_fdr10=sign_1gy_v_Mock_3h;
  else sign_1gy_v_Mock_3h_fdr10="N";

  if flag_1gy_v_Mock_24h_fdr10=1 then sign_1gy_v_Mock_24h_fdr10=sign_1gy_v_Mock_24h;
  else sign_1gy_v_Mock_24h_fdr10="N";

  if flag_1gy_v_Mock_72h_fdr10=1 then sign_1gy_v_Mock_72h_fdr10=sign_1gy_v_Mock_72h;
  else sign_1gy_v_Mock_72h_fdr10="N";


  /* FDR 20% */
  if flag_01gy_v_Mock_1h_fdr20=1 then sign_01gy_v_Mock_1h_fdr20=sign_01gy_v_Mock_1h;
  else sign_01gy_v_Mock_1h_fdr20="N";

  if flag_01gy_v_Mock_3h_fdr20=1 then sign_01gy_v_Mock_3h_fdr20=sign_01gy_v_Mock_3h;
  else sign_01gy_v_Mock_3h_fdr20="N";

  if flag_01gy_v_Mock_24h_fdr20=1 then sign_01gy_v_Mock_24h_fdr20=sign_01gy_v_Mock_24h;
  else sign_01gy_v_Mock_24h_fdr20="N";

  if flag_01gy_v_Mock_72h_fdr20=1 then sign_01gy_v_Mock_72h_fdr20=sign_01gy_v_Mock_72h;
  else sign_01gy_v_Mock_72h_fdr20="N";


  if flag_1gy_v_Mock_1h_fdr20=1 then sign_1gy_v_Mock_1h_fdr20=sign_1gy_v_Mock_1h;
  else sign_1gy_v_Mock_1h_fdr20="N";

  if flag_1gy_v_Mock_3h_fdr20=1 then sign_1gy_v_Mock_3h_fdr20=sign_1gy_v_Mock_3h;
  else sign_1gy_v_Mock_3h_fdr20="N";

  if flag_1gy_v_Mock_24h_fdr20=1 then sign_1gy_v_Mock_24h_fdr20=sign_1gy_v_Mock_24h;
  else sign_1gy_v_Mock_24h_fdr20="N";

  if flag_1gy_v_Mock_72h_fdr20=1 then sign_1gy_v_Mock_72h_fdr20=sign_1gy_v_Mock_72h;
  else sign_1gy_v_Mock_72h_fdr20="N";

  keep gene_id sign_: ;
  drop sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
       sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h ;
run;

data  rs.arab_sign_by_contrast_gene_fdr;
  set fdr_w_sign2;
run;


proc freq data=fdr_w_sign2;
   tables  sign_01gy_v_Mock_1h_fdr05  sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_24h_fdr05
           sign_01gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_3h_fdr05 
           sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_1h_fdr10
           sign_01gy_v_Mock_3h_fdr10 sign_01gy_v_Mock_24h_fdr10 sign_01gy_v_Mock_72h_fdr10
           sign_1gy_v_Mock_1h_fdr10 sign_1gy_v_Mock_3h_fdr10 sign_1gy_v_Mock_24h_fdr10
           sign_1gy_v_Mock_72h_fdr10 sign_01gy_v_Mock_1h_fdr20 sign_01gy_v_Mock_3h_fdr20
           sign_01gy_v_Mock_24h_fdr20 sign_01gy_v_Mock_72h_fdr20 sign_1gy_v_Mock_1h_fdr20
           sign_1gy_v_Mock_3h_fdr20 sign_1gy_v_Mock_24h_fdr20 sign_1gy_v_Mock_72h_fdr20 ;
 run;

/*
sign_01gy_
v_Mock_1h_                             Cumulative    Cumulative
fdr05         Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
D                 2247       10.34          2247        10.34
N                17523       80.61         19770        90.95
U                 1968        9.05         21738       100.00


sign_01gy_
v_Mock_3h_                             Cumulative    Cumulative
fdr05         Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
D                   76        0.35            76         0.35
N                21613       99.42         21689        99.77
U                   49        0.23         21738       100.00


sign_01gy_
v_Mock_                                Cumulative    Cumulative
24h_fdr05     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
D                  161        0.74           161         0.74
N                21440       98.63         21601        99.37
U                  137        0.63         21738       100.00

 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 72h_fdr05     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 D                  868        3.99           868         3.99
 N                20351       93.62         21219        97.61
 U                  519        2.39         21738       100.00


  sign_1gy_
  v_Mock_                               Cumulative    Cumulative
  1h_fdr05     Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------
  D                1313        6.04          1313         6.04
  N               18943       87.14         20256        93.18
  U                1482        6.82         21738       100.00


  sign_1gy_
  v_Mock_                               Cumulative    Cumulative
  3h_fdr05     Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------
  D                 429        1.97           429         1.97
  N               20995       96.58         21424        98.56
  U                 314        1.44         21738       100.00

   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   24h_fdr05    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   D                  15        0.07            15         0.07
   N               21698       99.82         21713        99.88
   U                  25        0.12         21738       100.00


   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   72h_fdr05    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   D                  16        0.07            16         0.07
   N               21681       99.74         21697        99.81
   U                  41        0.19         21738       100.00


  sign_01gy_
  v_Mock_1h_                             Cumulative    Cumulative
  fdr10         Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                 2919       13.43          2919        13.43
  N                16299       74.98         19218        88.41
  U                 2520       11.59         21738       100.00

 sign_01gy_
 v_Mock_3h_                             Cumulative    Cumulative
 fdr10         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 D                  180        0.83           180         0.83
 N                21422       98.55         21602        99.37
 U                  136        0.63         21738       100.00


 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 24h_fdr10     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 D                  303        1.39           303         1.39
 N                21157       97.33         21460        98.72
 U                  278        1.28         21738       100.00


 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 72h_fdr10     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 D                 1313        6.04          1313         6.04
 N                19534       89.86         20847        95.90
 U                  891        4.10         21738       100.00

sign_1gy_
v_Mock_                               Cumulative    Cumulative
1h_fdr10     Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
D                1763        8.11          1763         8.11
N               18019       82.89         19782        91.00
U                1956        9.00         21738       100.00


sign_1gy_
v_Mock_                               Cumulative    Cumulative
3h_fdr10     Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
D                 838        3.86           838         3.86
N               20194       92.90         21032        96.75
U                 706        3.25         21738       100.00


sign_1gy_
v_Mock_                               Cumulative    Cumulative
24h_fdr10    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
D                  88        0.40            88         0.40
N               21551       99.14         21639        99.54
U                  99        0.46         21738       100.00

   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   72h_fdr10    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   D                  91        0.42            91         0.42
   N               21443       98.64         21534        99.06
   U                 204        0.94         21738       100.00


  sign_01gy_
  v_Mock_1h_                             Cumulative    Cumulative
  fdr20         Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                 3639       16.74          3639        16.74
  N                15047       69.22         18686        85.96
  U                 3052       14.04         21738       100.00


  sign_01gy_
  v_Mock_3h_                             Cumulative    Cumulative
  fdr20         Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                  531        2.44           531         2.44
  N                20795       95.66         21326        98.10
  U                  412        1.90         21738       100.00

  sign_01gy_
  v_Mock_                                Cumulative    Cumulative
  24h_fdr20     Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                  718        3.30           718         3.30
  N                20336       93.55         21054        96.85
  U                  684        3.15         21738       100.00


  sign_01gy_
  v_Mock_                                Cumulative    Cumulative
  72h_fdr20     Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------
  D                 2012        9.26          2012         9.26
  N                18219       83.81         20231        93.07
  U                 1507        6.93         21738       100.00


   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   1h_fdr20     Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   D                2526       11.62          2526        11.62
   N               16548       76.12         19074        87.74
   U                2664       12.26         21738       100.00

 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 3h_fdr20     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 D                1539        7.08          1539         7.08
 N               18649       85.79         20188        92.87
 U                1550        7.13         21738       100.00


 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 24h_fdr20    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 D                 365        1.68           365         1.68
 N               20990       96.56         21355        98.24
 U                 383        1.76         21738       100.00


 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 72h_fdr20    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 D                 466        2.14           466         2.14
 N               20362       93.67         20828        95.81
 U                 910        4.19         21738       100.00


*/

