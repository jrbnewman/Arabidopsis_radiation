ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* update corresponding up/down indicators, given FDR-corrected P values

   If gene is DD, then we want to flag and indicate if Mock or treat only */

data dtct;
   set rs.arab_flag_gene_on_gt0;
   drop flag_gene_on_01gy_apn0 flag_gene_on_1gy_apn0 flag_gene_on_mock_apn0 
        flag_gene_on_1hr_apn0 flag_gene_on_3hr_apn0 flag_gene_on_34hr_apn0 flag_gene_on_72hr_apn0;
run;


data sign;
  set rs.arab_sign_by_contrast_gene;
  keep gene_id sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
       sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h;
run;

data fdr;
  set rs.fdr_by_gene_v2;
run;

proc sort data=sign;
  by gene_id;
proc sort data=dtct;
  by gene_id;
proc sort data=fdr;
  by gene_id;
run;

data fdr_w_sign;
  merge fdr (in=in1) sign (in=in2) dtct (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

/* Redo FC flags */

%macro signFDR(fdrsuffix);

data fdr_w_sign_&fdrsuffix.;
  set fdr_w_sign;
  length sign_01gy_v_Mock_1h_&fdrsuffix. $5.;
  length sign_01gy_v_Mock_3h_&fdrsuffix. $5.;
  length sign_01gy_v_Mock_24h_&fdrsuffix. $5.;
  length sign_01gy_v_Mock_72h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_1h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_3h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_24h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_72h_&fdrsuffix. $5.;

  /* 0.1gy vs Mock */
  if flag_gene_on_01gy_1hr_apn0=1 and flag_gene_on_Mock_1hr_apn0=1 then do;
     if flag_01gy_v_Mock_1h_&fdrsuffix.=1 then sign_01gy_v_Mock_1h_&fdrsuffix.=sign_01gy_v_Mock_1h;
     else sign_01gy_v_Mock_1h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_1hr_apn0=1 and flag_gene_on_Mock_1hr_apn0 ne 1 then sign_01gy_v_Mock_1h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_1hr_apn0 ne 1 and flag_gene_on_Mock_1hr_apn0=1 then sign_01gy_v_Mock_1h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_1h_&fdrsuffix.="";

  if flag_gene_on_01gy_3hr_apn0=1 and flag_gene_on_Mock_3hr_apn0=1 then do;
     if flag_01gy_v_Mock_3h_&fdrsuffix.=1 then sign_01gy_v_Mock_3h_&fdrsuffix.=sign_01gy_v_Mock_3h;
     else sign_01gy_v_Mock_3h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_3hr_apn0=1 and flag_gene_on_Mock_3hr_apn0 ne 1 then sign_01gy_v_Mock_3h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_3hr_apn0 ne 1 and flag_gene_on_Mock_3hr_apn0=1 then sign_01gy_v_Mock_3h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_3h_&fdrsuffix.="";

  if flag_gene_on_01gy_24hr_apn0=1 and flag_gene_on_Mock_24hr_apn0=1 then do;
     if flag_01gy_v_Mock_24h_&fdrsuffix.=1 then sign_01gy_v_Mock_24h_&fdrsuffix.=sign_01gy_v_Mock_24h;
     else sign_01gy_v_Mock_24h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_24hr_apn0=1 and flag_gene_on_Mock_24hr_apn0 ne 1 then sign_01gy_v_Mock_24h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_24hr_apn0 ne 1 and flag_gene_on_Mock_24hr_apn0=1 then sign_01gy_v_Mock_24h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_24h_&fdrsuffix.="";

  if flag_gene_on_01gy_72hr_apn0=1 and flag_gene_on_Mock_72hr_apn0=1 then do;
     if flag_01gy_v_Mock_72h_&fdrsuffix.=1 then sign_01gy_v_Mock_72h_&fdrsuffix.=sign_01gy_v_Mock_72h;
     else sign_01gy_v_Mock_72h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_72hr_apn0=1 and flag_gene_on_Mock_72hr_apn0=0 then sign_01gy_v_Mock_72h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_72hr_apn0 ne 1 and flag_gene_on_Mock_72hr_apn0=1 then sign_01gy_v_Mock_72h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_72h_&fdrsuffix.="";

  /* 1gy vs Mock */
  if flag_gene_on_1gy_1hr_apn0=1 and flag_gene_on_Mock_1hr_apn0=1 then do;
     if flag_1gy_v_Mock_1h_&fdrsuffix.=1 then sign_1gy_v_Mock_1h_&fdrsuffix.=sign_1gy_v_Mock_1h;
     else sign_1gy_v_Mock_1h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_1hr_apn0=1 and flag_gene_on_Mock_1hr_apn0 ne 1 then sign_1gy_v_Mock_1h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_1hr_apn0 ne 1 and flag_gene_on_Mock_1hr_apn0=1 then sign_1gy_v_Mock_1h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_1h_&fdrsuffix.="";

  if flag_gene_on_1gy_3hr_apn0=1 and flag_gene_on_Mock_3hr_apn0=1 then do;
     if flag_1gy_v_Mock_3h_&fdrsuffix.=1 then sign_1gy_v_Mock_3h_&fdrsuffix.=sign_1gy_v_Mock_3h;
     else sign_1gy_v_Mock_3h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_3hr_apn0=1 and flag_gene_on_Mock_3hr_apn0 ne 1 then sign_1gy_v_Mock_3h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_3hr_apn0 ne 1 and flag_gene_on_Mock_3hr_apn0=1 then sign_1gy_v_Mock_3h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_3h_&fdrsuffix.="";

  if flag_gene_on_1gy_24hr_apn0=1 and flag_gene_on_Mock_24hr_apn0=1 then do;
     if flag_1gy_v_Mock_24h_&fdrsuffix.=1 then sign_1gy_v_Mock_24h_&fdrsuffix.=sign_1gy_v_Mock_24h;
     else sign_1gy_v_Mock_24h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_24hr_apn0=1 and flag_gene_on_Mock_24hr_apn0 ne 1 then sign_1gy_v_Mock_24h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_24hr_apn0 ne 1 and flag_gene_on_Mock_24hr_apn0=1 then sign_1gy_v_Mock_24h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_24h_&fdrsuffix.="";

  if flag_gene_on_1gy_72hr_apn0=1 and flag_gene_on_Mock_72hr_apn0=1 then do;
     if flag_1gy_v_Mock_72h_&fdrsuffix.=1 then sign_1gy_v_Mock_72h_&fdrsuffix.=sign_1gy_v_Mock_72h;
     else sign_1gy_v_Mock_72h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_72hr_apn0=1 and flag_gene_on_Mock_72hr_apn0 ne 1 then sign_1gy_v_Mock_72h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_72hr_apn0 ne 1 and flag_gene_on_Mock_72hr_apn0=1 then sign_1gy_v_Mock_72h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_72h_&fdrsuffix.="";

  keep gene_id sign_: flag_gene_on_: ;
  drop sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
       sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h ;
run;

%mend;
%signFDR(fdr05);
%signFDR(fdr10);
%signFDR(fdr20);

proc sort data=fdr_w_sign_fdr05;
   by gene_id;
proc sort data=fdr_w_sign_fdr10;
   by gene_id;
proc sort data=fdr_w_sign_fdr20;
   by gene_id;
run;

data fdr_w_sign2;
  merge fdr_w_sign_fdr05 fdr_w_sign_fdr10 fdr_w_sign_fdr20;
  by gene_id;
  
if flag_gene_on_01gy_1hr_apn0=1 and flag_gene_on_mock_1hr_apn0=1 then flag_gene_01gy_v_Mock_1h_dd=0;
else if flag_gene_on_01gy_1hr_apn0=1 or flag_gene_on_mock_1hr_apn0=1 then flag_gene_01gy_v_Mock_1h_dd=1;
else flag_gene_01gy_v_Mock_1h_dd=0;

if flag_gene_on_01gy_3hr_apn0=1 and flag_gene_on_mock_3hr_apn0=1 then flag_gene_01gy_v_Mock_3h_dd=0;
else if flag_gene_on_01gy_3hr_apn0=1 or flag_gene_on_mock_3hr_apn0=1 then flag_gene_01gy_v_Mock_3h_dd=1;
else flag_gene_01gy_v_Mock_3h_dd=0;

if flag_gene_on_01gy_24hr_apn0=1 and flag_gene_on_mock_24hr_apn0=1 then flag_gene_01gy_v_Mock_24h_dd=0;
else if flag_gene_on_01gy_24hr_apn0=1 or flag_gene_on_mock_24hr_apn0=1 then flag_gene_01gy_v_Mock_24h_dd=1;
else flag_gene_01gy_v_Mock_24h_dd=0;

if flag_gene_on_01gy_72hr_apn0=1 and flag_gene_on_mock_72hr_apn0=1 then flag_gene_01gy_v_Mock_72h_dd=0;
else if flag_gene_on_01gy_72hr_apn0=1 or flag_gene_on_mock_72hr_apn0=1 then flag_gene_01gy_v_Mock_72h_dd=1;
else flag_gene_01gy_v_Mock_72h_dd=0;


if flag_gene_on_1gy_1hr_apn0=1 and flag_gene_on_mock_1hr_apn0=1 then flag_gene_1gy_v_Mock_1h_dd=0;
else if flag_gene_on_1gy_1hr_apn0=1 or flag_gene_on_mock_1hr_apn0=1 then flag_gene_1gy_v_Mock_1h_dd=1;
else flag_gene_1gy_v_Mock_1h_dd=0;

if flag_gene_on_1gy_3hr_apn0=1 and flag_gene_on_mock_3hr_apn0=1 then flag_gene_1gy_v_Mock_3h_dd=0;
else if flag_gene_on_1gy_3hr_apn0=1 or flag_gene_on_mock_3hr_apn0=1 then flag_gene_1gy_v_Mock_3h_dd=1;
else flag_gene_1gy_v_Mock_3h_dd=0;

if flag_gene_on_1gy_24hr_apn0=1 and flag_gene_on_mock_24hr_apn0=1 then flag_gene_1gy_v_Mock_24h_dd=0;
else if flag_gene_on_1gy_24hr_apn0=1 or flag_gene_on_mock_24hr_apn0=1 then flag_gene_1gy_v_Mock_24h_dd=1;
else flag_gene_1gy_v_Mock_24h_dd=0;

if flag_gene_on_1gy_72hr_apn0=1 and flag_gene_on_mock_72hr_apn0=1 then flag_gene_1gy_v_Mock_72h_dd=0;
else if flag_gene_on_1gy_72hr_apn0=1 or flag_gene_on_mock_72hr_apn0=1 then flag_gene_1gy_v_Mock_72h_dd=1;
else flag_gene_1gy_v_Mock_72h_dd=0;

run;


data  rs.arab_sign_by_contrast_gene_fdr2;
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
           sign_1gy_v_Mock_3h_fdr20 sign_1gy_v_Mock_24h_fdr20 sign_1gy_v_Mock_72h_fdr20 
	   flag_gene_01gy_v_Mock_1h_dd flag_gene_01gy_v_Mock_3h_dd 
           flag_gene_01gy_v_Mock_24h_dd flag_gene_01gy_v_Mock_72h_dd
	   flag_gene_1gy_v_Mock_1h_dd flag_gene_1gy_v_Mock_3h_dd 
           flag_gene_1gy_v_Mock_24h_dd flag_gene_1gy_v_Mock_72h_dd

	   flag_gene_01gy_v_Mock_1h_dd*sign_01gy_v_Mock_1h_fdr05 flag_gene_01gy_v_Mock_3h_dd*sign_01gy_v_Mock_3h_fdr05
           flag_gene_01gy_v_Mock_24h_dd*sign_01gy_v_Mock_24h_fdr05 flag_gene_01gy_v_Mock_72h_dd*sign_01gy_v_Mock_72h_fdr05
	   flag_gene_1gy_v_Mock_1h_dd*sign_1gy_v_Mock_1h_fdr05 flag_gene_1gy_v_Mock_3h_dd*sign_1gy_v_Mock_3h_fdr05
           flag_gene_1gy_v_Mock_24h_dd*sign_1gy_v_Mock_24h_fdr05 flag_gene_1gy_v_Mock_72h_dd*sign_1gy_v_Mock_72h_fdr05
           ;
 run;

/*
 sign_01gy_
 v_Mock_1h_                             Cumulative    Cumulative
 fdr05         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy              353        1.65           353         1.65
 D                 2275       10.64          2628        12.30
 Mock               420        1.97          3048        14.26
 N                16364       76.57         19412        90.83
 U                 1960        9.17         21372       100.00

                     Frequency Missing = 366

sign_01gy_
v_Mock_3h_                             Cumulative    Cumulative
fdr05         Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              366        1.71           366         1.71
D                   79        0.37           445         2.08
Mock               348        1.63           793         3.71
N                20555       96.06         21348        99.77
U                   50        0.23         21398       100.00

                    Frequency Missing = 340

sign_01gy_
v_Mock_                                Cumulative    Cumulative
24h_fdr05     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              377        1.76           377         1.76
D                  170        0.79           547         2.56
Mock               317        1.48           864         4.04
N                20399       95.31         21263        99.35
U                  140        0.65         21403       100.00

                    Frequency Missing = 335

 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 72h_fdr05     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy               97        0.46            97         0.46
 D                  889        4.21           986         4.67
 Mock               309        1.46          1295         6.14
 N                19261       91.31         20556        97.45
 U                  538        2.55         21094       100.00

                     Frequency Missing = 644

 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 1h_fdr05     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy               371        1.73           371         1.73
 D                1330        6.22          1701         7.95
 Mock              343        1.60          2044         9.56
 N               17845       83.43         19889        92.98
 U                1501        7.02         21390       100.00

                    Frequency Missing = 348

 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 3h_fdr05     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy               391        1.83           391         1.83
 D                 453        2.11           844         3.94
 Mock              285        1.33          1129         5.27
 N               19971       93.22         21100        98.49
 U                 323        1.51         21423       100.00

                    Frequency Missing = 315
sign_1gy_
v_Mock_                               Cumulative    Cumulative
24h_fdr05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               409        1.91           409         1.91
D                  14        0.07           423         1.97
Mock              293        1.37           716         3.34
N               20694       96.54         21410        99.88
U                  25        0.12         21435       100.00

                   Frequency Missing = 303

sign_1gy_
v_Mock_                               Cumulative    Cumulative
72h_fdr05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               472        2.20           472         2.20
D                  15        0.07           487         2.27
Mock              241        1.12           728         3.39
N               20703       96.43         21431        99.82
U                  38        0.18         21469       100.00

                   Frequency Missing = 269
sign_01gy_
v_Mock_1h_                             Cumulative    Cumulative
fdr10         Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              353        1.65           353         1.65
D                 2943       13.77          3296        15.42
Mock               420        1.97          3716        17.39
N                15139       70.84         18855        88.22
U                 2517       11.78         21372       100.00

                    Frequency Missing = 366
sign_01gy_
v_Mock_3h_                             Cumulative    Cumulative
fdr10         Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              366        1.71           366         1.71
D                  186        0.87           552         2.58
Mock               348        1.63           900         4.21
N                20362       95.16         21262        99.36
U                  136        0.64         21398       100.00

                    Frequency Missing = 340
sign_01gy_
v_Mock_                                Cumulative    Cumulative
24h_fdr10     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              377        1.76           377         1.76
D                  302        1.41           679         3.17
Mock               317        1.48           996         4.65
N                20138       94.09         21134        98.74
U                  269        1.26         21403       100.00

                    Frequency Missing = 335

sign_01gy_
v_Mock_                                Cumulative    Cumulative
72h_fdr10     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy               97        0.46            97         0.46
D                 1337        6.34          1434         6.80
Mock               309        1.46          1743         8.26
N                18448       87.46         20191        95.72
U                  903        4.28         21094       100.00

                    Frequency Missing = 644

sign_1gy_
v_Mock_                               Cumulative    Cumulative
1h_fdr10     Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               371        1.73           371         1.73
D                1789        8.36          2160        10.10
Mock              343        1.60          2503        11.70
N               16931       79.15         19434        90.86
U                1956        9.14         21390       100.00

                   Frequency Missing = 348

sign_1gy_
v_Mock_                               Cumulative    Cumulative
3h_fdr10     Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               391        1.83           391         1.83
D                 866        4.04          1257         5.87
Mock              285        1.33          1542         7.20
N               19150       89.39         20692        96.59
U                 731        3.41         21423       100.00

                   Frequency Missing = 315

sign_1gy_
v_Mock_                               Cumulative    Cumulative
24h_fdr10    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               409        1.91           409         1.91
D                  86        0.40           495         2.31
Mock              293        1.37           788         3.68
N               20552       95.88         21340        99.56
U                  95        0.44         21435       100.00

                   Frequency Missing = 303

sign_1gy_
v_Mock_                               Cumulative    Cumulative
72h_fdr10    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               472        2.20           472         2.20
D                  87        0.41           559         2.60
Mock              241        1.12           800         3.73
N               20480       95.39         21280        99.12
U                 189        0.88         21469       100.00

                   Frequency Missing = 269

sign_01gy_
v_Mock_1h_                             Cumulative    Cumulative
fdr20         Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              353        1.65           353         1.65
D                 3586       16.78          3939        18.43
Mock               420        1.97          4359        20.40
N                14026       65.63         18385        86.02
U                 2987       13.98         21372       100.00

                    Frequency Missing = 366

 sign_01gy_
 v_Mock_3h_                             Cumulative    Cumulative
 fdr20         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy              366        1.71           366         1.71
 D                  569        2.66           935         4.37
 Mock               348        1.63          1283         6.00
 N                19676       91.95         20959        97.95
 U                  439        2.05         21398       100.00

                     Frequency Missing = 340

sign_01gy_
v_Mock_                                Cumulative    Cumulative
24h_fdr20     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              377        1.76           377         1.76
D                  727        3.40          1104         5.16
Mock               317        1.48          1421         6.64
N                19300       90.17         20721        96.81
U                  682        3.19         21403       100.00

                    Frequency Missing = 335

sign_01gy_
v_Mock_                                Cumulative    Cumulative
72h_fdr20     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy               97        0.46            97         0.46
D                 2046        9.70          2143        10.16
Mock               309        1.46          2452        11.62
N                17120       81.16         19572        92.78
U                 1522        7.22         21094       100.00

                    Frequency Missing = 644

 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 1h_fdr20     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy               371        1.73           371         1.73
 D                2542       11.88          2913        13.62
 Mock              343        1.60          3256        15.22
 N               15483       72.38         18739        87.61
 U                2651       12.39         21390       100.00

                    Frequency Missing = 348

sign_1gy_
v_Mock_                               Cumulative    Cumulative
3h_fdr20     Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               391        1.83           391         1.83
D                1572        7.34          1963         9.16
Mock              285        1.33          2248        10.49
N               17606       82.18         19854        92.68
U                1569        7.32         21423       100.00

                   Frequency Missing = 315

sign_1gy_
v_Mock_                               Cumulative    Cumulative
24h_fdr20    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               409        1.91           409         1.91
D                 378        1.76           787         3.67
Mock              293        1.37          1080         5.04
N               19963       93.13         21043        98.17
U                 392        1.83         21435       100.00

                   Frequency Missing = 303

   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   72h_fdr20    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   1gy               472        2.20           472         2.20
   D                 459        2.14           931         4.34
   Mock              241        1.12          1172         5.46
   N               19438       90.54         20610        96.00
   U                 859        4.00         21469       100.00

                      Frequency Missing = 269


flag_gene_01gy_                             Cumulative    Cumulative
   v_Mock_1h_dd    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       20965       96.44         20965        96.44
              1         773        3.56         21738       100.00

  flag_gene_01gy_                             Cumulative    Cumulative
     v_Mock_3h_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       21024       96.72         21024        96.72
                1         714        3.28         21738       100.00


  flag_gene_01gy_                             Cumulative    Cumulative
    v_Mock_24h_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       21044       96.81         21044        96.81
                1         694        3.19         21738       100.00


  flag_gene_01gy_                             Cumulative    Cumulative
    v_Mock_72h_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       21051       96.84         21051        96.84
                1         687        3.16         21738       100.00

  flag_gene_1gy_                             Cumulative    Cumulative
    v_Mock_1h_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       21024       96.72         21024        96.72
               1         714        3.28         21738       100.00


  flag_gene_1gy_                             Cumulative    Cumulative
    v_Mock_3h_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       21062       96.89         21062        96.89
               1         676        3.11         21738       100.00


  flag_gene_1gy_                             Cumulative    Cumulative
   v_Mock_24h_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       21036       96.77         21036        96.77
               1         702        3.23         21738       100.00

  flag_gene_1gy_                             Cumulative    Cumulative
   v_Mock_72h_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       21025       96.72         21025        96.72
               1         713        3.28         21738       100.00
*/

















