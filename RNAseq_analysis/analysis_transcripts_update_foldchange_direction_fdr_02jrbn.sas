ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* update corresponding up/down indicators, given FDR-corrected P values

   If gene is DD, then we want to flag and indicate if Mock or treat only */

data dtct;
   set rs.arab_flag_transcript_on_gt0;
   drop flag_xscript_on_01gy_apn0 flag_xscript_on_1gy_apn0 flag_xscript_on_mock_apn0 
        flag_xscript_on_1hr_tpm0 flag_xscript_on_3hr_tpm0 flag_xscript_on_24hr_tpm0
        flag_xscript_on_72hr_tpm0;
run;


data sign;
  set rs.arab_sign_by_contrast_transcript;
  keep transcript_id sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
       sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h;
run;

data fdr;
  set rs.fdr_by_transcript;
run;

proc sort data=sign;
  by transcript_id;
proc sort data=dtct;
  by transcript_id;
proc sort data=fdr;
  by transcript_id;
run;

data fdr_w_sign;
  merge fdr (in=in1) sign (in=in2) dtct (in=in3);
  by transcript_id;
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
  if flag_xscript_on_01gy_1hr_tpm0=1 and flag_xscript_on_Mock_1hr_tpm0=1 then do;
     if flag_01gy_v_Mock_1h_&fdrsuffix.=1 then sign_01gy_v_Mock_1h_&fdrsuffix.=sign_01gy_v_Mock_1h;
     else sign_01gy_v_Mock_1h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_01gy_1hr_tpm0=1 and flag_xscript_on_Mock_1hr_tpm0 ne 1 then sign_01gy_v_Mock_1h_&fdrsuffix.="0.1gy";
  else if flag_xscript_on_01gy_1hr_tpm0 ne 1 and flag_xscript_on_Mock_1hr_tpm0=1 then sign_01gy_v_Mock_1h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_1h_&fdrsuffix.="";

  if flag_xscript_on_01gy_3hr_tpm0=1 and flag_xscript_on_Mock_3hr_tpm0=1 then do;
     if flag_01gy_v_Mock_3h_&fdrsuffix.=1 then sign_01gy_v_Mock_3h_&fdrsuffix.=sign_01gy_v_Mock_3h;
     else sign_01gy_v_Mock_3h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_01gy_3hr_tpm0=1 and flag_xscript_on_Mock_3hr_tpm0 ne 1 then sign_01gy_v_Mock_3h_&fdrsuffix.="0.1gy";
  else if flag_xscript_on_01gy_3hr_tpm0 ne 1 and flag_xscript_on_Mock_3hr_tpm0=1 then sign_01gy_v_Mock_3h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_3h_&fdrsuffix.="";

  if flag_xscript_on_01gy_24hr_tpm0=1 and flag_xscript_on_Mock_24hr_tpm0=1 then do;
     if flag_01gy_v_Mock_24h_&fdrsuffix.=1 then sign_01gy_v_Mock_24h_&fdrsuffix.=sign_01gy_v_Mock_24h;
     else sign_01gy_v_Mock_24h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_01gy_24hr_tpm0=1 and flag_xscript_on_Mock_24hr_tpm0 ne 1 then sign_01gy_v_Mock_24h_&fdrsuffix.="0.1gy";
  else if flag_xscript_on_01gy_24hr_tpm0 ne 1 and flag_xscript_on_Mock_24hr_tpm0=1 then sign_01gy_v_Mock_24h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_24h_&fdrsuffix.="";

  if flag_xscript_on_01gy_72hr_tpm0=1 and flag_xscript_on_Mock_72hr_tpm0=1 then do;
     if flag_01gy_v_Mock_72h_&fdrsuffix.=1 then sign_01gy_v_Mock_72h_&fdrsuffix.=sign_01gy_v_Mock_72h;
     else sign_01gy_v_Mock_72h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_01gy_72hr_tpm0=1 and flag_xscript_on_Mock_72hr_tpm0=0 then sign_01gy_v_Mock_72h_&fdrsuffix.="0.1gy";
  else if flag_xscript_on_01gy_72hr_tpm0 ne 1 and flag_xscript_on_Mock_72hr_tpm0=1 then sign_01gy_v_Mock_72h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_72h_&fdrsuffix.="";

  /* 1gy vs Mock */
  if flag_xscript_on_1gy_1hr_tpm0=1 and flag_xscript_on_Mock_1hr_tpm0=1 then do;
     if flag_1gy_v_Mock_1h_&fdrsuffix.=1 then sign_1gy_v_Mock_1h_&fdrsuffix.=sign_1gy_v_Mock_1h;
     else sign_1gy_v_Mock_1h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_1gy_1hr_tpm0=1 and flag_xscript_on_Mock_1hr_tpm0 ne 1 then sign_1gy_v_Mock_1h_&fdrsuffix.="1gy";
  else if flag_xscript_on_1gy_1hr_tpm0 ne 1 and flag_xscript_on_Mock_1hr_tpm0=1 then sign_1gy_v_Mock_1h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_1h_&fdrsuffix.="";

  if flag_xscript_on_1gy_3hr_tpm0=1 and flag_xscript_on_Mock_3hr_tpm0=1 then do;
     if flag_1gy_v_Mock_3h_&fdrsuffix.=1 then sign_1gy_v_Mock_3h_&fdrsuffix.=sign_1gy_v_Mock_3h;
     else sign_1gy_v_Mock_3h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_1gy_3hr_tpm0=1 and flag_xscript_on_Mock_3hr_tpm0 ne 1 then sign_1gy_v_Mock_3h_&fdrsuffix.="1gy";
  else if flag_xscript_on_1gy_3hr_tpm0 ne 1 and flag_xscript_on_Mock_3hr_tpm0=1 then sign_1gy_v_Mock_3h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_3h_&fdrsuffix.="";

  if flag_xscript_on_1gy_24hr_tpm0=1 and flag_xscript_on_Mock_24hr_tpm0=1 then do;
     if flag_1gy_v_Mock_24h_&fdrsuffix.=1 then sign_1gy_v_Mock_24h_&fdrsuffix.=sign_1gy_v_Mock_24h;
     else sign_1gy_v_Mock_24h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_1gy_24hr_tpm0=1 and flag_xscript_on_Mock_24hr_tpm0 ne 1 then sign_1gy_v_Mock_24h_&fdrsuffix.="1gy";
  else if flag_xscript_on_1gy_24hr_tpm0 ne 1 and flag_xscript_on_Mock_24hr_tpm0=1 then sign_1gy_v_Mock_24h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_24h_&fdrsuffix.="";

  if flag_xscript_on_1gy_72hr_tpm0=1 and flag_xscript_on_Mock_72hr_tpm0=1 then do;
     if flag_1gy_v_Mock_72h_&fdrsuffix.=1 then sign_1gy_v_Mock_72h_&fdrsuffix.=sign_1gy_v_Mock_72h;
     else sign_1gy_v_Mock_72h_&fdrsuffix.="N";
     end;
  else if flag_xscript_on_1gy_72hr_tpm0=1 and flag_xscript_on_Mock_72hr_tpm0 ne 1 then sign_1gy_v_Mock_72h_&fdrsuffix.="1gy";
  else if flag_xscript_on_1gy_72hr_tpm0 ne 1 and flag_xscript_on_Mock_72hr_tpm0=1 then sign_1gy_v_Mock_72h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_72h_&fdrsuffix.="";

  keep transcript_id sign_: flag_xscript_on_: ;
  drop sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
       sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h ;
run;

%mend;
%signFDR(fdr05);
%signFDR(fdr10);
%signFDR(fdr20);

proc sort data=fdr_w_sign_fdr05;
   by transcript_id;
proc sort data=fdr_w_sign_fdr10;
   by transcript_id;
proc sort data=fdr_w_sign_fdr20;
   by transcript_id;
run;

data fdr_w_sign2;
  merge fdr_w_sign_fdr05 fdr_w_sign_fdr10 fdr_w_sign_fdr20;
  by transcript_id;
  
if flag_xscript_on_01gy_1hr_tpm0=1 and flag_xscript_on_mock_1hr_tpm0=1 then flag_xscript_01gy_v_Mock_1h_dd=0;
else if flag_xscript_on_01gy_1hr_tpm0=1 or flag_xscript_on_mock_1hr_tpm0=1 then flag_xscript_01gy_v_Mock_1h_dd=1;
else flag_xscript_01gy_v_Mock_1h_dd=0;

if flag_xscript_on_01gy_3hr_tpm0=1 and flag_xscript_on_mock_3hr_tpm0=1 then flag_xscript_01gy_v_Mock_3h_dd=0;
else if flag_xscript_on_01gy_3hr_tpm0=1 or flag_xscript_on_mock_3hr_tpm0=1 then flag_xscript_01gy_v_Mock_3h_dd=1;
else flag_xscript_01gy_v_Mock_3h_dd=0;

if flag_xscript_on_01gy_24hr_tpm0=1 and flag_xscript_on_mock_24hr_tpm0=1 then flag_xscript_01gy_v_Mock_24h_dd=0;
else if flag_xscript_on_01gy_24hr_tpm0=1 or flag_xscript_on_mock_24hr_tpm0=1 then flag_xscript_01gy_v_Mock_24h_dd=1;
else flag_xscript_01gy_v_Mock_24h_dd=0;

if flag_xscript_on_01gy_72hr_tpm0=1 and flag_xscript_on_mock_72hr_tpm0=1 then flag_xscript_01gy_v_Mock_72h_dd=0;
else if flag_xscript_on_01gy_72hr_tpm0=1 or flag_xscript_on_mock_72hr_tpm0=1 then flag_xscript_01gy_v_Mock_72h_dd=1;
else flag_xscript_01gy_v_Mock_72h_dd=0;


if flag_xscript_on_1gy_1hr_tpm0=1 and flag_xscript_on_mock_1hr_tpm0=1 then flag_xscript_1gy_v_Mock_1h_dd=0;
else if flag_xscript_on_1gy_1hr_tpm0=1 or flag_xscript_on_mock_1hr_tpm0=1 then flag_xscript_1gy_v_Mock_1h_dd=1;
else flag_xscript_1gy_v_Mock_1h_dd=0;

if flag_xscript_on_1gy_3hr_tpm0=1 and flag_xscript_on_mock_3hr_tpm0=1 then flag_xscript_1gy_v_Mock_3h_dd=0;
else if flag_xscript_on_1gy_3hr_tpm0=1 or flag_xscript_on_mock_3hr_tpm0=1 then flag_xscript_1gy_v_Mock_3h_dd=1;
else flag_xscript_1gy_v_Mock_3h_dd=0;

if flag_xscript_on_1gy_24hr_tpm0=1 and flag_xscript_on_mock_24hr_tpm0=1 then flag_xscript_1gy_v_Mock_24h_dd=0;
else if flag_xscript_on_1gy_24hr_tpm0=1 or flag_xscript_on_mock_24hr_tpm0=1 then flag_xscript_1gy_v_Mock_24h_dd=1;
else flag_xscript_1gy_v_Mock_24h_dd=0;

if flag_xscript_on_1gy_72hr_tpm0=1 and flag_xscript_on_mock_72hr_tpm0=1 then flag_xscript_1gy_v_Mock_72h_dd=0;
else if flag_xscript_on_1gy_72hr_tpm0=1 or flag_xscript_on_mock_72hr_tpm0=1 then flag_xscript_1gy_v_Mock_72h_dd=1;
else flag_xscript_1gy_v_Mock_72h_dd=0;

run;


data  rs.arab_sign_by_contrast_xs_fdr;
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
	   flag_xscript_01gy_v_Mock_1h_dd flag_xscript_01gy_v_Mock_3h_dd 
           flag_xscript_01gy_v_Mock_24h_dd flag_xscript_01gy_v_Mock_72h_dd
	   flag_xscript_1gy_v_Mock_1h_dd flag_xscript_1gy_v_Mock_3h_dd 
           flag_xscript_1gy_v_Mock_24h_dd flag_xscript_1gy_v_Mock_72h_dd

	   flag_xscript_01gy_v_Mock_1h_dd*sign_01gy_v_Mock_1h_fdr05 flag_xscript_01gy_v_Mock_3h_dd*sign_01gy_v_Mock_3h_fdr05
           flag_xscript_01gy_v_Mock_24h_dd*sign_01gy_v_Mock_24h_fdr05 flag_xscript_01gy_v_Mock_72h_dd*sign_01gy_v_Mock_72h_fdr05
	   flag_xscript_1gy_v_Mock_1h_dd*sign_1gy_v_Mock_1h_fdr05 flag_xscript_1gy_v_Mock_3h_dd*sign_1gy_v_Mock_3h_fdr05
           flag_xscript_1gy_v_Mock_24h_dd*sign_1gy_v_Mock_24h_fdr05 flag_xscript_1gy_v_Mock_72h_dd*sign_1gy_v_Mock_72h_fdr05
           ;
 run;

/*


 sign_01gy_
 v_Mock_1h_                             Cumulative    Cumulative
 fdr05         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1412        3.81          1412         3.81
 D                 2137        5.76          3549         9.56
 Mock              1739        4.69          5288        14.25
 N                30223       81.45         35511        95.70
 U                 1595        4.30         37106       100.00

                    Frequency Missing = 1389


 sign_01gy_
 v_Mock_3h_                             Cumulative    Cumulative
 fdr05         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1434        3.86          1434         3.86
 D                   30        0.08          1464         3.94
 Mock              1478        3.97          2942         7.91
 N                34215       91.99         37157        99.90
 U                   37        0.10         37194       100.00

                    Frequency Missing = 1301


 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 24h_fdr05     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1561        4.20          1561         4.20
 D                  174        0.47          1735         4.67
 Mock              1484        3.99          3219         8.66
 N                33879       91.17         37098        99.83
 U                   64        0.17         37162       100.00

                    Frequency Missing = 1333


 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 72h_fdr05     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy              313        0.87           313         0.87
 D                  856        2.37          1169         3.24
 Mock              1466        4.06          2635         7.29
 N                33210       91.92         35845        99.21
 U                  286        0.79         36131       100.00

                    Frequency Missing = 2364

   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   1h_fdr05     Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   1gy              1509        4.06          1509         4.06
   D                1267        3.41          2776         7.46
   Mock             1484        3.99          4260        11.45
   N               31749       85.34         36009        96.79
   U                1194        3.21         37203       100.00

                      Frequency Missing = 1292


   sign_1gy_
   v_Mock_                               Cumulative    Cumulative
   3h_fdr05     Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------
   1gy              1508        4.05          1508         4.05
   D                 512        1.37          2020         5.42
   Mock             1378        3.70          3398         9.12
   N               33677       90.36         37075        99.48
   U                 193        0.52         37268       100.00

                      Frequency Missing = 1227
 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 24h_fdr05    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy              1571        4.23          1571         4.23
 D                   2        0.01          1573         4.23
 Mock             1298        3.49          2871         7.72
 N               34296       92.26         37167        99.99
 U                   5        0.01         37172       100.00

                    Frequency Missing = 1323


 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 72h_fdr05    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy              1677        4.47          1677         4.47
 D                   6        0.02          1683         4.49
 Mock             1217        3.25          2900         7.73
 N               34592       92.26         37492        99.99
 U                   3        0.01         37495       100.00

                    Frequency Missing = 1000
 sign_01gy_
 v_Mock_1h_                             Cumulative    Cumulative
 fdr10         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1412        3.81          1412         3.81
 D                 2865        7.72          4277        11.53
 Mock              1739        4.69          6016        16.21
 N                29057       78.31         35073        94.52
 U                 2033        5.48         37106       100.00

                    Frequency Missing = 1389


 sign_01gy_
 v_Mock_3h_                             Cumulative    Cumulative
 fdr10         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1434        3.86          1434         3.86
 D                   80        0.22          1514         4.07
 Mock              1478        3.97          2992         8.04
 N                34120       91.74         37112        99.78
 U                   82        0.22         37194       100.00

                    Frequency Missing = 1301

sign_01gy_
v_Mock_                                Cumulative    Cumulative
24h_fdr10     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy             1561        4.20          1561         4.20
D                  333        0.90          1894         5.10
Mock              1484        3.99          3378         9.09
N                33653       90.56         37031        99.65
U                  131        0.35         37162       100.00

                   Frequency Missing = 1333


sign_01gy_
v_Mock_                                Cumulative    Cumulative
72h_fdr10     Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              313        0.87           313         0.87
D                 1288        3.56          1601         4.43
Mock              1466        4.06          3067         8.49
N                32566       90.13         35633        98.62
U                  498        1.38         36131       100.00

                   Frequency Missing = 2364

  sign_1gy_
  v_Mock_                               Cumulative    Cumulative
  1h_fdr10     Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------
  1gy              1509        4.06          1509         4.06
  D                1790        4.81          3299         8.87
  Mock             1484        3.99          4783        12.86
  N               30822       82.85         35605        95.70
  U                1598        4.30         37203       100.00

                     Frequency Missing = 1292


  sign_1gy_
  v_Mock_                               Cumulative    Cumulative
  3h_fdr10     Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------
  1gy              1508        4.05          1508         4.05
  D                 908        2.44          2416         6.48
  Mock             1378        3.70          3794        10.18
  N               33042       88.66         36836        98.84
  U                 432        1.16         37268       100.00

                     Frequency Missing = 1227

  sign_1gy_
  v_Mock_                               Cumulative    Cumulative
  1h_fdr10     Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------
  1gy              1509        4.06          1509         4.06
  D                1790        4.81          3299         8.87
  Mock             1484        3.99          4783        12.86
  N               30822       82.85         35605        95.70
  U                1598        4.30         37203       100.00

                     Frequency Missing = 1292


  sign_1gy_
  v_Mock_                               Cumulative    Cumulative
  3h_fdr10     Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------
  1gy              1508        4.05          1508         4.05
  D                 908        2.44          2416         6.48
  Mock             1378        3.70          3794        10.18
  N               33042       88.66         36836        98.84
  U                 432        1.16         37268       100.00

                     Frequency Missing = 1227

 sign_01gy_
 v_Mock_1h_                             Cumulative    Cumulative
 fdr20         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1412        3.81          1412         3.81
 D                 4129       11.13          5541        14.93
 Mock              1739        4.69          7280        19.62
 N                26960       72.66         34240        92.28
 U                 2866        7.72         37106       100.00

                    Frequency Missing = 1389


 sign_01gy_
 v_Mock_3h_                             Cumulative    Cumulative
 fdr20         Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1434        3.86          1434         3.86
 D                  328        0.88          1762         4.74
 Mock              1478        3.97          3240         8.71
 N                33659       90.50         36899        99.21
 U                  295        0.79         37194       100.00

                    Frequency Missing = 1301

                          The SAS System          11:17 Thursday, A

 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 24h_fdr20     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy             1561        4.20          1561         4.20
 D                  725        1.95          2286         6.15
 Mock              1484        3.99          3770        10.14
 N                33037       88.90         36807        99.04
 U                  355        0.96         37162       100.00

                    Frequency Missing = 1333


 sign_01gy_
 v_Mock_                                Cumulative    Cumulative
 72h_fdr20     Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------
 0.1gy              313        0.87           313         0.87
 D                 1945        5.38          2258         6.25
 Mock              1466        4.06          3724        10.31
 N                31451       87.05         35175        97.35
 U                  956        2.65         36131       100.00

                    Frequency Missing = 2364.

 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 1h_fdr20     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy              1509        4.06          1509         4.06
 D                2635        7.08          4144        11.14
 Mock             1484        3.99          5628        15.13
 N               29355       78.90         34983        94.03
 U                2220        5.97         37203       100.00

                    Frequency Missing = 1292


 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 3h_fdr20     Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy              1508        4.05          1508         4.05
 D                1665        4.47          3173         8.51
 Mock             1378        3.70          4551        12.21
 N               31686       85.02         36237        97.23
 U                1031        2.77         37268       100.00

                    Frequency Missing = 1227

                         The SAS System          11:17 Thursday, Apri

 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 24h_fdr20    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy              1571        4.23          1571         4.23
 D                  98        0.26          1669         4.49
 Mock             1298        3.49          2967         7.98
 N               34112       91.77         37079        99.75
 U                  93        0.25         37172       100.00

                    Frequency Missing = 1323


 sign_1gy_
 v_Mock_                               Cumulative    Cumulative
 72h_fdr20    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------
 1gy              1677        4.47          1677         4.47
 D                 223        0.59          1900         5.07
 Mock             1217        3.25          3117         8.31
 N               33999       90.68         37116        98.99
 U                 379        1.01         37495       100.00

                    Frequency Missing = 1000


     flag_xscript_
      01gy_v_Mock_                             Cumulative    Cumulative
             1h_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       35344       91.81         35344        91.81
                 1        3151        8.19         38495       100.00


     flag_xscript_
      01gy_v_Mock_                             Cumulative    Cumulative
             3h_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       35583       92.44         35583        92.44
                 1        2912        7.56         38495       100.00


     flag_xscript_
      01gy_v_Mock_                             Cumulative    Cumulative
            24h_dd    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       35450       92.09         35450        92.09
                 1        3045        7.91         38495       100.00


    flag_xscript_
     01gy_v_Mock_                             Cumulative    Cumulative
           72h_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       35590       92.45         35590        92.45
                1        2905        7.55         38495       100.00


    flag_xscript_                             Cumulative    Cumulative
 1gy_v_Mock_1h_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       35502       92.22         35502        92.22
                1        2993        7.78         38495       100.00


    flag_xscript_                             Cumulative    Cumulative
 1gy_v_Mock_3h_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       35609       92.50         35609        92.50
                1        2886        7.50         38495       100.00

   flag_xscript_
     1gy_v_Mock_                             Cumulative    Cumulative
          24h_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       35626       92.55         35626        92.55
               1        2869        7.45         38495       100.00


   flag_xscript_
     1gy_v_Mock_                             Cumulative    Cumulative
          72h_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       35601       92.48         35601        92.48
               1        2894        7.52         38495       100.00



     Table of flag_xscript_01gy_v_Mock_1h_dd by sign_01gy_v_Mock_1h_fdr05

        flag_xscript_01gy_v_Mock_1h_dd     sign_01gy_v_Mock_1h_fdr05

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |0.1gy   |D       |Mock    |N       |U       |  Total
        ---------+--------+--------+--------+--------+--------+
               0 |      0 |   2137 |      0 |  30223 |   1595 |  33955
                 |   0.00 |   5.76 |   0.00 |  81.45 |   4.30 |  91.51
                 |   0.00 |   6.29 |   0.00 |  89.01 |   4.70 |
                 |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
        ---------+--------+--------+--------+--------+--------+
               1 |   1412 |      0 |   1739 |      0 |      0 |   3151
                 |   3.81 |   0.00 |   4.69 |   0.00 |   0.00 |   8.49
                 |  44.81 |   0.00 |  55.19 |   0.00 |   0.00 |
                 | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
        ---------+--------+--------+--------+--------+--------+
        Total        1412     2137     1739    30223     1595    37106
                     3.81     5.76     4.69    81.45     4.30   100.00

                           Frequency Missing = 1389

 Table of flag_xscript_01gy_v_Mock_3h_dd by sign_01gy_v_Mock_3h_fdr05

    flag_xscript_01gy_v_Mock_3h_dd     sign_01gy_v_Mock_3h_fdr05

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |0.1gy   |D       |Mock    |N       |U       |  Total
    ---------+--------+--------+--------+--------+--------+
           0 |      0 |     30 |      0 |  34215 |     37 |  34282
             |   0.00 |   0.08 |   0.00 |  91.99 |   0.10 |  92.17
             |   0.00 |   0.09 |   0.00 |  99.80 |   0.11 |
             |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
    ---------+--------+--------+--------+--------+--------+
           1 |   1434 |      0 |   1478 |      0 |      0 |   2912
             |   3.86 |   0.00 |   3.97 |   0.00 |   0.00 |   7.83
             |  49.24 |   0.00 |  50.76 |   0.00 |   0.00 |
             | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
    ---------+--------+--------+--------+--------+--------+
    Total        1434       30     1478    34215       37    37194
                 3.86     0.08     3.97    91.99     0.10   100.00

                       Frequency Missing = 1301
.

Table of flag_xscript_01gy_v_Mock_24h_dd by sign_01gy_v_Mock_24h_fdr05

    flag_xscript_01gy_v_Mock_24h_dd     sign_01gy_v_Mock_24h_fdr05

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |0.1gy   |D       |Mock    |N       |U       |  Total
    ---------+--------+--------+--------+--------+--------+
           0 |      0 |    174 |      0 |  33879 |     64 |  34117
             |   0.00 |   0.47 |   0.00 |  91.17 |   0.17 |  91.81
             |   0.00 |   0.51 |   0.00 |  99.30 |   0.19 |
             |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
    ---------+--------+--------+--------+--------+--------+
           1 |   1561 |      0 |   1484 |      0 |      0 |   3045
             |   4.20 |   0.00 |   3.99 |   0.00 |   0.00 |   8.19
             |  51.26 |   0.00 |  48.74 |   0.00 |   0.00 |
             | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
    ---------+--------+--------+--------+--------+--------+
    Total        1561      174     1484    33879       64    37162
                 4.20     0.47     3.99    91.17     0.17   100.00
Table of flag_xscript_01gy_v_Mock_72h_dd by sign_01gy_v_Mock_72h_fdr05

    flag_xscript_01gy_v_Mock_72h_dd     sign_01gy_v_Mock_72h_fdr05

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |0.1gy   |D       |Mock    |N       |U       |  Total
    ---------+--------+--------+--------+--------+--------+
           0 |      0 |    856 |      0 |  33210 |    286 |  34352
             |   0.00 |   2.37 |   0.00 |  91.92 |   0.79 |  95.08
             |   0.00 |   2.49 |   0.00 |  96.68 |   0.83 |
             |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
    ---------+--------+--------+--------+--------+--------+
           1 |    313 |      0 |   1466 |      0 |      0 |   1779
             |   0.87 |   0.00 |   4.06 |   0.00 |   0.00 |   4.92
             |  17.59 |   0.00 |  82.41 |   0.00 |   0.00 |
             | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
    ---------+--------+--------+--------+--------+--------+
    Total         313      856     1466    33210      286    36131
                 0.87     2.37     4.06    91.92     0.79   100.00

                       Frequency Missing = 2364
                       Frequency Missing = 1333


  Table of flag_xscript_1gy_v_Mock_1h_dd by sign_1gy_v_Mock_1h_fdr05

    flag_xscript_1gy_v_Mock_1h_dd     sign_1gy_v_Mock_1h_fdr05

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |1gy     |D       |Mock    |N       |U       |  Total
    ---------+--------+--------+--------+--------+--------+
           0 |      0 |   1267 |      0 |  31749 |   1194 |  34210
             |   0.00 |   3.41 |   0.00 |  85.34 |   3.21 |  91.95
             |   0.00 |   3.70 |   0.00 |  92.81 |   3.49 |
             |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
    ---------+--------+--------+--------+--------+--------+
           1 |   1509 |      0 |   1484 |      0 |      0 |   2993
             |   4.06 |   0.00 |   3.99 |   0.00 |   0.00 |   8.05
             |  50.42 |   0.00 |  49.58 |   0.00 |   0.00 |
             | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
    ---------+--------+--------+--------+--------+--------+
    Total        1509     1267     1484    31749     1194    37203
                 4.06     3.41     3.99    85.34     3.21   100.00

                       Frequency Missing = 1292


Table of flag_xscript_1gy_v_Mock_3h_dd by sign_1gy_v_Mock_3h_fdr05

  flag_xscript_1gy_v_Mock_3h_dd     sign_1gy_v_Mock_3h_fdr05

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |1gy     |D       |Mock    |N       |U       |  Total
  ---------+--------+--------+--------+--------+--------+
         0 |      0 |    512 |      0 |  33677 |    193 |  34382
           |   0.00 |   1.37 |   0.00 |  90.36 |   0.52 |  92.26
           |   0.00 |   1.49 |   0.00 |  97.95 |   0.56 |
           |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
  ---------+--------+--------+--------+--------+--------+
         1 |   1508 |      0 |   1378 |      0 |      0 |   2886
           |   4.05 |   0.00 |   3.70 |   0.00 |   0.00 |   7.74
           |  52.25 |   0.00 |  47.75 |   0.00 |   0.00 |
           | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
  ---------+--------+--------+--------+--------+--------+
  Total        1508      512     1378    33677      193    37268
               4.05     1.37     3.70    90.36     0.52   100.00

                     Frequency Missing = 1227



   Table of flag_xscript_1gy_v_Mock_24h_dd by sign_1gy_v_Mock_24h_fdr05

      flag_xscript_1gy_v_Mock_24h_dd     sign_1gy_v_Mock_24h_fdr05

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |1gy     |D       |Mock    |N       |U       |  Total
      ---------+--------+--------+--------+--------+--------+
             0 |      0 |      2 |      0 |  34296 |      5 |  34303
               |   0.00 |   0.01 |   0.00 |  92.26 |   0.01 |  92.28
               |   0.00 |   0.01 |   0.00 |  99.98 |   0.01 |
               |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
      ---------+--------+--------+--------+--------+--------+
             1 |   1571 |      0 |   1298 |      0 |      0 |   2869
               |   4.23 |   0.00 |   3.49 |   0.00 |   0.00 |   7.72
               |  54.76 |   0.00 |  45.24 |   0.00 |   0.00 |
               | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
      ---------+--------+--------+--------+--------+--------+
      Total        1571        2     1298    34296        5    37172
                   4.23     0.01     3.49    92.26     0.01   100.00

                         Frequency Missing = 1323

                              The SAS System          11:17 Thursday, April 5, 20



     Table of flag_xscript_1gy_v_Mock_72h_dd by sign_1gy_v_Mock_72h_fdr05

        flag_xscript_1gy_v_Mock_72h_dd     sign_1gy_v_Mock_72h_fdr05

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |1gy     |D       |Mock    |N       |U       |  Total
        ---------+--------+--------+--------+--------+--------+
               0 |      0 |      6 |      0 |  34592 |      3 |  34601
                 |   0.00 |   0.02 |   0.00 |  92.26 |   0.01 |  92.28
                 |   0.00 |   0.02 |   0.00 |  99.97 |   0.01 |
                 |   0.00 | 100.00 |   0.00 | 100.00 | 100.00 |
        ---------+--------+--------+--------+--------+--------+
               1 |   1677 |      0 |   1217 |      0 |      0 |   2894
                 |   4.47 |   0.00 |   3.25 |   0.00 |   0.00 |   7.72
                 |  57.95 |   0.00 |  42.05 |   0.00 |   0.00 |
                 | 100.00 |   0.00 | 100.00 |   0.00 |   0.00 |
        ---------+--------+--------+--------+--------+--------+
        Total        1677        6     1217    34592        3    37495
                     4.47     0.02     3.25    92.26     0.01   100.00

                           Frequency Missing = 1000




*/

















