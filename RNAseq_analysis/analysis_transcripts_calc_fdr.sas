ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* For the contrasts Mingqi wants to look at further,
   calc FDR and flag FDR 5%, 10%, 20%, and update corresponding up/down indicators */

%macro calcFDR(onflag1,onflag2,labelname,outname,fdrname);

data on_xs;
   set rs.arab_flag_transcript_on_gt0;
   where &onflag1.=1 and &onflag2.=1 ;
   keep transcript_id;
run;

data contrast_p;
  set rs.arab_xs_cntrs_constr_non;
  format ProbF best32. ;
  where label=&labelname.;
  keep transcript_id probf;
run;  

proc sort data=contrast_p;
  by transcript_id;
proc sort data=on_xs;
  by transcript_id;
run;

data contrast_p2;
  merge on_xs (in=in1) contrast_p (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc multtest inpvalues(ProbF)=contrast_p2 fdr
 out=fdr noprint;
run;
quit;

data &outname.;
  set fdr;
  if fdr_p = . then flag_&fdrname._fdr05=.;
     else if fdr_p <0.05 then flag_&fdrname._fdr05=1;
     else flag_&fdrname._fdr05=0;

  if fdr_p = . then flag_&fdrname._fdr10=.;
     else if fdr_p <0.10 then flag_&fdrname._fdr10=1;
     else flag_&fdrname._fdr10=0;

  if fdr_p = . then flag_&fdrname._fdr20=.;
     else if fdr_p <0.20 then flag_&fdrname._fdr20=1;
     else flag_&fdrname._fdr20=0;

  keep transcript_id ProbF fdr_p flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
  rename fdr_p=fdr_&fdrname. ProbF=p_&fdrname.;
run;

proc freq data=&outname.;
   tables flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
run;

%mend;

%calcFDR(flag_xscript_on_01gy_1hr_tpm0,flag_xscript_on_Mock_1hr_tpm0,"0.1gy-Mock: 1h",contrast1,01gy_v_Mock_1h);
%calcFDR(flag_xscript_on_01gy_3hr_tpm0,flag_xscript_on_Mock_3hr_tpm0,"0.1gy-Mock: 3h",contrast2,01gy_v_Mock_3h);
%calcFDR(flag_xscript_on_01gy_24hr_tpm0,flag_xscript_on_Mock_24hr_tpm0,"0.1gy-Mock: 24h",contrast3,01gy_v_Mock_24h);
%calcFDR(flag_xscript_on_01gy_72hr_tpm0,flag_xscript_on_Mock_72hr_tpm0,"0.1gy-Mock: 72h",contrast4,01gy_v_Mock_72h);

%calcFDR(flag_xscript_on_1gy_1hr_tpm0,flag_xscript_on_Mock_1hr_tpm0,"1gy-Mock: 1h",contrast5,1gy_v_Mock_1h);
%calcFDR(flag_xscript_on_1gy_3hr_tpm0,flag_xscript_on_Mock_3hr_tpm0,"1gy-Mock: 3h",contrast6,1gy_v_Mock_3h);
%calcFDR(flag_xscript_on_1gy_24hr_tpm0,flag_xscript_on_Mock_24hr_tpm0,"1gy-Mock: 24h",contrast7,1gy_v_Mock_24h);
%calcFDR(flag_xscript_on_1gy_72hr_tpm0,flag_xscript_on_Mock_72hr_tpm0,"1gy-Mock: 72h",contrast8,1gy_v_Mock_72h);



%macro calcFDRMain(labelname,outname,fdrname);

data main_p;
  set rs.arab_anova_xs_main_trt_tm_non;
  format ProbF best32. ;
  where effect=&labelname.;
  keep transcript_id probf;
run;  

proc multtest inpvalues(ProbF)=main_p fdr
 out=fdr noprint;
run;
quit;

data &outname.;
  set fdr;
  if fdr_p = . then flag_&fdrname._fdr05=.;
     else if fdr_p <0.05 then flag_&fdrname._fdr05=1;
     else flag_&fdrname._fdr05=0;

  if fdr_p = . then flag_&fdrname._fdr10=.;
     else if fdr_p <0.10 then flag_&fdrname._fdr10=1;
     else flag_&fdrname._fdr10=0;

  if fdr_p = . then flag_&fdrname._fdr20=.;
     else if fdr_p <0.20 then flag_&fdrname._fdr20=1;
     else flag_&fdrname._fdr20=0;

  keep transcript_id probf fdr_p flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
  rename fdr_p=fdr_&fdrname. ProbF=p_&fdrname.;
run;

proc freq data=&outname.;
   tables flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
run;

%mend;

%calcFDRMain("treatment",main1,treatment);
%calcFDRMain("time",main2,time);
%calcFDRMain("time*treatment",main3,trt_by_time);




/* Merge FDRs and make permenant */

proc sort data=contrast1;
  by transcript_id;
proc sort data=contrast2;
  by transcript_id;
proc sort data=contrast3;
  by transcript_id;
proc sort data=contrast4;
  by transcript_id;
proc sort data=contrast5;
  by transcript_id;
proc sort data=contrast6;
  by transcript_id;
proc sort data=contrast7;
  by transcript_id;
proc sort data=contrast8;
  by transcript_id;
proc sort data=main1;
  by transcript_id;
proc sort data=main2;
  by transcript_id;
proc sort data=main3;
  by transcript_id;
run;

data fdr_by_xscript;
  merge main1 main2 main3 contrast1 contrast2 contrast3 contrast4 contrast5 contrast6 contrast7 contrast8;
  by transcript_id;
run;

data rs.fdr_by_transcript;
  set fdr_by_xscript;
run;




/*
    flag_01gy_v_                             Cumulative    Cumulative
   Mock_1h_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       30223       89.01         30223        89.01
               1        3732       10.99         33955       100.00


    flag_01gy_v_                             Cumulative    Cumulative
   Mock_1h_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       29057       85.58         29057        85.58
               1        4898       14.42         33955       100.00


    flag_01gy_v_                             Cumulative    Cumulative
   Mock_1h_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       26960       79.40         26960        79.40
               1        6995       20.60         33955       100.00

      flag_01gy_v_                             Cumulative    Cumulative
     Mock_3h_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       34215       99.80         34215        99.80
                 1          67        0.20         34282       100.00


      flag_01gy_v_                             Cumulative    Cumulative
     Mock_3h_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       34120       99.53         34120        99.53
                 1         162        0.47         34282       100.00


      flag_01gy_v_                             Cumulative    Cumulative
     Mock_3h_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       33659       98.18         33659        98.18
                 1         623        1.82         34282       100.00

     flag_01gy_v_                             Cumulative    Cumulative
   Mock_24h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33879       99.30         33879        99.30
                1         238        0.70         34117       100.00


     flag_01gy_v_                             Cumulative    Cumulative
   Mock_24h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33653       98.64         33653        98.64
                1         464        1.36         34117       100.00


     flag_01gy_v_                             Cumulative    Cumulative
   Mock_24h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33037       96.83         33037        96.83
                1        1080        3.17         34117       100.00

     flag_01gy_v_                             Cumulative    Cumulative
   Mock_72h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33210       96.68         33210        96.68
                1        1142        3.32         34352       100.00


     flag_01gy_v_                             Cumulative    Cumulative
   Mock_72h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       32566       94.80         32566        94.80
                1        1786        5.20         34352       100.00


     flag_01gy_v_                             Cumulative    Cumulative
   Mock_72h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       31451       91.56         31451        91.56
                1        2901        8.44         34352       100.00

 flag_1gy_v_Mock_                             Cumulative    Cumulative
         1h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       31749       92.81         31749        92.81
                1        2461        7.19         34210       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
         1h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       30822       90.10         30822        90.10
                1        3388        9.90         34210       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
         1h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       29355       85.81         29355        85.81
                1        4855       14.19         34210       100.00

 flag_1gy_v_Mock_                             Cumulative    Cumulative
         3h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33677       97.95         33677        97.95
                1         705        2.05         34382       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
         3h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33042       96.10         33042        96.10
                1        1340        3.90         34382       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
         3h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       31686       92.16         31686        92.16
                1        2696        7.84         34382       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        24h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       34296       99.98         34296        99.98
                1           7        0.02         34303       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        24h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       34278       99.93         34278        99.93
                1          25        0.07         34303       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        24h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       34112       99.44         34112        99.44
                1         191        0.56         34303       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        72h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       34592       99.97         34592        99.97
                1           9        0.03         34601       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        72h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       34527       99.79         34527        99.79
                1          74        0.21         34601       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        72h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33999       98.26         33999        98.26
                1         602        1.74         34601       100.00


 flag_treatment_                             Cumulative    Cumulative
           fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       35106       91.20         35106        91.20
               1        3389        8.80         38495       100.00


 flag_treatment_                             Cumulative    Cumulative
           fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       33731       87.62         33731        87.62
               1        4764       12.38         38495       100.00


 flag_treatment_                             Cumulative    Cumulative
           fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       31223       81.11         31223        81.11
               1        7272       18.89         38495       100.00


                                             Cumulative    Cumulative
 flag_time_fdr05    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0       16552       43.00         16552        43.00
               1       21943       57.00         38495       100.00


                                             Cumulative    Cumulative
 flag_time_fdr10    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0       14405       37.42         14405        37.42
               1       24090       62.58         38495       100.00


                                             Cumulative    Cumulative
 flag_time_fdr20    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0       11580       30.08         11580        30.08
               1       26915       69.92         38495       100.00

     flag_trt_by_                             Cumulative    Cumulative
       time_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       36330       94.38         36330        94.38
                1        2165        5.62         38495       100.00


     flag_trt_by_                             Cumulative    Cumulative
       time_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       35250       91.57         35250        91.57
                1        3245        8.43         38495       100.00


     flag_trt_by_                             Cumulative    Cumulative
       time_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       33035       85.82         33035        85.82
                1        5460       14.18         38495       100.00

*/


