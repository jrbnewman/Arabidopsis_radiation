ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* For the contrasts Mingqi wants to look at further,
   calc FDR and flag FDR 5%, 10%, 20%, and update corresponding up/down indicators */

%macro calcFDR(onflag1,onflag2,labelname,outname,fdrname);

data on_genes;
   set rs.arab_flag_gene_on_gt0;
   where &onflag1.=1 and &onflag2.=1 ;
   keep gene_id;
run;

data contrast_p;
  set rs.arab_gene_cntrs_constr;
  format ProbF best32. ;
  where label=&labelname.;
  keep gene_id probf;
run;  

proc sort data=contrast_p;
  by gene_id;
proc sort data=on_genes;
  by gene_id;
run;

data contrast_p2;
  merge on_genes (in=in1) contrast_p (in=in2);
  by gene_id;
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

  keep gene_id ProbF fdr_p flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
  rename fdr_p=fdr_&fdrname. ProbF=p_&fdrname.;
run;

proc freq data=&outname.;
   tables flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
run;

%mend;

%calcFDR(flag_gene_on_01gy_1hr_apn0,flag_gene_on_Mock_1hr_apn0,"0.1gy-Mock: 1h",contrast1,01gy_v_Mock_1h);
%calcFDR(flag_gene_on_01gy_3hr_apn0,flag_gene_on_Mock_3hr_apn0,"0.1gy-Mock: 3h",contrast2,01gy_v_Mock_3h);
%calcFDR(flag_gene_on_01gy_24hr_apn0,flag_gene_on_Mock_24hr_apn0,"0.1gy-Mock: 24h",contrast3,01gy_v_Mock_24h);
%calcFDR(flag_gene_on_01gy_72hr_apn0,flag_gene_on_Mock_72hr_apn0,"0.1gy-Mock: 72h",contrast4,01gy_v_Mock_72h);

%calcFDR(flag_gene_on_1gy_1hr_apn0,flag_gene_on_Mock_1hr_apn0,"1gy-Mock: 1h",contrast5,1gy_v_Mock_1h);
%calcFDR(flag_gene_on_1gy_3hr_apn0,flag_gene_on_Mock_3hr_apn0,"1gy-Mock: 3h",contrast6,1gy_v_Mock_3h);
%calcFDR(flag_gene_on_1gy_24hr_apn0,flag_gene_on_Mock_24hr_apn0,"1gy-Mock: 24h",contrast7,1gy_v_Mock_24h);
%calcFDR(flag_gene_on_1gy_72hr_apn0,flag_gene_on_Mock_72hr_apn0,"1gy-Mock: 72h",contrast8,1gy_v_Mock_72h);



%macro calcFDRMain(labelname,outname,fdrname);

data main_p;
  set rs.arab_anova_main_trt_tm_genes;
  format ProbF best32. ;
  where effect=&labelname.;
  keep gene_id probf;
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

  keep gene_id probf fdr_p flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
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
  by gene_id;
proc sort data=contrast2;
  by gene_id;
proc sort data=contrast3;
  by gene_id;
proc sort data=contrast4;
  by gene_id;
proc sort data=contrast5;
  by gene_id;
proc sort data=contrast6;
  by gene_id;
proc sort data=contrast7;
  by gene_id;
proc sort data=contrast8;
  by gene_id;
proc sort data=main1;
  by gene_id;
proc sort data=main2;
  by gene_id;
proc sort data=main3;
  by gene_id;
run;

data fdr_by_gene;
  merge main1 main2 main3 contrast1 contrast2 contrast3 contrast4 contrast5 contrast6 contrast7 contrast8;
  by gene_id;
run;

data rs.fdr_by_gene_v2;
  set fdr_by_gene;
run;




/*
    flag_01gy_v_                             Cumulative    Cumulative
   Mock_1h_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       16364       79.44         16364        79.44
               1        4235       20.56         20599       100.00


    flag_01gy_v_                             Cumulative    Cumulative
   Mock_1h_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       15139       73.49         15139        73.49
               1        5460       26.51         20599       100.00


    flag_01gy_v_                             Cumulative    Cumulative
   Mock_1h_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       13291       64.52         13291        64.52
               1        7308       35.48         20599       100.00

     flag_01gy_v_                             Cumulative    Cumulative
    Mock_3h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       20555       99.38         20555        99.38
                1         129        0.62         20684       100.00


     flag_01gy_v_                             Cumulative    Cumulative
    Mock_3h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       20362       98.44         20362        98.44
                1         322        1.56         20684       100.00


     flag_01gy_v_                             Cumulative    Cumulative
    Mock_3h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       19676       95.13         19676        95.13
                1        1008        4.87         20684       100.00

    flag_01gy_v_                             Cumulative    Cumulative
  Mock_24h_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       20399       98.50         20399        98.50
               1         310        1.50         20709       100.00


    flag_01gy_v_                             Cumulative    Cumulative
  Mock_24h_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       20138       97.24         20138        97.24
               1         571        2.76         20709       100.00


    flag_01gy_v_                             Cumulative    Cumulative
  Mock_24h_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       19300       93.20         19300        93.20
               1        1409        6.80         20709       100.00

    flag_01gy_v_                             Cumulative    Cumulative
  Mock_72h_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       19261       93.10         19261        93.10
               1        1427        6.90         20688       100.00


    flag_01gy_v_                             Cumulative    Cumulative
  Mock_72h_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       18448       89.17         18448        89.17
               1        2240       10.83         20688       100.00


    flag_01gy_v_                             Cumulative    Cumulative
  Mock_72h_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       17120       82.75         17120        82.75
               1        3568       17.25         20688       100.00

flag_1gy_v_Mock_                             Cumulative    Cumulative
        1h_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       17845       86.31         17845        86.31
               1        2831       13.69         20676       100.00


flag_1gy_v_Mock_                             Cumulative    Cumulative
        1h_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       16931       81.89         16931        81.89
               1        3745       18.11         20676       100.00


flag_1gy_v_Mock_                             Cumulative    Cumulative
        1h_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       15469       74.82         15469        74.82
               1        5207       25.18         20676       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
          3h_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       19971       96.26         19971        96.26
                 1         776        3.74         20747       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
          3h_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       19150       92.30         19150        92.30
                 1        1597        7.70         20747       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
          3h_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       17606       84.86         17606        84.86
                 1        3141       15.14         20747       100.00

 flag_1gy_v_Mock_                             Cumulative    Cumulative
        24h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       20694       99.81         20694        99.81
                1          39        0.19         20733       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        24h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       20552       99.13         20552        99.13
                1         181        0.87         20733       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
        24h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       19963       96.29         19963        96.29
                1         770        3.71         20733       100.00

  flag_1gy_v_Mock_                             Cumulative    Cumulative
         72h_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       20703       99.74         20703        99.74
                 1          53        0.26         20756       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
         72h_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       20480       98.67         20480        98.67
                 1         276        1.33         20756       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
         72h_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       19438       93.65         19438        93.65
                 1        1318        6.35         20756       100.00



   flag_treatment_                             Cumulative    Cumulative
             fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       17627       81.09         17627        81.09
                 1        4111       18.91         21738       100.00


   flag_treatment_                             Cumulative    Cumulative
             fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       16062       73.89         16062        73.89
                 1        5676       26.11         21738       100.00


   flag_treatment_                             Cumulative    Cumulative
             fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       13748       63.24         13748        63.24
                 1        7990       36.76         21738       100.00

                                            Cumulative    Cumulative
flag_time_fdr05    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0        4741       21.81          4741        21.81
              1       16997       78.19         21738       100.00


                                            Cumulative    Cumulative
flag_time_fdr10    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0        3893       17.91          3893        17.91
              1       17845       82.09         21738       100.00


                                            Cumulative    Cumulative
flag_time_fdr20    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0        3001       13.81          3001        13.81
              1       18737       86.19         21738       100.00


      flag_trt_by_                             Cumulative    Cumulative
        time_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       18781       86.40         18781        86.40
                 1        2957       13.60         21738       100.00


      flag_trt_by_                             Cumulative    Cumulative
        time_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       17429       80.18         17429        80.18
                 1        4309       19.82         21738       100.00


      flag_trt_by_                             Cumulative    Cumulative
        time_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       15024       69.11         15024        69.11
                 1        6714       30.89         21738       100.00


*/


