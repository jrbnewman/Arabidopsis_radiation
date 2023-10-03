ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* For the contrasts Mingqi wants to look at further,
   calc FDR and flag FDR 5%, 10%, 20%, and update corresponding up/down indicators */

%macro calcFDR(labelname,outname,fdrname);

data contrast_p;
  set rs.arab_gene_cntrs_constr;
  format ProbF best32. ;
  where label=&labelname.;
  keep gene_id probf;
run;  

proc multtest inpvalues(ProbF)=contrast_p fdr
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

  keep gene_id fdr_p flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
  rename fdr_p=fdr_&fdrname.;
run;

proc freq data=&outname.;
   tables flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
run;

%mend;

%calcFDR("0.1gy-Mock: 1h",contrast1,01gy_v_Mock_1h);
%calcFDR("0.1gy-Mock: 3h",contrast2,01gy_v_Mock_3h);
%calcFDR("0.1gy-Mock: 24h",contrast3,01gy_v_Mock_24h);
%calcFDR("0.1gy-Mock: 72h",contrast4,01gy_v_Mock_72h);

%calcFDR("1gy-Mock: 1h",contrast5,1gy_v_Mock_1h);
%calcFDR("1gy-Mock: 3h",contrast6,1gy_v_Mock_3h);
%calcFDR("1gy-Mock: 24h",contrast7,1gy_v_Mock_24h);
%calcFDR("1gy-Mock: 72h",contrast8,1gy_v_Mock_72h);



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

  keep gene_id fdr_p flag_&fdrname._fdr05 flag_&fdrname._fdr10  flag_&fdrname._fdr20;
  rename fdr_p=fdr_&fdrname.;
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

data rs.fdr_by_gene;
  set fdr_by_gene;
run;




/*
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



     flag_01gy_v_                             Cumulative    Cumulative
    Mock_1h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       17523       80.61         17523        80.61
                1        4215       19.39         21738       100.00


     flag_01gy_v_                             Cumulative    Cumulative
    Mock_1h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       16299       74.98         16299        74.98
                1        5439       25.02         21738       100.00


     flag_01gy_v_                             Cumulative    Cumulative
    Mock_1h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       14410       66.29         14410        66.29
                1        7328       33.71         21738       100.00

       flag_01gy_v_                             Cumulative    Cumulative
     Mock_3h_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       21613       99.42         21613        99.42
                 1         125        0.58         21738       100.00


      flag_01gy_v_                             Cumulative    Cumulative
     Mock_3h_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       21422       98.55         21422        98.55
                 1         316        1.45         21738       100.00


      flag_01gy_v_                             Cumulative    Cumulative
     Mock_3h_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       20795       95.66         20795        95.66
                 1         943        4.34         21738       100.00
      flag_01gy_v_                             Cumulative    Cumulative
    Mock_24h_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       21440       98.63         21440        98.63
                 1         298        1.37         21738       100.00


      flag_01gy_v_                             Cumulative    Cumulative
    Mock_24h_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       21157       97.33         21157        97.33
                 1         581        2.67         21738       100.00


      flag_01gy_v_                             Cumulative    Cumulative
    Mock_24h_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       20336       93.55         20336        93.55
                 1        1402        6.45         21738       100.00

     flag_01gy_v_                             Cumulative    Cumulative
   Mock_72h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       20351       93.62         20351        93.62
                1        1387        6.38         21738       100.00


     flag_01gy_v_                             Cumulative    Cumulative
   Mock_72h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       19534       89.86         19534        89.86
                1        2204       10.14         21738       100.00


     flag_01gy_v_                             Cumulative    Cumulative
   Mock_72h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18219       83.81         18219        83.81
                1        3519       16.19         21738       100.00



 flag_1gy_v_Mock_                             Cumulative    Cumulative
         1h_fdr05    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18943       87.14         18943        87.14
                1        2795       12.86         21738       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
         1h_fdr10    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18019       82.89         18019        82.89
                1        3719       17.11         21738       100.00


 flag_1gy_v_Mock_                             Cumulative    Cumulative
         1h_fdr20    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       16548       76.12         16548        76.12
                1        5190       23.88         21738       100.00

   flag_1gy_v_Mock_                             Cumulative    Cumulative
           3h_fdr05    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0       20995       96.58         20995        96.58
                  1         743        3.42         21738       100.00


   flag_1gy_v_Mock_                             Cumulative    Cumulative
           3h_fdr10    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0       20194       92.90         20194        92.90
                  1        1544        7.10         21738       100.00


   flag_1gy_v_Mock_                             Cumulative    Cumulative
           3h_fdr20    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0       18649       85.79         18649        85.79
                  1        3089       14.21         21738       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
         24h_fdr05    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       21698       99.82         21698        99.82
                 1          40        0.18         21738       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
         24h_fdr10    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       21551       99.14         21551        99.14
                 1         187        0.86         21738       100.00


  flag_1gy_v_Mock_                             Cumulative    Cumulative
         24h_fdr20    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       20990       96.56         20990        96.56
                 1         748        3.44         21738       100.00


flag_1gy_v_Mock_                             Cumulative    Cumulative
       72h_fdr05    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       21681       99.74         21681        99.74
               1          57        0.26         21738       100.00


flag_1gy_v_Mock_                             Cumulative    Cumulative
       72h_fdr10    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       21443       98.64         21443        98.64
               1         295        1.36         21738       100.00


flag_1gy_v_Mock_                             Cumulative    Cumulative
       72h_fdr20    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       20362       93.67         20362        93.67
               1        1376        6.33         21738       100.00



*/


