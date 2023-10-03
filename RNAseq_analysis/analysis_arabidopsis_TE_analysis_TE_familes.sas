ods listing; ods html close;
libname arabRNA "!HOME/concannon/DTRA/arabidopsis/sas_data";

proc datasets lib=work kill noprint;
run;
quit;

/* TE local counts only!!! */

/* (1) FLAG ON/OFF */

/* Flag genes on/off, CPM>0 in at least 50% of samples per treatment */

%macro flagOnOff(countType);

data counts_w_key;
  set arabRNA.cpm_by_te_&countType.;
  where flag_TE=1;
  if cpm > 0 then flag_cpm_gt0=1;
  else flag_cpm_gt0=0;

  if cpm ge 5 then flag_cpm_ge5=1;
  else flag_cpm_ge5=0;
run;


/* Flag on per treatment */

proc sort data=counts_w_key;
   by treatment gene_id;
proc means data=counts_w_key noprint;
   by treatment gene_id;
   var flag_cpm_gt0 flag_cpm_ge5;
   output out=mean_on_by_trt mean=;
run;

data flag_iso_on_by_trt;
  set mean_on_by_trt;

  /* For CPM > 0 */
  if flag_cpm_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_cpm_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  /* For CPM ge 5 */
  if flag_cpm_ge5 ge 0.5 then flag_gene_on_ge5=1;
  else if flag_cpm_ge5 = 0 then flag_gene_on_ge5=0;
  else flag_gene_on_ge5=.;

  keep gene_id treatment flag_gene_on_gt0 flag_gene_on_ge5;
run;

proc sort data=flag_iso_on_by_trt;
   by gene_id treatment;
proc transpose data=flag_iso_on_by_trt out=flag_iso_trt_sbys_cpm0;
   by gene_id;
   id treatment;
   var flag_gene_on_gt0;
run;

proc transpose data=flag_iso_on_by_trt out=flag_iso_trt_sbys_cpm5;
   by gene_id;
   id treatment;
   var flag_gene_on_ge5;
run;

data treat_cpm0;
   set flag_iso_trt_sbys_cpm0;
   drop _NAME_;
   rename _0_1gy=flag_gene_on_01gy_cpm0
          _1gy=flag_gene_on_1gy_cpm0
          Mock=flag_gene_on_mock_cpm0;
run;

data treat_cpm5;
   set flag_iso_trt_sbys_cpm5;
   drop _NAME_;
   rename _0_1gy=flag_gene_on_01gy_cpm5
          _1gy=flag_gene_on_1gy_cpm5
          Mock=flag_gene_on_mock_cpm5;
run;

/* Flag on per time */

proc sort data=counts_w_key;
   by time gene_id;
proc means data=counts_w_key noprint;
   by time gene_id;
   var flag_cpm_gt0 flag_cpm_ge5;
   output out=mean_on_by_time mean=;
run;

data flag_iso_on_by_time;
  set mean_on_by_time;

  /* For APN > 0 */
  if flag_cpm_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_cpm_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  /* For APN ge 5 */
  if flag_cpm_ge5 ge 0.5 then flag_gene_on_ge5=1;
  else if flag_cpm_ge5 = 0 then flag_gene_on_ge5=0;
  else flag_gene_on_ge5=.;

  keep gene_id time flag_gene_on_gt0 flag_gene_on_ge5;
run;

proc sort data=flag_iso_on_by_time;
   by gene_id time;
proc transpose data=flag_iso_on_by_time out=flag_iso_time_sbys_cpm0;
   by gene_id;
   id time;
   var flag_gene_on_gt0;
run;

proc transpose data=flag_iso_on_by_time out=flag_iso_time_sbys_cpm5;
   by gene_id;
   id time;
   var flag_gene_on_ge5;
run;

data time_cpm0;
   set flag_iso_time_sbys_cpm0;
   drop _NAME_;
   rename _1=flag_gene_on_1hr_cpm0
          _3=flag_gene_on_3hr_cpm0
          _24=flag_gene_on_24hr_cpm0
          _72=flag_gene_on_72hr_cpm0;
run;

data time_cpm5;
   set flag_iso_time_sbys_cpm5;
   drop _NAME_;
   rename _1=flag_gene_on_1hr_cpm5
          _3=flag_gene_on_3hr_cpm5
          _24=flag_gene_on_24hr_cpm5
          _72=flag_gene_on_72hr_cpm5;
run;

/* Flag on per treatment-by-time */

proc sort data=counts_w_key;
   by treatment time gene_id;
proc means data=counts_w_key noprint;
   by treatment time gene_id;
   var flag_cpm_gt0 flag_cpm_ge5;
   output out=mean_on_by_trttime mean=;
run;

data flag_iso_on_by_trttime;
  length condition $10.;
  set mean_on_by_trttime;
  condition=catx("_",treatment,time);

  /* For TPM > 0 */
  if flag_cpm_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_cpm_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  /* For TPM ge 5 */
  if flag_cpm_ge5 ge 0.5 then flag_gene_on_ge5=1;
  else if flag_cpm_ge5 = 0 then flag_gene_on_ge5=0;
  else flag_gene_on_ge5=.;

  keep gene_id condition flag_gene_on_gt0 flag_gene_on_ge5;
run;

proc sort data=flag_iso_on_by_trttime;
   by gene_id condition;
proc transpose data=flag_iso_on_by_trttime out=flag_iso_trttime_sbys_cpm0;
   by gene_id;
   id condition;
   var flag_gene_on_gt0;
run;

proc transpose data=flag_iso_on_by_trttime out=flag_iso_trttime_sbys_cpm5;
   by gene_id;
   id condition;
   var flag_gene_on_ge5;
run;

data trttime_cpm0;
   set flag_iso_trttime_sbys_cpm0;
   drop _NAME_;
   rename _0_1gy_1=flag_gene_on_01gy_1hr_cpm0
          _0_1gy_3=flag_gene_on_01gy_3hr_cpm0
          _0_1gy_24=flag_gene_on_01gy_24hr_cpm0
          _0_1gy_72=flag_gene_on_01gy_72hr_cpm0
          _1gy_1=flag_gene_on_1gy_1hr_cpm0
          _1gy_3=flag_gene_on_1gy_3hr_cpm0
          _1gy_24=flag_gene_on_1gy_24hr_cpm0
          _1gy_72=flag_gene_on_1gy_72hr_cpm0
          Mock_1=flag_gene_on_mock_1hr_cpm0
          Mock_3=flag_gene_on_mock_3hr_cpm0
          Mock_24=flag_gene_on_mock_24hr_cpm0
          Mock_72=flag_gene_on_mock_72hr_cpm0;
run;

data trttime_cpm5;
   set flag_iso_trttime_sbys_cpm5;
   drop _NAME_;
   rename _0_1gy_1=flag_gene_on_01gy_1hr_cpm5
          _0_1gy_3=flag_gene_on_01gy_3hr_cpm5
          _0_1gy_24=flag_gene_on_01gy_24hr_cpm5
          _0_1gy_72=flag_gene_on_01gy_72hr_cpm5
          _1gy_1=flag_gene_on_1gy_1hr_cpm5
          _1gy_3=flag_gene_on_1gy_3hr_cpm5
          _1gy_24=flag_gene_on_1gy_24hr_cpm5
          _1gy_72=flag_gene_on_1gy_72hr_cpm5
          Mock_1=flag_gene_on_mock_1hr_cpm5
          Mock_3=flag_gene_on_mock_3hr_cpm5
          Mock_24=flag_gene_on_mock_24hr_cpm5
          Mock_72=flag_gene_on_mock_72hr_cpm5;
run;

/* Merge and make permenant */

proc sort data=treat_cpm0;
   by gene_id;
proc sort data=time_cpm0;
   by gene_id;
proc sort data=trttime_cpm0;
   by gene_id;
proc sort data=treat_cpm5;
   by gene_id;
proc sort data=time_cpm5;
   by gene_id;
proc sort data=trttime_cpm5;
   by gene_id;
run;




data arabRNA.arab_flag_TE_&countType._on_cpm_gt0_2;
  merge treat_cpm0 (in=in1) time_cpm0 (in=in2) trttime_cpm0 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

data arabRNA.arab_flag_TE_&countType._on_cpm_ge5_2;
  merge treat_cpm5 (in=in1) time_cpm5 (in=in2) trttime_cpm5 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

/* Count CPM > 0 */

ods listing;
proc freq data=arabRNA.arab_flag_TE_&countType._on_cpm_gt0_2;
   tables flag_gene_on_01gy_cpm0 flag_gene_on_1gy_cpm0 flag_gene_on_Mock_cpm0;
run;

%mend;

%flagOnOff(uniq);
%flagOnOff(multi);

/*
     flag_gene_on_                             Cumulative    Cumulative
         01gy_cpm0    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       27900       96.79         27900        96.79
                 1         925        3.21         28825       100.00

                         Frequency Missing = 2364


     flag_gene_on_                             Cumulative    Cumulative
          1gy_cpm0    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       27765       96.69         27765        96.69
                 1         951        3.31         28716       100.00

                         Frequency Missing = 2473


     flag_gene_on_                             Cumulative    Cumulative
         mock_cpm0    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0       27909       96.74         27909        96.74
                 1         941        3.26         28850       100.00

                         Frequency Missing = 2339




    flag_gene_on_                             Cumulative    Cumulative
        01gy_cpm0    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18727       88.25         18727        88.25
                1        2494       11.75         21221       100.00

                        Frequency Missing = 9968


    flag_gene_on_                             Cumulative    Cumulative
         1gy_cpm0    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18333       87.53         18333        87.53
                1        2613       12.47         20946       100.00

                       Frequency Missing = 10243


    flag_gene_on_                             Cumulative    Cumulative
        mock_cpm0    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       18363       87.48         18363        87.48
                1        2627       12.52         20990       100.00

                       Frequency Missing = 10199


*/


/* (2) MODELS AND FDR */


%macro runModels(countType);

data on_gene;
   set arabRNA.arab_flag_TE_&countType._on_cpm_gt0_2;
   if flag_gene_on_mock_cpm0=1 or flag_gene_on_01gy_cpm0=1 or flag_gene_on_1gy_cpm0=1 ;
   keep gene_id;
run;


/* Merge in with count data */
 

proc sort data=on_gene;
   by gene_id;
proc sort data=arabRNA.cpm_by_te_&countType.;
   by gene_id;
run;

data gene_counts_on;
   merge arabRNA.cpm_by_te_&countType. (in=in1) on_gene (in=in2);
   by gene_id;
   if in1 and in2;
   log_cpm = log(cpm+1);
run;

/* Merge in design file
   Want variables: treatment, time */

/* Model for contrasts */

*Cat together genotype, time and treatment, then find order for contrasts;

data gene_counts_on2;
  set gene_counts_on;
  where flag_TE=1;
  length time_trt $20.;
  time_trt = catx('_', time, treatment);
run;

proc freq data=gene_counts_on2 noprint;
  tables time_trt / out=order;
run;

proc print data=order;
run;


*12 levels;
/* ORDER IS:

1_0.1gy
1_1gy
1_Mock
24_0.1gy
24_1gy
24_Mock
3_0.1gy
3_1gy
3_Mock
72_0.1gy
72_1gy
72_Mock

*/

proc sort data=gene_counts_on2;
   by gene_id time treatment;
   run;

ods listing close;

proc mixed data=gene_counts_on2;
  by gene_id;
  class time treatment time_trt;
  model log_cpm = time_trt / htype=1;

  lsmeans time_trt;
                  /* Order of contrasts:     1hr              24hr             3h               72h  
                                             0.1gy 1gy Mock   0.1gy 1gy Mock   0.1gy 1gy Mock   0.1gy 1gy Mock */
  contrast '0.1gy-Mock: 1h'  time_trt        1     0    -1    0     0   0      0     0   0      0     0   0  ;
  contrast '1gy-Mock: 1h'    time_trt        0     1    -1    0     0   0      0     0   0      0     0   0  ;
  contrast '0.1gy-Mock: 24h' time_trt        0     0     0    1     0   -1     0     0   0      0     0   0  ;
  contrast '1gy-Mock: 24h'   time_trt        0     0     0    0     1   -1     0     0   0      0     0   0  ;
   contrast '0.1gy-Mock: 3h'  time_trt        0     0     0    0     0   0      1     0  -1      0     0   0     ;   
  contrast '1gy-Mock: 3h'    time_trt        0     0     0    0     0   0      0     1  -1      0     0   0     ;   

  contrast '0.1gy-Mock: 72h' time_trt        0     0     0    0     0   0      0     0   0      1     0  -1    ;
  contrast '1gy-Mock: 72h'   time_trt        0     0     0    0     0   0      0     0   0      0     1  -1  ;
 estimate '0.1gy-Mock: 1h'  time_trt        1     0    -1    0     0   0      0     0   0      0     0   0  ;
  estimate '1gy-Mock: 1h'    time_trt        0     1    -1    0     0   0      0     0   0      0     0   0  ;

  estimate '0.1gy-Mock: 24h' time_trt        0     0     0    1     0   -1     0     0   0      0     0   0  ;
  estimate '1gy-Mock: 24h'   time_trt        0     0     0    0     1   -1     0     0   0      0     0   0  ;

  estimate '0.1gy-Mock: 3h'  time_trt        0     0     0    0     0   0      1     0  -1      0     0   0     ;   
  estimate '1gy-Mock: 3h'    time_trt        0     0     0    0     0   0      0     1  -1      0     0   0     ;   

  estimate '0.1gy-Mock: 72h' time_trt        0     0     0    0     0   0      0     0   0      1     0  -1    ;
  estimate '1gy-Mock: 72h'   time_trt        0     0     0    0     0   0      0     0   0      0     1  -1  ;

ods output tests1=tests1 lsmeans=lsmeans contrasts=contrasts estimates=estimates ;
run;
quit;


/* Make permenant */


data arabRNA.arab_TE_cntrs_tests1_&countType._2 ;
  set tests1 ;
  run ;

data arabRNA.arab_TE_cntrs_lsmeans_&countType._2 ;
  set lsmeans ;
  run ;

data arabRNA.arab_TE_cntrs_constr_&countType._2 ;
  set contrasts ;
  run ;

data arabRNA.arab_TE_cntrs_estim_&countType._2 ;
  set estimates ;
  run ;


%mend;

%runModels(uniq);
%runModels(multi);
   
%macro calcFDR(countType,onflag1,onflag2,labelname,outname,fdrname);

data on_gene;
   set arabRNA.arab_flag_TE_&countType._on_cpm_gt0_2;
   where &onflag1.=1 and &onflag2.=1 ;
   keep gene_id;
run;

data contrast_p;
  set arabRNA.arab_TE_cntrs_constr_&countType._2;
  format ProbF best32. ;
  where label=&labelname.;
  keep gene_id probf;
run;  

proc sort data=contrast_p;
  by gene_id;
proc sort data=on_gene;
  by gene_id;
run;

data contrast_p2;
  merge on_gene (in=in1) contrast_p (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc multtest inpvalues(ProbF)=contrast_p2 fdr
 out=fdr noprint;
run;
quit;

data fdr_&countType._&outname.;
  set fdr;
  if fdr_p = . then flag_&fdrname._fdr05=.;
     else if fdr_p <0.05 then flag_&fdrname._fdr05=1;
     else flag_&fdrname._fdr05=0;
  keep gene_id ProbF fdr_p flag_&fdrname._fdr05;
  rename fdr_p=fdr_&fdrname. ProbF=p_&fdrname.;
run;

proc freq data=fdr_&countType._&outname.;
   tables flag_&fdrname._fdr05 ;
run;

%mend;


%calcFDR(uniq,flag_gene_on_01gy_1hr_cpm0,flag_gene_on_Mock_1hr_cpm0,"0.1gy-Mock: 1h",contrast1,01gy_v_Mock_1h);
%calcFDR(uniq,flag_gene_on_01gy_3hr_cpm0,flag_gene_on_Mock_3hr_cpm0,"0.1gy-Mock: 3h",contrast2,01gy_v_Mock_3h);
%calcFDR(uniq,flag_gene_on_01gy_24hr_cpm0,flag_gene_on_Mock_24hr_cpm0,"0.1gy-Mock: 24h",contrast3,01gy_v_Mock_24h);
%calcFDR(uniq,flag_gene_on_01gy_72hr_cpm0,flag_gene_on_Mock_72hr_cpm0,"0.1gy-Mock: 72h",contrast4,01gy_v_Mock_72h);
%calcFDR(uniq,flag_gene_on_1gy_1hr_cpm0,flag_gene_on_Mock_1hr_cpm0,"1gy-Mock: 1h",contrast5,1gy_v_Mock_1h);
%calcFDR(uniq,flag_gene_on_1gy_3hr_cpm0,flag_gene_on_Mock_3hr_cpm0,"1gy-Mock: 3h",contrast6,1gy_v_Mock_3h);
%calcFDR(uniq,flag_gene_on_1gy_24hr_cpm0,flag_gene_on_Mock_24hr_cpm0,"1gy-Mock: 24h",contrast7,1gy_v_Mock_24h);
%calcFDR(uniq,flag_gene_on_1gy_72hr_cpm0,flag_gene_on_Mock_72hr_cpm0,"1gy-Mock: 72h",contrast8,1gy_v_Mock_72h);


%calcFDR(multi,flag_gene_on_01gy_1hr_cpm0,flag_gene_on_Mock_1hr_cpm0,"0.1gy-Mock: 1h",contrast1,01gy_v_Mock_1h);
%calcFDR(multi,flag_gene_on_01gy_3hr_cpm0,flag_gene_on_Mock_3hr_cpm0,"0.1gy-Mock: 3h",contrast2,01gy_v_Mock_3h);
%calcFDR(multi,flag_gene_on_01gy_24hr_cpm0,flag_gene_on_Mock_24hr_cpm0,"0.1gy-Mock: 24h",contrast3,01gy_v_Mock_24h);
%calcFDR(multi,flag_gene_on_01gy_72hr_cpm0,flag_gene_on_Mock_72hr_cpm0,"0.1gy-Mock: 72h",contrast4,01gy_v_Mock_72h);
%calcFDR(multi,flag_gene_on_1gy_1hr_cpm0,flag_gene_on_Mock_1hr_cpm0,"1gy-Mock: 1h",contrast5,1gy_v_Mock_1h);
%calcFDR(multi,flag_gene_on_1gy_3hr_cpm0,flag_gene_on_Mock_3hr_cpm0,"1gy-Mock: 3h",contrast6,1gy_v_Mock_3h);
%calcFDR(multi,flag_gene_on_1gy_24hr_cpm0,flag_gene_on_Mock_24hr_cpm0,"1gy-Mock: 24h",contrast7,1gy_v_Mock_24h);
%calcFDR(multi,flag_gene_on_1gy_72hr_cpm0,flag_gene_on_Mock_72hr_cpm0,"1gy-Mock: 72h",contrast8,1gy_v_Mock_72h);

proc sort data=fdr_uniq_contrast1; by gene_id;
proc sort data=fdr_uniq_contrast2; by gene_id;
proc sort data=fdr_uniq_contrast3; by gene_id;
proc sort data=fdr_uniq_contrast4; by gene_id;
proc sort data=fdr_uniq_contrast5; by gene_id;
proc sort data=fdr_uniq_contrast6; by gene_id;
proc sort data=fdr_uniq_contrast7; by gene_id;
proc sort data=fdr_uniq_contrast8; by gene_id;
run;

proc sort data=fdr_multi_contrast1; by gene_id;
proc sort data=fdr_multi_contrast2; by gene_id;
proc sort data=fdr_multi_contrast3; by gene_id;
proc sort data=fdr_multi_contrast4; by gene_id;
proc sort data=fdr_multi_contrast5; by gene_id;
proc sort data=fdr_multi_contrast6; by gene_id;
proc sort data=fdr_multi_contrast7; by gene_id;
proc sort data=fdr_multi_contrast8; by gene_id;
run;



data fdr_all_uniq;
  merge fdr_uniq_contrast1 fdr_uniq_contrast2 fdr_uniq_contrast3 fdr_uniq_contrast4
        fdr_uniq_contrast5 fdr_uniq_contrast6 fdr_uniq_contrast7 fdr_uniq_contrast8 ;
  by gene_id;
run;


data fdr_all_multi;
  merge fdr_multi_contrast1 fdr_multi_contrast2 fdr_multi_contrast3 fdr_multi_contrast4
        fdr_multi_contrast5 fdr_multi_contrast6 fdr_multi_contrast7 fdr_multi_contrast8 ;
  by gene_id;
run;


data arabRNA.arab_TE_fdr_uniq_2 ;
   set fdr_all_uniq;
run;

data arabRNA.arab_TE_fdr_multi_2 ;
   set fdr_all_multi;
run;



/* (3) CALC FOLDCHANGE */


%macro calcFC(countType);


data counts_w_key;
  set arabRNA.cpm_by_te_&countType.;
  where flag_TE=1;
run;

/* Flag on per treatment */

proc sort data=counts_w_key;
   by treatment time gene_id;
proc means data=counts_w_key noprint;
   by treatment time gene_id;
   var cpm;
   output out=mean_cpm mean=;
run;

data mean_cpm2;
  set mean_cpm;
  trt_time=catx('_',treatment,time);
run;

proc sort data=mean_cpm2;
  by gene_id  trt_time;
proc transpose data=mean_cpm2 out=mean_cpm_sbys;
   by gene_id ;
   var cpm;
   id trt_time;
run;

data logFC;
  set mean_cpm_sbys;
  log2FC_01Gy_1h=log2(_0_1gy_1 +1) - log2(Mock_1 +1);
  log2FC_01Gy_3h=log2(_0_1gy_3 +1) - log2(Mock_3 +1);
  log2FC_01Gy_24h=log2(_0_1gy_24 +1) - log2(Mock_24 +1);
  log2FC_01Gy_72h=log2(_0_1gy_72 +1) - log2(Mock_72 +1);

  log2FC_1Gy_1h=log2(_1gy_1 +1) - log2(Mock_1 +1);
  log2FC_1Gy_3h=log2(_1gy_3 +1) - log2(Mock_3 +1);
  log2FC_1Gy_24h=log2(_1gy_24 +1) - log2(Mock_24 +1);
  log2FC_1Gy_72h=log2(_1gy_72 +1) - log2(Mock_72 +1);
  keep gene_id log2FC_: ;
run;

data arabRNA.logFC_by_te_&countType.;
  set logFC;
run;

%mend;


%calcFC(uniq);
%calcFC(multi);

/* (4) COUNT SIG BY TIME AND DOSE */

data fdr_uniq;
   set arabRNA.arab_te_fdr_uniq_2;
run;

data fdr_multi;
   set arabRNA.arab_te_fdr_multi_2;
run;

data fc_uniq;
  set arabRNA.logFC_by_te_uniq;
run;

data fc_multi;
  set arabRNA.logFC_by_te_multi;
run;

proc sort data=fdr_uniq;
  by gene_id;
proc sort data=fdr_multi;
  by gene_id;
proc sort data=fc_uniq;
  by gene_id;
proc sort data=fc_multi;
  by gene_id;
run;

data fdr_fc_uniq;
  merge fdr_uniq (in=in1) fc_uniq (in=in2);
  by gene_id;
  if in1 and in2;
run;


data fdr_fc_multi;
  merge fdr_multi (in=in1) fc_multi (in=in2);
  by gene_id;
  if in1 and in2;
run;

data uniq_01_1_up uniq_01_3_up uniq_01_24_up uniq_01_72_up   
     uniq_01_1_dn uniq_01_3_dn uniq_01_24_dn uniq_01_72_dn   
     uniq_1_1_up  uniq_1_3_up  uniq_1_24_up  uniq_1_72_up
     uniq_1_1_dn  uniq_1_3_dn  uniq_1_24_dn  uniq_1_72_dn;
   set fdr_fc_uniq;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h <= -1 then output uniq_01_1_dn;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h >=  1 then output uniq_01_1_up;

   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h <= -1 then output uniq_01_3_dn;
   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h >=  1 then output uniq_01_3_up;

   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h <= -1 then output uniq_01_24_dn;
   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h >=  1 then output uniq_01_24_up;

   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h <= -1 then output uniq_01_72_dn;
   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h >=  1 then output uniq_01_72_up;

   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h <= -1 then output uniq_1_1_dn;
   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h >=  1 then output uniq_1_1_up;

   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h <= -1 then output uniq_1_3_dn;
   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h >=  1 then output uniq_1_3_up;

   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h <= -1 then output uniq_1_24_dn;
   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h >=  1 then output uniq_1_24_up;

   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h <= -1 then output uniq_1_72_dn;
   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h >=  1 then output uniq_1_72_up;
run;


data multi_01_1_up multi_01_3_up multi_01_24_up multi_01_72_up   
     multi_01_1_dn multi_01_3_dn multi_01_24_dn multi_01_72_dn   
     multi_1_1_up  multi_1_3_up  multi_1_24_up  multi_1_72_up
     multi_1_1_dn  multi_1_3_dn  multi_1_24_dn  multi_1_72_dn;
   set fdr_fc_multi;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h <= -1 then output multi_01_1_dn;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h >=  1 then output multi_01_1_up;

   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h <= -1 then output multi_01_3_dn;
   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h >=  1 then output multi_01_3_up;

   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h <= -1 then output multi_01_24_dn;
   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h >=  1 then output multi_01_24_up;

   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h <= -1 then output multi_01_72_dn;
   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h >=  1 then output multi_01_72_up;

   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h <= -1 then output multi_1_1_dn;
   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h >=  1 then output multi_1_1_up;

   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h <= -1 then output multi_1_3_dn;
   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h >=  1 then output multi_1_3_up;

   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h <= -1 then output multi_1_24_dn;
   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h >=  1 then output multi_1_24_up;

   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h <= -1 then output multi_1_72_dn;
   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h >=  1 then output multi_1_72_up;
run;

/*
UNIQ:
    10cGy  100cGy
    UP  DN  UP  DN
1h  0   0   0   2
3h  0   0   0   0
24h 1   0   0   0
72h 1   0   0   1


MULTI:
    10cGy   100cGy
    UP  DN
1h  0   0   0   1
3h  0   0   0   0
24h 1   0   0   0
72h 1   0   0   1

*/


data uniq_01_1_up uniq_01_3_up uniq_01_24_up uniq_01_72_up   
     uniq_01_1_dn uniq_01_3_dn uniq_01_24_dn uniq_01_72_dn   
     uniq_1_1_up  uniq_1_3_up  uniq_1_24_up  uniq_1_72_up
     uniq_1_1_dn  uniq_1_3_dn  uniq_1_24_dn  uniq_1_72_dn;
   set fdr_fc_uniq;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h < 0 then output uniq_01_1_dn;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h > 0 then output uniq_01_1_up;

   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h < 0 then output uniq_01_3_dn;
   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h > 0 then output uniq_01_3_up;

   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h < 0 then output uniq_01_24_dn;
   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h > 0 then output uniq_01_24_up;

   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h < 0 then output uniq_01_72_dn;
   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h > 0 then output uniq_01_72_up;

   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h < 0 then output uniq_1_1_dn;
   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h > 0 then output uniq_1_1_up;

   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h < 0 then output uniq_1_3_dn;
   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h > 0 then output uniq_1_3_up;

   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h < 0 then output uniq_1_24_dn;
   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h > 0 then output uniq_1_24_up;

   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h < 0 then output uniq_1_72_dn;
   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h > 0 then output uniq_1_72_up;
run;


data multi_01_1_up multi_01_3_up multi_01_24_up multi_01_72_up   
     multi_01_1_dn multi_01_3_dn multi_01_24_dn multi_01_72_dn   
     multi_1_1_up  multi_1_3_up  multi_1_24_up  multi_1_72_up
     multi_1_1_dn  multi_1_3_dn  multi_1_24_dn  multi_1_72_dn;
   set fdr_fc_multi;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h < 0 then output multi_01_1_dn;
   if flag_01gy_v_Mock_1h_fdr05 and log2FC_01Gy_1h > 0 then output multi_01_1_up;

   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h < 0 then output multi_01_3_dn;
   if flag_01gy_v_Mock_3h_fdr05 and log2FC_01Gy_3h > 0 then output multi_01_3_up;

   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h < 0 then output multi_01_24_dn;
   if flag_01gy_v_Mock_24h_fdr05 and log2FC_01Gy_24h > 0 then output multi_01_24_up;

   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h < 0 then output multi_01_72_dn;
   if flag_01gy_v_Mock_72h_fdr05 and log2FC_01Gy_72h > 0 then output multi_01_72_up;

   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h < 0 then output multi_1_1_dn;
   if flag_1gy_v_Mock_1h_fdr05 and log2FC_1Gy_1h > 0 then output multi_1_1_up;

   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h < 0 then output multi_1_3_dn;
   if flag_1gy_v_Mock_3h_fdr05 and log2FC_1Gy_3h > 0 then output multi_1_3_up;

   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h < 0 then output multi_1_24_dn;
   if flag_1gy_v_Mock_24h_fdr05 and log2FC_1Gy_24h > 0 then output multi_1_24_up;

   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h < 0 then output multi_1_72_dn;
   if flag_1gy_v_Mock_72h_fdr05 and log2FC_1Gy_72h > 0 then output multi_1_72_up;
run;

/*
UNIQ:
    10cGy  100cGy
    UP  DN  UP  DN
1h  12   6   8   2
3h  0   0   11   1
24h 1   0   0   0
72h 4   0   1   0


MULTI:
    10cGy   100cGy
    UP  DN
1h  9   3   4   0
3h  0   1   4   0
24h 1   0   0   0
72h 5   0   0   0

*/

