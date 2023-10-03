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
  set arabRNA.cpm_by_te_local_&countType.;
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




data arabRNA.arab_flag_TE_&countType._on_cpm_gt0;
  merge treat_cpm0 (in=in1) time_cpm0 (in=in2) trttime_cpm0 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

data arabRNA.arab_flag_TE_&countType._on_cpm_ge5;
  merge treat_cpm5 (in=in1) time_cpm5 (in=in2) trttime_cpm5 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

/* Count CPM > 0 */

ods listing;
proc freq data=arabRNA.arab_flag_TE_&countType._on_cpm_gt0;
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
   set arabRNA.arab_flag_TE_&countType._on_cpm_gt0;
   if flag_gene_on_mock_cpm0=1 or flag_gene_on_01gy_cpm0=1 or flag_gene_on_1gy_cpm0=1 ;
   keep gene_id;
run;


/* Merge in with count data */
 

proc sort data=on_gene;
   by gene_id;
proc sort data=arabRNA.cpm_by_te_local_&countType.;
   by gene_id;
run;

data gene_counts_on;
   merge arabRNA.cpm_by_te_local_&countType. (in=in1) on_gene (in=in2);
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


data arabRNA.arab_TE_cntrs_tests1_&countType. ;
  set tests1 ;
  run ;

data arabRNA.arab_TE_cntrs_lsmeans_&countType. ;
  set lsmeans ;
  run ;

data arabRNA.arab_TE_cntrs_constr_&countType. ;
  set contrasts ;
  run ;

data arabRNA.arab_TE_cntrs_estim_&countType. ;
  set estimates ;
  run ;


%mend;

%runModels(uniq);
%runModels(multi);
   
%macro calcFDR(countType,onflag1,onflag2,labelname,outname,fdrname);

data on_gene;
   set arabRNA.arab_flag_TE_&countType._on_cpm_gt0;
   where &onflag1.=1 and &onflag2.=1 ;
   keep gene_id;
run;

data contrast_p;
  set arabRNA.arab_TE_cntrs_constr_&countType.;
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


data arabRNA.arab_TE_fdr_uniq ;
   set fdr_all_uniq;
run;

data arabRNA.arab_TE_fdr_multi ;
   set fdr_all_multi;
run;



/* (3) CALC FOLDCHANGE */


%macro calcFC(countType);


data counts_w_key;
  set arabRNA.cpm_by_te_local_&countType.;
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

data arabRNA.logFC_by_te_local_&countType.;
  set logFC;
run;

%mend;


%calcFC(uniq);
%calcFC(multi);

/* (4) COUNT SIG BY TIME AND DOSE */

data fdr_uniq;
   set arabRNA.arab_te_fdr_uniq;
run;

data fdr_multi;
   set arabRNA.arab_te_fdr_multi;
run;

data fc_uniq;
  set arabRNA.logFC_by_te_local_uniq;
run;

data fc_multi;
  set arabRNA.logFC_by_te_local_multi;
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
   length geneID $50.;
   geneID=compress(scan(gene_id,1,":"));
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

   keep geneID;

run;


data multi_01_1_up multi_01_3_up multi_01_24_up multi_01_72_up   
     multi_01_1_dn multi_01_3_dn multi_01_24_dn multi_01_72_dn   
     multi_1_1_up  multi_1_3_up  multi_1_24_up  multi_1_72_up
     multi_1_1_dn  multi_1_3_dn  multi_1_24_dn  multi_1_72_dn;
   set fdr_fc_multi;
   length geneID $50.;
   geneID=compress(scan(gene_id,1,":"));
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

   keep gene_id geneID;

run;

ods listing;
proc print data=multi_01_1_up;
proc print data=multi_01_1_dn;
run;


AT1TE04835:ATGP1:Gypsy:LTR
AT1TE04980:ATREP3:Helitron:RC
AT1TE10830:ATLINE2:L1:LINE
AT1TE41790:SADHU:SADHU:Unassigned
AT1TE54835:ATCOPIA25:Copia:LTR
AT1TE71770:ATHATN3:HAT:DNA
AT1TE97605:HELITRONY1B:Helitron:RC
AT1TE99345:BRODYAGA1A:MuDR:DNA
AT2TE42485:ATCOPIA75:Copia:LTR
AT2TE57160:Unassigned:Unassigned:Unassigned
AT3TE05245:Unassigned:Unassigned:Unassigned
AT3TE91850:ATLINEIII:L1:LINE
AT4TE42860:ATCOPIA4:Copia:LTR
AT5TE48385:ATREP10D:Helitron:RC
AT5TE60290:ATLINE1_6:L1:LINE
AT5TE80015:ATREP3:Helitron:RC
AT5TE81425:Unassigned:Unassigned:Unassigned

  AT1TE00225:HELITRONY1D:Helitron:RC
  AT1TE59005:BRODYAGA1A:MuDR:DNA
  AT2TE11255:VANDAL5A:MuDR:DNA
  AT3TE43425:ATHPOGON2:Pogo:DNA
  AT5TE29585:ATREP5:Helitron:RC
  AT5TE29590:SIMPLEGUY1:Harbinger:DNA


proc print data=multi_01_24_up;
proc print data=multi_01_24_dn;
run;


 AT2TE11560:ATHILA6A:Gypsy:LTR
 AT5TE36215:HELITRONY1B:Helitron:RC
 AT3TE67135:TNAT1A:TNAT1A:DNA



proc print data=multi_01_72_up;
proc print data=multi_01_72_dn;
run;

 AT5TE36215:HELITRONY1B:Helitron:RC
 AT5TE36220:ATREP5:Helitron:RC




proc print data=multi_1_1_up;
proc print data=multi_1_1_dn;
run;

 AT1TE10830:ATLINE2:L1:LINE
 AT1TE54835:ATCOPIA25:Copia:LTR
 AT1TE71770:ATHATN3:HAT:DNA
 AT2TE42485:ATCOPIA75:Copia:LTR
 AT4TE04415:ATIS112A:Harbinger:DNA
 AT5TE37615:ATREP3:Helitron:RC
 AT5TE38720:SADHU:SADHU:Unassigned
 AT5TE48385:ATREP10D:Helitron:RC
 AT5TE81425:Unassigned:Unassigned:Unassigned
AT1TE59005:BRODYAGA1A:MuDR:DNA
AT4TE39650:BRODYAGA1A:MuDR:DNA
AT5TE29585:ATREP5:Helitron:RC
AT5TE29590:SIMPLEGUY1:Harbinger:DNA






proc print data=multi_1_72_up;
proc print data=multi_1_72_dn;
run;

 AT5TE36215:HELITRONY1B:Helitron:RC
 AT2TE24525:ATENSPM1A:En-Spm:DNA



 
/*
UNIQ:
    10cGy  100cGy
    UP  DN  UP  DN
1h  24   8   14   5
3h  0   0   2   3
24h 1   1   0   0
72h 5   0   1   1


MULTI:
    10cGy   100cGy
    UP  DN
1h  17   6   9   4
3h  0   0   0   0
24h 2   1   0   0
72h 2   0   1   1

*/



libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

data cg_up_dmr_01_72 cg_up_dmr_1_72 cg_dn_dmr_01_72 cg_dn_dmr_1_72
     chg_up_dmr_01_72 chg_up_dmr_1_72 chg_dn_dmr_01_72 chg_dn_dmr_1_72
     chh_up_dmr_01_72 chh_up_dmr_1_72 chh_dn_dmr_01_72 chh_dn_dmr_1_72;
     set arabMAP.results_by_dmr_annot_TE;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     if feature = "Intergenic" then delete;
     /* FET */
     if site_type="CG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72;
     end;


     /* FET */
     if site_type="CHG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72;
     end;

     /* FET */
     if site_type="CHH" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72;
     end;


     keep  geneID;
run;

proc sort data=cg_up_dmr_01_72 nodup; by _all_; run;
proc sort data=cg_dn_dmr_01_72 nodup; by _all_; run;
proc sort data=cg_up_dmr_1_72 nodup; by _all_; run;
proc sort data=cg_dn_dmr_1_72 nodup; by _all_; run;

proc sort data=chg_up_dmr_01_72 nodup; by _all_; run;
proc sort data=chg_dn_dmr_01_72 nodup; by _all_; run;
proc sort data=chg_up_dmr_1_72 nodup; by _all_; run;
proc sort data=chg_dn_dmr_1_72 nodup; by _all_; run;

proc sort data=chh_up_dmr_01_72 nodup; by _all_; run;
proc sort data=chh_dn_dmr_01_72 nodup; by _all_; run;
proc sort data=chh_up_dmr_1_72 nodup; by _all_; run;
proc sort data=chh_dn_dmr_1_72 nodup; by _all_; run;




data up_dar_01_72 up_dar_1_72 dn_dar_01_72 dn_dar_1_72;
     set arabMAP.results_by_dar_annot_TE;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     if feature = "Intergenic" then delete;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;
     keep geneID;
run;



proc sort data=up_dar_01_72 nodup; by _all_; run;
proc sort data=dn_dar_01_72 nodup; by _all_; run;
proc sort data=up_dar_1_72 nodup; by _all_; run;
proc sort data=dn_dar_1_72 nodup; by _all_; run;


%macro compareGenes(DEG,DMG);

proc sort data=&DEG.;
  by geneID;
proc sort data=&DMG.;
  by geneID;
run;


data compare;
  merge &DEG. (in=in1) &DMG. (in=in2);
  by geneID;
  if in1 and in2;
run;

%mend;


%compareGenes(uniq_01_1_up,cg_up_dmr_01_72); %compareGenes(uniq_01_1_up,cg_dn_dmr_01_72); %compareGenes(uniq_01_1_up,chg_up_dmr_01_72); %compareGenes(uniq_01_1_up,chg_dn_dmr_01_72);
%compareGenes(uniq_01_1_up,chh_up_dmr_01_72); %compareGenes(uniq_01_1_up,chh_dn_dmr_01_72); %compareGenes(uniq_01_1_up,up_dar_01_72); %compareGenes(uniq_01_1_up,dn_dar_01_72);

%compareGenes(uniq_01_3_up,cg_up_dmr_01_72); %compareGenes(uniq_01_3_up,cg_dn_dmr_01_72); %compareGenes(uniq_01_3_up,chg_up_dmr_01_72); %compareGenes(uniq_01_3_up,chg_dn_dmr_01_72); 
%compareGenes(uniq_01_3_up,chh_up_dmr_01_72); %compareGenes(uniq_01_3_up,chh_dn_dmr_01_72); %compareGenes(uniq_01_3_up,up_dar_01_72); %compareGenes(uniq_01_3_up,dn_dar_01_72); 

%compareGenes(uniq_01_24_up,cg_up_dmr_01_72); %compareGenes(uniq_01_24_up,cg_dn_dmr_01_72); %compareGenes(uniq_01_24_up,chg_up_dmr_01_72); %compareGenes(uniq_01_24_up,chg_dn_dmr_01_72); 
%compareGenes(uniq_01_24_up,chh_up_dmr_01_72); %compareGenes(uniq_01_24_up,chh_dn_dmr_01_72); %compareGenes(uniq_01_24_up,up_dar_01_72); %compareGenes(uniq_01_24_up,dn_dar_01_72); 

%compareGenes(uniq_01_72_up,cg_up_dmr_01_72 );   %compareGenes(uniq_01_72_up, cg_dn_dmr_01_72);   %compareGenes(uniq_01_72_up,chg_up_dmr_01_72);   %compareGenes(uniq_01_72_up,chg_dn_dmr_01_72 );   
%compareGenes(uniq_01_72_up,chh_up_dmr_01_72 );   %compareGenes(uniq_01_72_up,chh_dn_dmr_01_72 );   %compareGenes(uniq_01_72_up,up_dar_01_72);   %compareGenes(uniq_01_72_up,dn_dar_01_72 );   

%compareGenes(uniq_01_1_dn,cg_up_dmr_01_72); %compareGenes(uniq_01_1_dn,cg_dn_dmr_01_72); %compareGenes(uniq_01_1_dn,chg_up_dmr_01_72); %compareGenes(uniq_01_1_dn,chg_dn_dmr_01_72); 
%compareGenes(uniq_01_1_dn,chh_up_dmr_01_72); %compareGenes(uniq_01_1_dn,chh_dn_dmr_01_72); %compareGenes(uniq_01_1_dn,up_dar_01_72); %compareGenes(uniq_01_1_dn,dn_dar_01_72); 

%compareGenes(uniq_01_3_dn,cg_up_dmr_01_72); %compareGenes(uniq_01_3_dn,cg_dn_dmr_01_72); %compareGenes(uniq_01_3_dn,chg_up_dmr_01_72); %compareGenes(uniq_01_3_dn,chg_dn_dmr_01_72); 
%compareGenes(uniq_01_3_dn,chh_up_dmr_01_72); %compareGenes(uniq_01_3_dn,chh_dn_dmr_01_72); %compareGenes(uniq_01_3_dn,up_dar_01_72); %compareGenes(uniq_01_3_dn,dn_dar_01_72); 

%compareGenes(uniq_01_24_dn,cg_up_dmr_01_72); %compareGenes(uniq_01_24_dn,cg_dn_dmr_01_72); %compareGenes(uniq_01_24_dn,chg_up_dmr_01_72); %compareGenes(uniq_01_24_dn,chg_dn_dmr_01_72); %compareGenes(uniq_01_24_dn,chh_up_dmr_01_72); %compareGenes(uniq_01_24_dn,chh_dn_dmr_01_72); %compareGenes(uniq_01_24_dn,up_dar_01_72); %compareGenes(uniq_01_24_dn,dn_dar_01_72); 

%compareGenes(uniq_01_72_dn,cg_up_dmr_01_72); %compareGenes(uniq_01_72_dn,cg_dn_dmr_01_72); %compareGenes(uniq_01_72_dn,chg_up_dmr_01_72); %compareGenes(uniq_01_72_dn,chg_dn_dmr_01_72); %compareGenes(uniq_01_72_dn,chh_up_dmr_01_72); %compareGenes(uniq_01_72_dn,chh_dn_dmr_01_72); %compareGenes(uniq_01_72_dn,up_dar_01_72); %compareGenes(uniq_01_72_dn,dn_dar_01_72); 

%compareGenes(uniq_1_1_up,cg_up_dmr_1_72); %compareGenes(uniq_1_1_up,cg_dn_dmr_1_72); %compareGenes(uniq_1_1_up,chg_up_dmr_1_72); %compareGenes(uniq_1_1_up,chg_dn_dmr_1_72); %compareGenes(uniq_1_1_up,chh_up_dmr_1_72); %compareGenes(uniq_1_1_up,chh_dn_dmr_1_72); %compareGenes(uniq_1_1_up,up_dar_1_72); %compareGenes(uniq_1_1_up,dn_dar_1_72); 

%compareGenes(uniq_1_3_up,cg_up_dmr_1_72); %compareGenes(uniq_1_3_up,cg_dn_dmr_1_72); %compareGenes(uniq_1_3_up,chg_up_dmr_1_72); %compareGenes(uniq_1_3_up,chg_dn_dmr_1_72); %compareGenes(uniq_1_3_up,chh_up_dmr_1_72); %compareGenes(uniq_1_3_up,chh_dn_dmr_1_72); %compareGenes(uniq_1_3_up,up_dar_1_72); %compareGenes(uniq_1_3_up,dn_dar_1_72); 

%compareGenes(uniq_1_24_up,cg_up_dmr_1_72); %compareGenes(uniq_1_24_up,cg_dn_dmr_1_72); %compareGenes(uniq_1_24_up,chg_up_dmr_1_72); %compareGenes(uniq_1_24_up,chg_dn_dmr_1_72); %compareGenes(uniq_1_24_up,chh_up_dmr_1_72); %compareGenes(uniq_1_24_up,chh_dn_dmr_1_72); %compareGenes(uniq_1_24_up,up_dar_1_72); %compareGenes(uniq_1_24_up,dn_dar_1_72); 

%compareGenes(uniq_1_72_up,cg_up_dmr_1_72); %compareGenes(uniq_1_72_up,cg_dn_dmr_1_72); %compareGenes(uniq_1_72_up,chg_up_dmr_1_72); %compareGenes(uniq_1_72_up,chg_dn_dmr_1_72); %compareGenes(uniq_1_72_up,chh_up_dmr_1_72); %compareGenes(uniq_1_72_up,chh_dn_dmr_1_72); %compareGenes(uniq_1_72_up,up_dar_1_72); %compareGenes(uniq_1_72_up,dn_dar_1_72); 

%compareGenes(uniq_1_1_dn,cg_up_dmr_1_72); %compareGenes(uniq_1_1_dn,cg_dn_dmr_1_72); %compareGenes(uniq_1_1_dn,chg_up_dmr_1_72); %compareGenes(uniq_1_1_dn,chg_dn_dmr_1_72); %compareGenes(uniq_1_1_dn,chh_up_dmr_1_72); %compareGenes(uniq_1_1_dn,chh_dn_dmr_1_72); %compareGenes(uniq_1_1_dn,up_dar_1_72); %compareGenes(uniq_1_1_dn,dn_dar_1_72); 

%compareGenes(uniq_1_3_dn,cg_up_dmr_1_72); %compareGenes(uniq_1_3_dn,cg_dn_dmr_1_72); %compareGenes(uniq_1_3_dn,chg_up_dmr_1_72); %compareGenes(uniq_1_3_dn,chg_dn_dmr_1_72); %compareGenes(uniq_1_3_dn,chh_up_dmr_1_72); %compareGenes(uniq_1_3_dn,chh_dn_dmr_1_72); %compareGenes(uniq_1_3_dn,up_dar_1_72); %compareGenes(uniq_1_3_dn,dn_dar_1_72); 

%compareGenes(uniq_1_24_dn,cg_up_dmr_1_72); %compareGenes(uniq_1_24_dn,cg_dn_dmr_1_72); %compareGenes(uniq_1_24_dn,chg_up_dmr_1_72); %compareGenes(uniq_1_24_dn,chg_dn_dmr_1_72);  %compareGenes(uniq_1_24_dn,chh_up_dmr_1_72); %compareGenes(uniq_1_24_dn,chh_dn_dmr_1_72); %compareGenes(uniq_1_24_dn,up_dar_1_72); %compareGenes(uniq_1_24_dn,dn_dar_1_72); 

%compareGenes(uniq_1_72_dn,cg_up_dmr_1_72); %compareGenes(uniq_1_72_dn,cg_dn_dmr_1_72); %compareGenes(uniq_1_72_dn,chg_up_dmr_1_72); %compareGenes(uniq_1_72_dn,chg_dn_dmr_1_72); %compareGenes(uniq_1_72_dn,chh_up_dmr_1_72); %compareGenes(uniq_1_72_dn,chh_dn_dmr_1_72); %compareGenes(uniq_1_72_dn,up_dar_1_72); %compareGenes(uniq_1_72_dn,dn_dar_1_72); 



%compareGenes(multi_01_1_up,cg_up_dmr_01_72); %compareGenes(multi_01_1_up,cg_dn_dmr_01_72); %compareGenes(multi_01_1_up,chg_up_dmr_01_72); %compareGenes(multi_01_1_up,chg_dn_dmr_01_72);%compareGenes(multi_01_1_up,chh_up_dmr_01_72); %compareGenes(multi_01_1_up,chh_dn_dmr_01_72); %compareGenes(multi_01_1_up,up_dar_01_72); %compareGenes(multi_01_1_up,dn_dar_01_72);

%compareGenes(multi_01_3_up,cg_up_dmr_01_72); %compareGenes(multi_01_3_up,cg_dn_dmr_01_72); %compareGenes(multi_01_3_up,chg_up_dmr_01_72); %compareGenes(multi_01_3_up,chg_dn_dmr_01_72); %compareGenes(multi_01_3_up,chh_up_dmr_01_72); %compareGenes(multi_01_3_up,chh_dn_dmr_01_72); %compareGenes(multi_01_3_up,up_dar_01_72); %compareGenes(multi_01_3_up,dn_dar_01_72); 

%compareGenes(multi_01_24_up,cg_up_dmr_01_72); %compareGenes(multi_01_24_up,cg_dn_dmr_01_72); %compareGenes(multi_01_24_up,chg_up_dmr_01_72); %compareGenes(multi_01_24_up,chg_dn_dmr_01_72); %compareGenes(multi_01_24_up,chh_up_dmr_01_72); %compareGenes(multi_01_24_up,chh_dn_dmr_01_72); %compareGenes(multi_01_24_up,up_dar_01_72); %compareGenes(multi_01_24_up,dn_dar_01_72); 

%compareGenes(multi_01_72_up,cg_up_dmr_01_72 );   %compareGenes(multi_01_72_up, cg_dn_dmr_01_72);   %compareGenes(multi_01_72_up,chg_up_dmr_01_72);   %compareGenes(multi_01_72_up,chg_dn_dmr_01_72 );   %compareGenes(multi_01_72_up,chh_up_dmr_01_72 );   %compareGenes(multi_01_72_up,chh_dn_dmr_01_72 );   %compareGenes(multi_01_72_up,up_dar_01_72);   %compareGenes(multi_01_72_up,dn_dar_01_72 );   

%compareGenes(multi_01_1_dn,cg_up_dmr_01_72); %compareGenes(multi_01_1_dn,cg_dn_dmr_01_72); %compareGenes(multi_01_1_dn,chg_up_dmr_01_72); %compareGenes(multi_01_1_dn,chg_dn_dmr_01_72); %compareGenes(multi_01_1_dn,chh_up_dmr_01_72); %compareGenes(multi_01_1_dn,chh_dn_dmr_01_72); %compareGenes(multi_01_1_dn,up_dar_01_72); %compareGenes(multi_01_1_dn,dn_dar_01_72); 

%compareGenes(multi_01_3_dn,cg_up_dmr_01_72); %compareGenes(multi_01_3_dn,cg_dn_dmr_01_72); %compareGenes(multi_01_3_dn,chg_up_dmr_01_72); %compareGenes(multi_01_3_dn,chg_dn_dmr_01_72); %compareGenes(multi_01_3_dn,chh_up_dmr_01_72); %compareGenes(multi_01_3_dn,chh_dn_dmr_01_72); %compareGenes(multi_01_3_dn,up_dar_01_72); %compareGenes(multi_01_3_dn,dn_dar_01_72); 

%compareGenes(multi_01_24_dn,cg_up_dmr_01_72); %compareGenes(multi_01_24_dn,cg_dn_dmr_01_72); %compareGenes(multi_01_24_dn,chg_up_dmr_01_72); %compareGenes(multi_01_24_dn,chg_dn_dmr_01_72); %compareGenes(multi_01_24_dn,chh_up_dmr_01_72); %compareGenes(multi_01_24_dn,chh_dn_dmr_01_72); %compareGenes(multi_01_24_dn,up_dar_01_72); %compareGenes(multi_01_24_dn,dn_dar_01_72); 

%compareGenes(multi_01_72_dn,cg_up_dmr_01_72); %compareGenes(multi_01_72_dn,cg_dn_dmr_01_72); %compareGenes(multi_01_72_dn,chg_up_dmr_01_72); %compareGenes(multi_01_72_dn,chg_dn_dmr_01_72); %compareGenes(multi_01_72_dn,chh_up_dmr_01_72); %compareGenes(multi_01_72_dn,chh_dn_dmr_01_72); %compareGenes(multi_01_72_dn,up_dar_01_72); %compareGenes(multi_01_72_dn,dn_dar_01_72); 

%compareGenes(multi_1_1_up,cg_up_dmr_1_72); %compareGenes(multi_1_1_up,cg_dn_dmr_1_72); %compareGenes(multi_1_1_up,chg_up_dmr_1_72); %compareGenes(multi_1_1_up,chg_dn_dmr_1_72); %compareGenes(multi_1_1_up,chh_up_dmr_1_72); %compareGenes(multi_1_1_up,chh_dn_dmr_1_72); %compareGenes(multi_1_1_up,up_dar_1_72); %compareGenes(multi_1_1_up,dn_dar_1_72); 

%compareGenes(multi_1_3_up,cg_up_dmr_1_72); %compareGenes(multi_1_3_up,cg_dn_dmr_1_72); %compareGenes(multi_1_3_up,chg_up_dmr_1_72); %compareGenes(multi_1_3_up,chg_dn_dmr_1_72); %compareGenes(multi_1_3_up,chh_up_dmr_1_72); %compareGenes(multi_1_3_up,chh_dn_dmr_1_72); %compareGenes(multi_1_3_up,up_dar_1_72); %compareGenes(multi_1_3_up,dn_dar_1_72); 

%compareGenes(multi_1_24_up,cg_up_dmr_1_72); %compareGenes(multi_1_24_up,cg_dn_dmr_1_72); %compareGenes(multi_1_24_up,chg_up_dmr_1_72); %compareGenes(multi_1_24_up,chg_dn_dmr_1_72); %compareGenes(multi_1_24_up,chh_up_dmr_1_72); %compareGenes(multi_1_24_up,chh_dn_dmr_1_72); %compareGenes(multi_1_24_up,up_dar_1_72); %compareGenes(multi_1_24_up,dn_dar_1_72); 

%compareGenes(multi_1_72_up,cg_up_dmr_1_72); %compareGenes(multi_1_72_up,cg_dn_dmr_1_72); %compareGenes(multi_1_72_up,chg_up_dmr_1_72); %compareGenes(multi_1_72_up,chg_dn_dmr_1_72); %compareGenes(multi_1_72_up,chh_up_dmr_1_72); %compareGenes(multi_1_72_up,chh_dn_dmr_1_72); %compareGenes(multi_1_72_up,up_dar_1_72); %compareGenes(multi_1_72_up,dn_dar_1_72); 

%compareGenes(multi_1_1_dn,cg_up_dmr_1_72); %compareGenes(multi_1_1_dn,cg_dn_dmr_1_72); %compareGenes(multi_1_1_dn,chg_up_dmr_1_72); %compareGenes(multi_1_1_dn,chg_dn_dmr_1_72); %compareGenes(multi_1_1_dn,chh_up_dmr_1_72); %compareGenes(multi_1_1_dn,chh_dn_dmr_1_72); %compareGenes(multi_1_1_dn,up_dar_1_72); %compareGenes(multi_1_1_dn,dn_dar_1_72); 

%compareGenes(multi_1_3_dn,cg_up_dmr_1_72); %compareGenes(multi_1_3_dn,cg_dn_dmr_1_72); %compareGenes(multi_1_3_dn,chg_up_dmr_1_72); %compareGenes(multi_1_3_dn,chg_dn_dmr_1_72); %compareGenes(multi_1_3_dn,chh_up_dmr_1_72); %compareGenes(multi_1_3_dn,chh_dn_dmr_1_72); %compareGenes(multi_1_3_dn,up_dar_1_72); %compareGenes(multi_1_3_dn,dn_dar_1_72); 

%compareGenes(multi_1_24_dn,cg_up_dmr_1_72); %compareGenes(multi_1_24_dn,cg_dn_dmr_1_72); %compareGenes(multi_1_24_dn,chg_up_dmr_1_72); %compareGenes(multi_1_24_dn,chg_dn_dmr_1_72);  %compareGenes(multi_1_24_dn,chh_up_dmr_1_72); %compareGenes(multi_1_24_dn,chh_dn_dmr_1_72); %compareGenes(multi_1_24_dn,up_dar_1_72); %compareGenes(multi_1_24_dn,dn_dar_1_72); 

%compareGenes(multi_1_72_dn,cg_up_dmr_1_72); %compareGenes(multi_1_72_dn,cg_dn_dmr_1_72); %compareGenes(multi_1_72_dn,chg_up_dmr_1_72); %compareGenes(multi_1_72_dn,chg_dn_dmr_1_72); %compareGenes(multi_1_72_dn,chh_up_dmr_1_72); %compareGenes(multi_1_72_dn,chh_dn_dmr_1_72); %compareGenes(multi_1_72_dn,up_dar_1_72); %compareGenes(multi_1_72_dn,dn_dar_1_72); 





/*

UNIQ            TOTAL   CG-UP   CG-DN   CHG-UP  CHG-DN  CHH-UP  CHH-DN  GC-UP   GC-DN
 10cGy  1h UP   24      0       0       0       0       1       6       8       1
 10cGy  1h DN   8       0       0       0       0       1       3       2       0
 10cGy  3h UP   0       0       0       0       0       0       0       0       0
 10cGy  3h DN   0       0       0       0       0       0       0       0       0
 10cGy 24h UP   1       0       0       0       0       0       0       0       0
 10cGy 24h DN   1       0       0       0       0       0       0       0       0
 10cGy 72h UP   5       0       0       0       0       0       3       1       0
 10cGy 72h DN   0       0       0       0       0       0       0       0       0
100cGy  1h UP   14      0       0       0       0       3       1       4       0
100cGy  1h DN   5       0       0       0       0       1       1       0       0
100cGy  3h UP   2       0       0       0       0       0       0       0       0
100cGy  3h DN   3       0       0       0       0       1       1       2       0
100cGy 24h UP   0       0       0       0       0       0       0       0       0
100cGy 24h DN   0       0       0       0       0       0       0       0       0
100cGy 72h UP   1       0       0       0       0       0       0       0       0
100cGy 72h DN   1       0       0       0       0       1       0       1       0
    

MULTI:          TOTAL   CG-UP   CG-DN   CHG-UP  CHG-DN  CHH-UP  CHH-DN  GC-UP   GC-DN
 10cGy  1h UP   17      0       0       0       0       1       5       7       1
 10cGy  1h DN   6       0       0       0       0       0       4       1       0
 10cGy  3h UP   0       0       0       0       0       0       0       0       0
 10cGy  3h DN   0       0       0       0       0       0       0       0       0
 10cGy 24h UP   2       0       0       0       0       0       1       1       0
 10cGy 24h DN   1       0       0       0       0       0       0       0       0
 10cGy 72h UP   2       0       0       0       0       0       0       0       0
 10cGy 72h DN   0       0       0       0       0       0       0       0       0
100cGy  1h UP   9       0       0       0       0       2       0       2       0
100cGy  1h DN   4       0       0       0       0       1       1       0       0
100cGy  3h UP   0       0       0       0       0       0       0       0       0
100cGy  3h DN   0       0       0       0       0       0       0       0       0
100cGy 24h UP   0       0       0       0       0       0       0       0       0
100cGy 24h DN   0       0       0       0       0       0       0       0       0
100cGy 72h UP   1       0       0       0       0       0       0       0       0
100cGy 72h DN   1       0       0       0       0       1       0       1       0


*/











/*
UNIQ            TOTAL   CG-UP   CG-DN   CHG-UP  CHG-DN  CHH-UP  CHH-DN  GC-UP   GC-DN
 10cGy  1h UP   24      0       0       0       0       1       6       8       1
 10cGy  1h DN   8       0       0       0       0       1       3       2       0
 10cGy  3h UP   0       0       0       0       0       0       0       0       0
 10cGy  3h DN   0       0       0       0       0       0       0       0       0
 10cGy 24h UP   1       0       0       0       0       0       0       0       0
 10cGy 24h DN   1       0       0       0       0       0       0       0       0
 10cGy 72h UP   5       0       0       0       0       0       3       1       0
 10cGy 72h DN   0       0       0       0       0       0       0       0       0
100cGy  1h UP   14      0       0       0       0       3       1       4       0
100cGy  1h DN   5       0       0       0       0       1       1       0       0
100cGy  3h UP   2       0       0       0       0       0       0       0       0
100cGy  3h DN   3       0       0       0       0       1       1       2       0
100cGy 24h UP   0       0       0       0       0       0       0       0       0
100cGy 24h DN   0       0       0       0       0       0       0       0       0
100cGy 72h UP   1       0       0       0       0       0       0       0       0
100cGy 72h DN   1       0       0       0       0       1       0       1       0
    
    

MULTI:          TOTAL   CG-UP   CG-DN   CHG-UP  CHG-DN   CHH-UP  CHH-DN   GC-UP  GC-DN
 10cGy  1h UP   17
 10cGy  1h DN   6
 10cGy  3h UP   0       0       0       0       0       0       0       0       0
 10cGy  3h DN   0       0       0       0       0       0       0       0       0
 10cGy 24h UP   2
 10cGy 24h DN   1
 10cGy 72h UP   2
 10cGy 72h DN   0       0       0       0       0       0       0       0       0
100cGy  1h UP   9
100cGy  1h DN   4
100cGy  3h UP   0       0       0       0       0       0       0       0       0
100cGy  3h DN   0       0       0       0       0       0       0       0       0
100cGy 24h UP   0       0       0       0       0       0       0       0       0
100cGy 24h DN   0       0       0       0       0       0       0       0       0
100cGy 72h UP   1
100cGy 72h DN   1


