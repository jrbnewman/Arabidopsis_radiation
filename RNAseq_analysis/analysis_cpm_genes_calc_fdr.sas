ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* For the contrasts Mingqi wants to look at further,
   calc FDR and flag FDR 5%, 10%, 20%, and update corresponding up/down indicators */

%macro calcFDR(onflag1,onflag2,labelname,outname,fdrname);

data on_gene;
   set rs.arab_flag_gene_on_cpm_gt0;
   where &onflag1.=1 and &onflag2.=1 ;
   keep gene_id;
run;

data contrast_p;
  set rs.arab_gene_cntrs_constr_lcpm;
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

%calcFDR(flag_gene_on_01gy_1hr_cpm0,flag_gene_on_Mock_1hr_cpm0,"0.1gy-Mock: 1h",contrast1,01gy_v_Mock_1h);
%calcFDR(flag_gene_on_01gy_3hr_cpm0,flag_gene_on_Mock_3hr_cpm0,"0.1gy-Mock: 3h",contrast2,01gy_v_Mock_3h);
%calcFDR(flag_gene_on_01gy_24hr_cpm0,flag_gene_on_Mock_24hr_cpm0,"0.1gy-Mock: 24h",contrast3,01gy_v_Mock_24h);
%calcFDR(flag_gene_on_01gy_72hr_cpm0,flag_gene_on_Mock_72hr_cpm0,"0.1gy-Mock: 72h",contrast4,01gy_v_Mock_72h);

%calcFDR(flag_gene_on_1gy_1hr_cpm0,flag_gene_on_Mock_1hr_cpm0,"1gy-Mock: 1h",contrast5,1gy_v_Mock_1h);
%calcFDR(flag_gene_on_1gy_3hr_cpm0,flag_gene_on_Mock_3hr_cpm0,"1gy-Mock: 3h",contrast6,1gy_v_Mock_3h);
%calcFDR(flag_gene_on_1gy_24hr_cpm0,flag_gene_on_Mock_24hr_cpm0,"1gy-Mock: 24h",contrast7,1gy_v_Mock_24h);
%calcFDR(flag_gene_on_1gy_72hr_cpm0,flag_gene_on_Mock_72hr_cpm0,"1gy-Mock: 72h",contrast8,1gy_v_Mock_72h);

%macro calcFDRMain(labelname,outname,fdrname);

data main_p;
  set rs.arab_anova_gene_main_trt_tm_lcpm;
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

data rs.fdr_by_gene_log_cpm;
  set fdr_by_gene;
run;

