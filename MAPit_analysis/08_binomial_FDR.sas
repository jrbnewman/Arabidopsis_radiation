/* Binomial models -- FDR */

libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


/* DMRs */

data binomal_dmr;
  set wgbslocA.bin_dmr_anova_1 wgbslocA.bin_dmr_anova_2
      wgbslocA.bin_dmr_anova_3 wgbslocA.bin_dmr_anova_4
      wgbslocA.bin_dmr_anova_5;
  where effect="group";
  keep comparison site_type chr region_num ProbF;
run;

proc sort data=binomal_dmr;
  by comparison site_type ProbF;
proc multtest inpvalues(ProbF)=binomal_dmr fdr out=bin_dmr_fdr noprint;
  by comparison site_type;
run;

data bin_dmr_fdr2;
  set bin_dmr_fdr;
  if ProbF = . then flag_p05=.;
  else if ProbF < 0.05 then flag_p05=1;
  else  flag_p05=0;
  if fdr_P = . then flag_fdr05=.;
  else if fdr_P < 0.05 then flag_fdr05=1;
  else  flag_fdr05=0;
run;

proc freq data=bin_dmr_fdr2 noprint;
 tables comparison*site_type*flag_p05 / out=p_check;
 tables comparison*site_type*flag_fdr05 / out=fdr_check;
run;

data wgbslocA.binomial_dmr_fdr;
   set bin_dmr_fdr2;
   rename probf=binomial_P
          fdr_p=binomial_FDR_P
          flag_p05=flag_binomial_p05
          flag_fdr05=flag_binomial_fdr05;
run;




/* DARs */

data binomal_dar;
  set wgbslocA.bin_dar_contrasts2_all;
  where label="TRT_v_CTL" and probf ne .;
  keep comparison site_type chr region_num ProbF;
run;

proc sort data=binomal_dar;
  by comparison site_type ProbF;
proc multtest inpvalues(ProbF)=binomal_dar fdr out=bin_dar_fdr noprint;
  by comparison site_type;
run;

data bin_dar_fdr2;
  set bin_dar_fdr;
  if ProbF = . then flag_p05=.;
  else if ProbF < 0.05 then flag_p05=1;
  else  flag_p05=0;
  if fdr_P = . then flag_fdr05=.;
  else if fdr_P < 0.05 then flag_fdr05=1;
  else  flag_fdr05=0;
run;

proc freq data=bin_dar_fdr2 noprint;
 tables comparison*site_type*flag_p05 / out=p_check;
 tables comparison*site_type*flag_fdr05 / out=fdr_check;
run;

data wgbslocA.binomial_dar_fdr;
   set bin_dar_fdr2;
   rename probf=binomial_P
          fdr_p=binomial_FDR_P
          flag_p05=flag_binomial_p05
          flag_fdr05=flag_binomial_fdr05;
run;


