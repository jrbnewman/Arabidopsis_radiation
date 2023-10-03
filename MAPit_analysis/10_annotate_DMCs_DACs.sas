libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


/* Annotate DMCs */

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_sites_1_HOMER_annotated.txt"
    out=annot_1 dbms=tab replace;
    guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_sites_2_HOMER_annotated.txt"
    out=annot_2 dbms=tab replace;
    guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_sites_3_HOMER_annotated.txt"
    out=annot_3 dbms=tab replace;
    guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_sites_4_HOMER_annotated.txt"
    out=annot_4 dbms=tab replace;
    guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_sites_5_HOMER_annotated.txt"
    out=annot_5 dbms=tab replace;
    guessingrows=max;
run;

data annot_all;
  length gene_alias $358.;
  length gene_description $133.;
  set annot_1 annot_2 annot_3 annot_4 annot_5;
  length chrom $3.;
  length geneID $100.;
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'_'));
  start_pos=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'_'))+0;
  stop_pos=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'_'))+0;
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr start end;
   rename chrom=chr;
run;

data dmc;
  set wgbslocA.results_By_dmc_W_meth;
  drop XPL_FISH XPR_FISH;
  if FET_P = . then flag_p05=.;
  else if FET_P < 0.05 then flag_p05=1;
  else flag_p05=0;

  if FET_FDR_P = . then flag_fdr05=.;
  else if FET_FDR_P < 0.05 then flag_fdr05=1;
  else flag_fdr05=0;

  if flag_fdr05=1 and flag_meth_diff_10perc=1 then flag_fdr05_10perc=1;
  else flag_fdr05_10perc=0;
run;


proc sort data=dmc;
  by chr start_pos stop_pos;
proc sort data=annot_all;
  by chr start_pos stop_pos;
run;

data dmc_annot;
  merge dmc (in=in1) annot_all (in=in2);
  by chr start_pos stop_pos;
  if in1 and in2;
run;


data wgbslocA.results_by_dmc_annot;
   set dmc_annot;
run;


/* Annotate DACs */

data dac;
  set wgbslocA.results_by_dac_w_meth;
  if FET_FDR_P_TRT_100U_0U=. then flag_TRT_fdr05=.;
  else if FET_FDR_P_TRT_100U_0U < 0.05 then flag_TRT_fdr05=1;
  else flag_TRT_fdr05=0;

  if FET_FDR_P_CTL_100U_0U=. then flag_CTL_fdr05=.;
  else if FET_FDR_P_CTL_100U_0U < 0.05 then flag_CTL_fdr05=1;
  else flag_CTL_fdr05=0;

   if flag_TRT_fdr05=1 and flag_CTL_fdr05=1
     then flag_FDR05_CTL_and_TRT=1;
     else flag_FDR05_CTL_and_TRT=0;
  if flag_FDR05_CTL_or_TRT=1 and flag_meth_diff_20perc=1 then flag_fdr05_20perc=1;
  else flag_fdr05_20perc=0;
run;


proc sort data=dac;
  by chr start_pos stop_pos;
proc sort data=annot_all;
  by chr start_pos stop_pos;
run;

data dac_annot;
  merge dac (in=in1) annot_all (in=in2);
  by chr start_pos stop_pos;
  if in1 and in2;
run;


data wgbslocA.results_by_dac_annot;
   set dac_annot;
run;


 
