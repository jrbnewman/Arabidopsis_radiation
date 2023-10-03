libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Summarize DMRs, with the following columns:

Comparison
site_type
chr
region_start
region_stop
num_sites_in_region
mean_meth_CTL
mean_meth_TRT
mean_diff_TRT_CTL
num_FDR05
num_sites_10perc_diff
num_DMC_FDR05_10perc_diff
binomial_P
binomial_fdr_p
flag_binomial_P05
flag_binomial_FDR05

complete HOMER annotation

geneID
flag_promoter
flag_genebody
*/

/* Summarize DMCs and get mean methylation for region */

data dmc;
  set wgbslocA.results_by_dmc_w_meth;
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


data site2region;
   set wgbslocA.cytosine_to_meth_region_index;
   keep comparison region_num site_type chr start_pos stop_pos;
run;

proc sort data=dmc;
  by comparison site_type chr start_pos stop_pos;
proc sort data=site2region;
  by comparison site_type chr start_pos stop_pos;
run;

data dmc2region;
  merge dmc (in=in1) site2region (in=in2);
  by comparison site_type chr start_pos stop_pos;
  if in1 and in2;
run;

proc sort data=dmc2region;
  by comparison site_type chr region_num;
proc means data=dmc2region noprint;
  by comparison site_type chr region_num;
  var start_pos stop_pos methyl_TRT methyl_CTL methyl_diff
      flag_meth_diff_10perc flag_fdr05 flag_fdr05_10perc;
  output out=meth_summary_by_region 
         min(start_pos)=region_start
         max(stop_pos)=region_stop
         mean(methyl_TRT)=mean_methyl_TRT
         mean(methyl_CTL)=mean_methyl_CTL
         mean(methyl_diff)=mean_methyl_diff
         sum(flag_meth_diff_10perc)=num_sites_diff_10perc
         sum(flag_fdr05)=num_sites_FDR05
         sum(flag_fdr05_10perc)=num_sites_FDR05_diff_10perc;
run;


data meth_summary_by_region_ge2;
  set meth_summary_by_region;
  where _FREQ_ >= 2;
  rename _FREQ_ = num_sites_in_region;
  drop _TYPE_;
run;

/* Add binomial results */

data binomial;
   set wgbslocA.binomial_dmr_fdr;
  run;

proc sort data=binomial;
  by comparison site_type chr region_num;
proc sort data=meth_summary_by_region_ge2;
  by comparison site_type chr region_num;
run;

data region_meth_binom no_meth;
  merge meth_summary_by_region_ge2 (in=in1) binomial (in=in2);
  by comparison site_type chr region_num;
  if in1 then output region_meth_binom;
  else if in2 then output no_meth;
run;

/* Import HOMER annotations */

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_regions_HOMER_annotated.txt"
   out=annot dbms=tab replace;
   guessingrows=all;
run;

data annot2;
  set annot;
  length comparison $100.;
  length site_type $3.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'.'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'.'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'.'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'.')) + 0;
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;
 
proc sort data=annot2;
   by comparison site_type chr region_num;
proc sort data=region_meth_binom;
   by comparison site_type chr region_num;
run;


data region_meth_binom_annot;
   merge region_meth_binom (in=in1) annot2 (in=in2);
   by comparison site_type chr region_num;
   if in1 and in2;
run;

/* Make permanent */

data wgbslocA.results_by_DMR_annot;
  set region_meth_binom_annot;
run;


/* DAR */






/* Summarize DACs and get mean methylation for region */

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


data site2region;
   set wgbslocA.cytosine_to_acc_region_index;
   keep comparison region_num site_type chr start_pos stop_pos;
run;

proc sort data=dac;
  by comparison site_type chr start_pos stop_pos;
proc sort data=site2region;
  by comparison site_type chr start_pos stop_pos;
run;

data dac2region;
  merge dac (in=in1) site2region (in=in2);
  by comparison site_type chr start_pos stop_pos;
  if in1 and in2;
run;


proc sort data=dac2region;
  by comparison site_type chr region_num;
proc means data=dac2region noprint;
  by comparison site_type chr region_num;
  var start_pos stop_pos
  flag_CTL_fdr05 flag_FDR05_CTL_and_TRT
  flag_FDR05_CTL_or_TRT flag_TRT_fdr05
  flag_fdr05_20perc flag_meth_diff_20perc
  methyl_CTL_0U methyl_CTL_100U methyl_TRT_0U methyl_TRT_100U
  methyl_diff_CTL methyl_diff_TRT methyl_diff_TRT_CTL;
  output out=meth_summary_by_region 
         min(start_pos)=region_start
         max(stop_pos)=region_stop
         sum(flag_CTL_fdr05)=num_sites_FDR05_CTL
         sum(flag_FDR05_CTL_and_TRT)=num_sites_FDR05_CTL_and_TRT
         sum(flag_FDR05_CTL_or_TRT)=num_sites_FDR05_CTL_or_TRT
         sum(flag_TRT_fdr05)=num_sites_FDR05_TRT
         sum(flag_fdr05_20perc)=num_sites_FDR05_diff_20perc
         sum(flag_meth_diff_20perc)=num_sites_diff_20perc
         mean(methyl_CTL_0U)=mean_methyl_CTL_0U
         mean(methyl_CTL_100U)=mean_methyl_CTL_100U
         mean(methyl_TRT_0U)=mean_methyl_TRT_0U
         mean(methyl_TRT_100U)=mean_methyl_TRT_100U
         mean(methyl_diff_CTL)=mean_methyl_diff_CTL_100U_0U
         mean(methyl_diff_TRT)=mean_methyl_diff_TRT_100U_0U
         mean(methyl_diff_TRT_CTL)=mean_methyl_diff_TRT_CTL;
run;


data meth_summary_by_region_ge2;
  set meth_summary_by_region;
  where _FREQ_ >= 2;
  rename _FREQ_ = num_sites_in_region;
  drop _TYPE_;
run;

/* Add binomial results */

data binomial;
   set wgbslocA.binomial_dar_fdr;
  run;

proc sort data=binomial;
  by comparison site_type chr region_num;
proc sort data=meth_summary_by_region_ge2;
  by comparison site_type chr region_num;
run;

data region_meth_binom no_meth;
  merge meth_summary_by_region_ge2 (in=in1) binomial (in=in2);
  by comparison site_type chr region_num;
  if in1 then output region_meth_binom;
  else if in2 then output no_meth;
run;

/* Import HOMER annotations */

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_radiation_regions_HOMER_annotated.txt"
   out=annot dbms=tab replace;
   guessingrows=all;
run;

data annot2;
  set annot;
  length comparison $100.;
  length site_type $3.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'.'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'.'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'.'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'.')) + 0;
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;
 
proc sort data=annot2;
   by comparison site_type chr region_num;
proc sort data=region_meth_binom;
   by comparison site_type chr region_num;
run;


data region_meth_binom_annot;
   merge region_meth_binom (in=in1) annot2 (in=in2);
   by comparison site_type chr region_num;
   if in1 and in2;
run;

/* Make permanent */

data wgbslocA.results_by_DAR_annot;
  set region_meth_binom_annot;
run;

