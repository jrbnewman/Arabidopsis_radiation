libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';


ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DMRs_min_5_sites_for_HOMER_annotation.txt"
   out=dmr_annot dbms=tab replace;
   guessingrows=max;
run;


data dmr_annot2;
  set dmr_annot;
  length comparison $12.;
  length site_type $4.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  length feature $20.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'|'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'|'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'|'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'|')) + 0;
  feature = scan(annotation, 1, " ");
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;


data results_by_dmr;
   set wgbsA.results_by_dmr_5sites;
run;



proc sort data=dmr_annot2;
  by comparison site_type chr  region_num;
proc sort data=results_by_dmr;
  by comparison site_type chr  region_num;
run;

data dmr_w_annot;
  merge dmr_annot2 (in=in1) results_by_dmr (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;
data dmr_w_annot2;
  set dmr_w_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;


/* Count DMRs */

data dmg_Genes;
     set dmr_w_annot2;
     /* all DMRs by site type and comparison */
     if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output ;
     if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output ;
     if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output ;
     if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output ;

     if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output ;
     if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output ;
     if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output ;
     if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output ;
  keep geneid;
run;

proc sort data=dmg_genes nodup;
 by geneid;
run;


%macro concordance(siteType,treatment,units);

%if &siteType.=GC %then %do;

data meth_data_rep1;
   set arabMAP.methylation_data_gc;
   where rep=1 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   keep chr stop_pos perc_methyl_norm total_C;
   rename stop_pos=pos perc_methyl_norm=perc_methyl_&treatment._&units._1 total_C=total_C_1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_gc;
   where rep=2 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   keep chr stop_pos perc_methyl_norm total_C;
   rename stop_pos=pos perc_methyl_norm=perc_methyl_&treatment._&units._2 total_C=total_C_2;
run;

%end;

%else %do;

data meth_data_rep1;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=1 and site_type="&siteType." and treatment="&treatment."  ;
   keep chr stop_pos perc_methyl total_C;
   rename stop_pos=pos perc_methyl=perc_methyl_&treatment._&units._1 total_C=total_C_1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=2 and site_type="&siteType." and treatment="&treatment."  ;
   keep chr stop_pos perc_methyl total_C;
   rename stop_pos=pos perc_methyl=perc_methyl_&treatment._&units._2 total_C=total_C_2;
run;

%end;

proc sort data=meth_data_rep1;
  by chr pos;
proc sort data=meth_data_rep2;
  by chr pos;
run;

data meth_data_compare;
   merge meth_data_rep1 (in=in1) meth_Data_rep2 (in=in2);
  by chr pos;
  if in1 then flag_rep1_&treatment._&units.=1;
  else do; flag_rep1_&treatment._&units.=0; perc_methyl_&treatment._&units._1=0; total_c_1=0; end;
  if in2 then flag_rep2_&treatment._&units.=1;
  else do; flag_rep2_&treatment._&units.=0; perc_methyl_&treatment._&units._2=0; total_c_2=0; end;
run;

data meth_data_compare2;
  set meth_data_compare;
  max_C = total_C_1 + total_C_2;
  if max_C >= 10 then flag_keep=1; else flag_keep=0;
  run;

proc freq data=meth_data_compare2;
   where flag_keep=1;
   tables flag_rep1_&treatment._&units. * flag_rep2_&treatment._&units. ; 
run;

proc corr data=meth_data_compare2 pearson;
  where total_C_1 > 0 and total_C_2 > 0;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;

proc corr data=meth_data_compare2 pearson;
  where flag_keep=1 and flag_rep1_&treatment._&units. =1 and flag_rep2_&treatment._&units. = 1 ;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;


proc corr data=meth_data_compare2 pearson;
  where total_C_1 > 0 and total_C_2 > 0;
  var  total_C_1 total_C_2;
run;

proc corr data=meth_data_compare2 pearson;
  where flag_keep=1 and flag_rep1_&treatment._&units. =1 and flag_rep2_&treatment._&units. = 1 ;
  var  total_C_1 total_C_2;
run;

data export_Data;
   set meth_data_compare2;
   keep chr pos perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_percMeth_&siteType._&treatment._&units..csv"
     dbms=csv replace;
run;
data export_Data;
   set meth_data_compare2;
   keep chr pos total_C_1 total_C_2 ;
   rename total_C_1=total_C_&treatment._&units._1 total_C_2=total_C_&treatment._&units._2;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_totalC_&siteType._&treatment._&units..csv"
     dbms=csv replace;
run;

%mend;


%concordance(CG, 0Gy, 0U);
%concordance(CG, 01Gy, 0U);
%concordance(CG, 1Gy, 0U);

%concordance(CHG, 0Gy, 0U);
%concordance(CHG, 01Gy, 0U);
%concordance(CHG, 1Gy, 0U);

%concordance(CHH, 0Gy, 0U);
%concordance(CHH, 01Gy, 0U);
%concordance(CHH, 1Gy, 0U);

%concordance(GC, 0Gy, 0U);
%concordance(GC, 01Gy, 0U);
%concordance(GC, 1Gy, 0U);

%concordance(GC, 0Gy, 100U);
%concordance(GC, 01Gy, 100U);
%concordance(GC, 1Gy, 100U);



%macro concordance(siteType,treatment,units);

%if &siteType.=GC %then %do;

data meth_data_rep1;
   set arabMAP.methylation_data_gc;
   where rep=1 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads total_C; 
   rename  perc_methyl2=perc_methyl_&treatment._&units._1 total_C=total_C_1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_gc;
   where rep=2 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads total_C;
   rename perc_methyl2=perc_methyl_&treatment._&units._2 total_C=total_C_2;
run;

%end;

%else %do;

data meth_data_rep1;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=1 and site_type="&siteType." and treatment="&treatment."  ;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads total_C;
   rename perc_methyl2=perc_methyl_&treatment._&units._1  total_C=total_C_1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=2 and site_type="&siteType." and treatment="&treatment."  ;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads total_C;
   rename perc_methyl2=perc_methyl_&treatment._&units._2 total_C=total_C_2;
run;

%end;



proc sort data=meth_data_rep1;
  by chr chr_bin;
proc sort data=meth_data_rep2;
  by chr chr_bin;
run;

proc means data=meth_data_rep1 noprint;
 by chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl_&treatment._&units._1 total_C_1;
  output out=mean_methyl_by_bin_rep1 (drop=_TYPE_ _FREQ_) mean=;
run;

proc means data=meth_data_rep2 noprint;
 by chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl_&treatment._&units._2 total_C_2; 
  output out=mean_methyl_by_bin_rep2 (drop=_TYPE_ _FREQ_) mean=;
run;

proc sort data=mean_methyl_by_bin_rep1;
  by chr chr_bin;
proc sort data=mean_methyl_by_bin_rep2;
  by chr chr_bin;
run;


data meth_data_compare;
   merge mean_methyl_by_bin_rep1 (in=in1) mean_methyl_by_bin_rep2 (in=in2);
  by chr chr_bin;
  if in1 then flag_rep1_&treatment._&units.=1;
  else do; flag_rep1_&treatment._&units.=0; perc_methyl_&treatment._&units._1=0; total_c_1=0; end;
  if in2 then flag_rep2_&treatment._&units.=1;
  else do; flag_rep2_&treatment._&units.=0; perc_methyl_&treatment._&units._2=0; total_c_2=0; end;
run;

proc freq data=meth_data_compare;
   tables flag_rep1_&treatment._&units. * flag_rep2_&treatment._&units. ; 
run;

proc corr data=meth_data_compare pearson;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;


proc corr data=meth_data_compare pearson;
  where flag_rep1_&treatment._&units. =1 and flag_rep2_&treatment._&units. = 1 ;
  var  perc_methyl_&treatment._&units._1 perc_methyl_&treatment._&units._2;
run;
proc corr data=meth_data_compare pearson;
  var  total_c_1 total_c_2;
run;


proc corr data=meth_data_compare pearson;
  where flag_rep1_&treatment._&units. =1 and flag_rep2_&treatment._&units. = 1 ;
  var  total_c_1 total_c_2;
run;

data export_Data;
   set meth_data_compare2;
   keep chr pos perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_percMeth_&siteType._&treatment._&units._100X_100bp_binned.csv"
     dbms=csv replace;
run;
data export_Data;
   set meth_data_compare2;
   keep chr pos total_C_1 total_C_2 ;
   rename total_C_1=total_C_&treatment._&units._1 total_C_2=total_C_&treatment._&units._2;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_totalC_&siteType._&treatment._&units._100X_100bp_binned.csv"
     dbms=csv replace;
run;



%mend;


%concordance(CG, 0Gy, 0U);
%concordance(CG, 01Gy, 0U);
%concordance(CG, 1Gy, 0U);

%concordance(CHG, 0Gy, 0U);
%concordance(CHG, 01Gy, 0U);
%concordance(CHG, 1Gy, 0U);

%concordance(CHH, 0Gy, 0U);
%concordance(CHH, 01Gy, 0U);
%concordance(CHH, 1Gy, 0U);

%concordance(GC, 0Gy, 0U);
%concordance(GC, 01Gy, 0U);
%concordance(GC, 1Gy, 0U);

%concordance(GC, 0Gy, 100U);
%concordance(GC, 01Gy, 100U);
%concordance(GC, 1Gy, 100U);






%macro concordance(siteType,trt1,unit1,trt2,unit2);

%if &siteType.=GC %then %do;

data meth_data_trt1;
   set arabMAP.methylation_data_gc;
   where site_type="&siteType." and treatment="&trt1." and units="&unit1." and  flag_normalized=1;;
   keep site_type chr stop_pos perc_methyl_norm total_C;
   rename perc_methyl_norm=perc_methyl_&trt1._&unit1. total_C=total_C_&trt1._&unit1.;
run;


data meth_data_trt2;
   set arabMAP.methylation_data_gc;
   where site_type="&siteType." and treatment="&trt2." and units="&unit2." and  flag_normalized=1;;
   keep site_type chr stop_pos perc_methyl_norm total_C;
   rename perc_methyl_norm=perc_methyl_&trt1._&unit1. total_C=total_C_&trt2._&unit2.;
run;

%end;

%else %do;

data meth_data_trt1;
   set arabMAP.methylation_data_cg_chg_chh;
   where site_type="&siteType." and treatment="&trt1."  ;
   keep site_type chr stop_pos perc_methyl total_C;
   rename  perc_methyl=perc_methyl_&trt1._&unit1. total_C=total_C_&trt1._&unit1.;
run;


data meth_data_trt2;
   set arabMAP.methylation_data_cg_chg_chh;
   where site_type="&siteType." and treatment="&trt2."  ;
   keep site_type chr stop_pos perc_methyl total_C;
   rename perc_methyl=perc_methyl_&trt2._&unit2. total_C=total_C_&trt2._&unit2.;
run;

%end;



proc sort data=meth_Data_trt1;
  by site_type chr stop_pos;
proc sort data=meth_Data_trt2;
  by site_type chr stop_pos;
run;

proc means data=meth_data_trt1 noprint;
  by site_type chr stop_pos;
  var perc_methyl_&trt1._&unit1. total_C_&trt1._&unit1.;
  output out=meth_data_trt1_A mean(perc_methyl_&trt1._&unit1.)=perc_methyl_&trt1._&unit1.
                              sum(total_C_&trt1._&unit1.)=total_C_&trt1._&unit1.;
run;

proc means data=meth_data_trt2 noprint;
  by site_type chr stop_pos;
  var perc_methyl_&trt2._&unit2. total_C_&trt2._&unit2.;
  output out=meth_data_trt2_A mean(perc_methyl_&trt2._&unit2.)=perc_methyl_&trt2._&unit2.
                              sum(total_C_&trt2._&unit2.)=total_C_&trt2._&unit2.;
run;


proc sort data=meth_data_trt1_A;
  by chr stop_pos;
proc sort data=meth_data_trt2_A;
  by chr stop_pos;
run;

data meth_data_compare;
   merge meth_data_trt1_A (in=in1) meth_data_trt2_A (in=in2);
  by site_type chr stop_pos;
  if in1 then flag_&trt1._&unit1.=1;
  else do;  flag_&trt1._&unit1.=0; perc_methyl_&trt1._&unit1.=0; total_C_&trt1._&unit1.=0; end;
  if in2 then flag_&trt2._&unit2.=1;
  else do;  flag_&trt2._&unit2.=0; perc_methyl_&trt2._&unit2.=0; total_C_&trt2._&unit2.=0; end;
run;

data meth_data_compare2;
  set meth_data_compare;
  if total_C_&trt1._&unit1. >= 10 and  total_C_&trt2._&unit2. >= 10 then flag_keep=1; else flag_keep=0;
  run;

proc freq data=meth_data_compare2;
   where flag_keep=1;
   tables flag_&trt1._&unit1. * flag_&trt2._&unit2. ; 
run;

proc corr data=meth_data_compare2 pearson;
  where total_C_&trt1._&unit1. > 0 and total_C_&trt2._&unit2. > 0;
  var  perc_methyl_&trt1._&unit1. perc_methyl_&trt2._&unit2.;
run;

proc corr data=meth_data_compare2 pearson;
  where flag_keep=1 and flag_&trt1._&unit1. =1 and flag_&trt2._&unit2. = 1 ;
  var  perc_methyl_&trt1._&unit1. perc_methyl_&trt2._&unit2.;
run;

proc corr data=meth_data_compare2 pearson;
  where total_C_&trt1._&unit1. > 0 and total_C_&trt2._&unit2. > 0;
  var  total_C_&trt1._&unit1. total_C_&trt2._&unit2.;
run;

proc corr data=meth_data_compare2 pearson;
  where flag_keep=1 and flag_&trt1._&unit1. =1 and flag_&trt2._&unit2. = 1 ;
  var  total_C_&trt1._&unit1. total_C_&trt2._&unit2.;
run;

data export_Data;
   set meth_data_compare2;
   keep chr stop_pos perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_percMeth_&siteType._&trt1._&unit1._&trt2._&unit2..csv"
     dbms=csv replace;
run;

data export_Data2;
   set meth_data_compare2;
   keep chr stop_pos total_C_: ;
run;

proc export data=export_data2
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_totalC_&siteType._&trt1._&unit1._&trt2._&unit2..csv"
     dbms=csv replace;
run;

%mend;


%concordance(CG, 0Gy, 0U, 01Gy, 0U);
%concordance(CG, 0Gy, 0U, 1Gy, 0U);
%concordance(CHG, 0Gy, 0U, 01Gy, 0U);
%concordance(CHG, 0Gy, 0U, 1Gy, 0U);
%concordance(CHH, 0Gy, 0U, 01Gy, 0U);
%concordance(CHH, 0Gy, 0U, 1Gy, 0U);
%concordance(GC, 0Gy, 0U, 01Gy, 0U);
%concordance(GC, 0Gy, 0U, 1Gy, 0U);
%concordance(GC, 0Gy, 100U, 01Gy, 100U);
%concordance(GC, 0Gy, 100U, 1Gy, 100U);
%concordance(GC, 0Gy, 0U, 0Gy, 100U);
%concordance(GC, 01Gy, 0U, 01Gy, 100U);
%concordance(GC, 1Gy, 0U, 1Gy, 100U);






/* 100X sites only and bin into 100bp windows */



%macro concordance(siteType,treatment,units);





%if &siteType.=GC %then %do;

data meth_data_trt1;
   set arabMAP.methylation_data_gc;
   where site_type="&siteType." and treatment="&trt1." and units="&unit1." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads_&trt1._&unit1.=1;
   else flag_lt100_reads_&trt1._&unit1.=0;
   keep site_type chr stop_pos perc_methyl2 total_C flag_lt100_reads_&trt1._&unit1. chr_bin;
   rename  perc_methyl2=perc_methyl_&trt1._&unit1. total_C=total_C_&trt1._&unit1.;
run;


data meth_data_trt2;
   set arabMAP.methylation_data_gc;
   where site_type="&siteType." and treatment="&trt2." and units="&unit2." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads_&trt2._&unit2.=1;
   else flag_lt100_reads_&trt2._&unit2.=0;
   keep site_type chr stop_pos perc_methyl2 total_C flag_lt100_reads_&trt2._&unit2. chr_bin;
   rename perc_methyl_norm=perc_methyl_&trt2._&unit2. total_C=total_C_&trt2._&unit2.;
run;

%end;

%else %do;

data meth_data_trt1;
   set arabMAP.methylation_data_cg_chg_chh;
   where site_type="&siteType." and treatment="&trt1."  ;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads_&trt1._&unit1.=1;
   else flag_lt100_reads_&trt1._&unit1.=0;
   keep site_type chr stop_pos perc_methyl2 total_C flag_lt100_reads_&trt1._&unit1. chr_bin;
   rename  perc_methyl2=perc_methyl_&trt1._&unit1. total_C=total_C_&trt1._&unit1.;
run;


data meth_data_trt2;
   set arabMAP.methylation_data_cg_chg_chh;
   where site_type="&siteType." and treatment="&trt2."  ;
    chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads_&trt2._&unit2.=1;
   else flag_lt100_reads_&trt2._&unit2.=0;
   keep site_type chr stop_pos perc_methyl2 total_C flag_lt100_reads_&trt2._&unit2. chr_bin;
   rename perc_methyl2=perc_methyl_&trt2._&unit2. total_C=total_C_&trt2._&unit2.;
run;

%end;



proc sort data=meth_Data_trt1;
  by site_type chr stop_pos chr_bin flag_lt100_reads_&trt1._&unit1.;
proc sort data=meth_Data_trt2;
  by site_type chr stop_pos chr_bin flag_lt100_reads_&trt2._&unit2.;
run;

proc means data=meth_data_trt1 noprint;
  by site_type chr stop_pos chr_bin ;
  var perc_methyl_&trt1._&unit1. total_C_&trt1._&unit1. flag_lt100_reads_&trt1._&unit1.;
  output out=meth_data_trt1_A mean(perc_methyl_&trt1._&unit1.)=perc_methyl_&trt1._&unit1.
                              sum(total_C_&trt1._&unit1.)=total_C_&trt1._&unit1.
                              max(flag_lt100_reads_&trt1._&unit1.)=flag_lt100_reads_&trt1._&unit1.;
run;

proc means data=meth_data_trt2 noprint;
  by site_type chr stop_pos chr_bin ;
  var perc_methyl_&trt2._&unit2. total_C_&trt2._&unit2. flag_lt100_reads_&trt2._&unit2.;
  output out=meth_data_trt2_A mean(perc_methyl_&trt2._&unit2.)=perc_methyl_&trt2._&unit2.
                              sum(total_C_&trt2._&unit2.)=total_C_&trt2._&unit2.
                              max(flag_lt100_reads_&trt2._&unit2.)=flag_lt100_reads_&trt2._&unit2.;
run;


data meth_data_trt1_A2;
 set meth_data_trt1_A;
 if total_C_&trt1._&unit1. >= 100 then flag_lt100_reads_&trt1._&unit1.=1;
 else  flag_lt100_reads_&trt1._&unit1.=0;
run;

data meth_data_trt2_A2;
 set meth_data_trt2_A;
 if total_C_&trt2._&unit2. >= 100 then flag_lt100_reads_&trt2._&unit2.=1;
 else  flag_lt100_reads_&trt2._&unit2.=0;
run;


proc sort data=meth_data_trt1_A2;
  by chr chr_bin;
proc sort data=meth_data_trt2_A2;
  by chr chr_bin;
run;

proc means data=meth_data_trt1_A2 noprint;
 by chr chr_bin;
  where flag_lt100_reads_&trt1._&unit1.=0;
  var perc_methyl_&trt1._&unit1. total_C_&trt1._&unit1.;
  output out=mean_methyl_by_bin_&trt1._&unit1. (drop=_TYPE_ _FREQ_) mean=;
run;

proc means data=meth_data_trt2_A2 noprint;
 by chr chr_bin;
  where flag_lt100_reads_&trt2._&unit2.=0;
  var perc_methyl_&trt2._&unit2. total_C_&trt2._&unit2.;
  output out=mean_methyl_by_bin_&trt2._&unit2. (drop=_TYPE_ _FREQ_) mean=;
run;

proc sort data=mean_methyl_by_bin_&trt1._&unit1.;
  by chr chr_bin;
proc sort data=mean_methyl_by_bin_&trt2._&unit2.;
  by chr chr_bin;
run;


data meth_data_compare;
   merge mean_methyl_by_bin_&trt1._&unit1. (in=in1) mean_methyl_by_bin_&trt2._&unit2. (in=in2);
  by chr chr_bin;
  if in1 then flag_&trt1._&unit1.=1;
  else do; flag_&trt1._&unit1.=0;
           total_C_&trt1._&unit1. = 0;
           perc_methyl_&trt1._&unit1.=0;
           end;
  if in2 then flag_&trt2._&unit2.=1;
  else do; flag_&trt2._&unit2.=0;
           total_C_&trt2._&unit2. = 0;
           perc_methyl_&trt2._&unit2.=0;
           end;
run;

proc freq data=meth_data_compare;
   tables  flag_&trt1._&unit1. * flag_&trt2._&unit2.; 
run;

proc corr data=meth_data_compare pearson;
  var  perc_methyl_&trt1._&unit1. perc_methyl_&trt2._&unit2.;
run;

proc corr data=meth_data_compare pearson;
  where flag_&trt1._&unit1. =1 and flag_&trt2._&unit2. = 1 ;
  var  perc_methyl_&trt1._&unit1. perc_methyl_&trt2._&unit2.;
run;




proc corr data=meth_data_compare pearson;
  var  total_C_&trt1._&unit1. total_C_&trt2._&unit2.;
run;

proc corr data=meth_data_compare pearson;
  where flag_&trt1._&unit1. =1 and flag_&trt2._&unit2. = 1 ;
  var  total_C_&trt1._&unit1. total_C_&trt2._&unit2.;
run;





data export_Data;
   set meth_data_compare;
   keep chr chr_bin perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_percMeth_&siteType._&trt1._&unit1._&trt2._&unit2._100X_100bp_binned.csv"
     dbms=csv replace;
run;

data export_Data;
   set meth_data_compare;
   keep chr  chr_bin total_C_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_totalC_&siteType._&trt1._&unit1._&trt2._&unit2._100X_100bp_binned.csv"
     dbms=csv replace;
run;

%mend;






%concordance(CG, 0Gy, 0U, 01Gy, 0U);
%concordance(CG, 0Gy, 0U, 1Gy, 0U);
%concordance(CHG, 0Gy, 0U, 01Gy, 0U);
%concordance(CHG, 0Gy, 0U, 1Gy, 0U);
%concordance(CHH, 0Gy, 0U, 01Gy, 0U);
%concordance(CHH, 0Gy, 0U, 1Gy, 0U);
%concordance(GC, 0Gy, 0U, 01Gy, 0U);
%concordance(GC, 0Gy, 0U, 1Gy, 0U);
%concordance(GC, 0Gy, 100U, 01Gy, 100U);
%concordance(GC, 0Gy, 100U, 1Gy, 100U);
%concordance(GC, 0Gy, 0U, 0Gy, 100U);
%concordance(GC, 01Gy, 0U, 01Gy, 100U);
%concordance(GC, 1Gy, 0U, 1Gy, 100U);




