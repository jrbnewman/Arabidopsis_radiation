/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

/* Import and set up genes */



/* Get all regions */



proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DMRs_min_5_sites_for_HOMER_annotation.txt"
   out=dmr_annot dbms=tab replace;
   guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER_annotation.txt"
   out=dar_annot dbms=tab replace;
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
   geneID=compress(tranwrd(upcase(scan(scan(annotation,2,"("),1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;

data dar_annot2;
  set dar_annot;
  length comparison $12.;
  length site_type $4.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'|'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'|'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'|'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'|')) + 0;
  feature = scan(annotation, 1, " ");
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(scan(annotation,2,"("),1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;


data results_by_dmr;
   set wgbsA.results_by_dmr_5sites;
run;

data results_by_dar;
   set wgbsA.results_by_dar_5sites;
run;

proc freq data=results_by_dar;
  tables comparison*flag_fdr05;
run;


proc sort data=dmr_annot2;
  by comparison site_type chr  region_num;
proc sort data=dar_annot2;
  by comparison site_type chr  region_num;
proc sort data=results_by_dmr;
  by comparison site_type chr  region_num;
proc sort data=results_by_dar;
  by comparison site_type chr  region_num;
run;

data dmr_w_annot;
  merge dmr_annot2 (in=in1) results_by_dmr (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

data dar_w_annot;
  merge dar_annot2 (in=in1) results_by_dar (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


data dmr_w_annot2;
  set dmr_w_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;


data dar_w_annot2;
  set dar_w_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;

/* Count DMRs */

data all_01_up all_01_dn all_1_up all_1_dn
     cg_01_up cg_01_dn  cg_1_up cg_1_dn
     chg_01_up chg_01_dn  chg_1_up chg_1_dn
     chh_01_up chh_01_dn  chh_1_up chh_1_dn;
     set dmr_w_annot2;
     /* all DMRs by site type and comparison */
     if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output all_01_up;
     if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output all_01_dn;
     if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output all_1_up;
     if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output all_1_dn;

     if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output all_01_up;
     if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output all_01_dn;
     if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output all_1_up;
     if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output all_1_dn;
     if site_type="CG" then do;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output cg_01_up;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output cg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output cg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output cg_1_dn;

         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output cg_01_up;
         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output cg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output cg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output cg_1_dn;
     end;
     if site_type="CHG" then do;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chg_01_up;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chg_1_dn;

         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output chg_01_up;
         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output chg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output chg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output chg_1_dn;
     end;
     if site_type="CHH" then do;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chh_01_up;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chh_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chh_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chh_1_dn;

         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output chh_01_up;
         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output chh_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output chh_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output chh_1_dn;
     end;
run;

proc sort data=all_01_up nodup; by _all_;
proc sort data=all_01_dn nodup; by _all_;
proc sort data=all_1_up nodup; by _all_;
proc sort data=all_1_dn nodup; by _all_;
proc sort data=cg_01_up nodup; by _all_;
proc sort data=cg_01_dn nodup; by _all_;
proc sort data=cg_1_up nodup; by _all_;
proc sort data=cg_1_dn nodup; by _all_;
proc sort data=chg_01_up nodup; by _all_;
proc sort data=chg_01_dn nodup; by _all_;
proc sort data=chg_1_up nodup; by _all_;
proc sort data=chg_1_dn nodup; by _all_;
proc sort data=chh_01_up nodup; by _all_;
proc sort data=chh_01_dn nodup; by _all_;
proc sort data=chh_1_up nodup; by _all_;
proc sort data=chh_1_dn nodup; by _all_;
run;


/* Count DARs */

data gc_01_up gc_01_dn gc_1_up gc_1_dn;
     set dar_w_annot2;
     /* all DMRs by site type and comparison */
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output gc_01_up;
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output gc_01_dn;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output gc_1_up;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output gc_1_dn;

     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output gc_01_up;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output gc_01_dn;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output gc_1_up;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output gc_1_dn;
run;

proc sort data=gc_01_up nodup; by _all_;
proc sort data=gc_01_dn nodup; by _all_;
proc sort data=gc_1_up nodup; by _all_;
proc sort data=gc_1_dn nodup; by _all_;
run;


/* Import DAR-only regions for me to subset */

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_cg_vs_gc_01/Intervene_results/sets/01_all_GC_DMR_01.bed"
out=cg01_dar2remove dbms=tab replace;
getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chg_vs_gc_01/Intervene_results/sets/01_all_GC_DMR_01.bed"
out=chg01_dar2remove dbms=tab replace;
getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chh_vs_gc_01/Intervene_results/sets/01_all_GC_DMR_01.bed"
out=chh01_dar2remove dbms=tab replace;
getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_cg_vs_gc_1/Intervene_results/sets/01_all_GC_DMR_1.bed"
out=cg1_dar2remove dbms=tab replace;
getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chg_vs_gc_1/Intervene_results/sets/01_all_GC_DMR_1.bed"
out=chg1_dar2remove dbms=tab replace;
getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chh_vs_gc_1/Intervene_results/sets/01_all_GC_DMR_1.bed"
out=chh1_dar2remove dbms=tab replace;
getnames=no;
run;



proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_cg_vs_gc_01/Intervene_results/sets/10_all_CG_DMR_01.bed"
out=cg01_dmr2remove dbms=tab replace;
getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chg_vs_gc_01/Intervene_results/sets/10_all_CHG_DMR_01.bed"
out=chg01_dmr2remove dbms=tab replace;
getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chh_vs_gc_01/Intervene_results/sets/10_all_CHH_DMR_01.bed"
out=chh01_dmr2remove dbms=tab replace;
getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_cg_vs_gc_1/Intervene_results/sets/10_all_CG_DMR_1.bed"
out=cg1_dmr2remove dbms=tab replace;
getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chg_vs_gc_1/Intervene_results/sets/10_all_CHG_DMR_1.bed"
out=chg1_dmr2remove dbms=tab replace;
getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/intervene_chh_vs_gc_1/Intervene_results/sets/10_all_CHH_DMR_1.bed"
out=chh1_dmr2remove dbms=tab replace;
getnames=no;
run;








data cg01_dar2remove2;
  set cg01_dar2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data cg1_dar2remove2;
  set cg1_dar2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data chg01_dar2remove2;
  set chg01_dar2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data chg1_dar2remove2;
  set chg1_dar2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data chh01_dar2remove2;
  set chh01_dar2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;
data chh1_dar2remove2;
  set chh1_dar2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;


proc sort data=cg01_dar2remove2;
  by chr region_start region_stop;
proc sort data=chg01_dar2remove2;
  by chr region_start region_stop;
proc sort data=chh01_dar2remove2;
  by chr region_start region_stop;
proc sort data=gc_01_up;
  by chr region_start region_stop;
proc sort data=gc_01_dn;
  by chr region_start region_stop;
run;




data cg01_dmr2remove2;
  set cg01_dmr2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data cg1_dmr2remove2;
  set cg1_dmr2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data chg01_dmr2remove2;
  set chg01_dmr2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data chg1_dmr2remove2;
  set chg1_dmr2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;

data chh01_dmr2remove2;
  set chh01_dmr2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;
data chh1_dmr2remove2;
  set chh1_dmr2remove;
  length chr $4.;
  chr=compress(scan(VAR4,1,"_"));
  region_start=compress(scan(VAR4,2,"_")) + 0;
  region_stop=compress(scan(VAR4,3,"_")) + 0;
run;


proc sort data=cg01_dmr2remove2;
  by chr region_start region_stop;
proc sort data=chg01_dmr2remove2;
  by chr region_start region_stop;
proc sort data=chh01_dmr2remove2;
  by chr region_start region_stop;
proc sort data=cg1_dmr2remove2;
  by chr region_start region_stop;
proc sort data=chg1_dmr2remove2;
  by chr region_start region_stop;
proc sort data=chh1_dmr2remove2;
  by chr region_start region_stop;
proc sort data=cg_01_up;
  by chr region_start region_stop;
proc sort data=cg_1_dn;
  by chr region_start region_stop;
proc sort data=chg_01_up;
  by chr region_start region_stop;
proc sort data=chg_1_dn;
  by chr region_start region_stop;
proc sort data=chh_01_up;
  by chr region_start region_stop;
proc sort data=chh_1_dn;
  by chr region_start region_stop;

proc sort data=cg_1_up;
  by chr region_start region_stop;
proc sort data=cg_01_dn;
  by chr region_start region_stop;
proc sort data=chg_1_up;
  by chr region_start region_stop;
proc sort data=chg_01_dn;
  by chr region_start region_stop;
proc sort data=chh_1_up;
  by chr region_start region_stop;
proc sort data=chh_01_dn;
  by chr region_start region_stop;
run;


data cg2gc_01_up ;
  merge cg_01_up (in=in1) cg01_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data cg2gc_01_dn ;
  merge cg_01_dn (in=in1) cg01_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;
data cg2gc_1_up ;
  merge cg_1_up (in=in1) cg1_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data cg2gc_1_dn ;
  merge cg_1_dn (in=in1) cg1_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;




data chg2gc_01_up ;
  merge chg_01_up (in=in1) chg01_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data chg2gc_01_dn ;
  merge chg_01_dn (in=in1) chg01_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data chg2gc_1_up ;
  merge chg_1_up (in=in1) chg1_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data chg2gc_1_dn ;
  merge chg_1_dn (in=in1) chg1_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;






data chh2gc_01_up ;
  merge chh_01_up (in=in1) chh01_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data chh2gc_01_dn ;
  merge chh_01_dn (in=in1) chh01_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;
data chh2gc_1_up ;
  merge chh_1_up (in=in1) chh1_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;

data chh2gc_1_dn ;
  merge chh_1_dn (in=in1) chh1_dmr2remove2 (in=in2);
  by chr region_start region_stop;
     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
  if feature = "Intergenic" then flag_intergenic=1; else flag_intergenic=0;
  if in1 and not in2 then output;
  keep chr region_start region_stop flag_intergenic dar_center plot_start plot_stop;
run;





proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/Arabidopsis_thaliana.TAIR10.37.gtf"
  out=gtf dbms=tab replace;
  guessingrows=max;
  getnames=no;
run;


data gene;
  set gtf;
  where VAR3 = "gene";
  length gene_id $15.;
  gene_id=compress(tranwrd(tranwrd(scan(VAR9,2," "), ";", ""), '"', ''));
  keep gene_id VAR1 VAR4 VAR5 VAR7; 
  rename VAR7=strand VAR1=chr VAR4=gene_start VAR5=gene_stop;
run;




data gc_01_up_gn gc_01_dn_gn gc_1_up_gn gc_1_dn_gn;
     set dar_w_annot2;
     /* all DMRs by site type and comparison */
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output gc_01_up_gn;
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output gc_01_dn_gn;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output gc_1_up_gn;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output gc_1_dn_gn;

     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output gc_01_up_gn;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output gc_01_dn_gn;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output gc_1_up_gn;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output gc_1_dn_gn;
  keep geneID;
  rename geneID=gene_id;
run;

proc sort data=gc_01_up_gn nodup; by  gene_id;
proc sort data=gc_01_dn_gn nodup; by  gene_id;
proc sort data=gc_1_up_gn nodup; by  gene_id;
proc sort data=gc_1_dn_gn nodup; by  gene_id;
run;

data da_Gene_tabs;
  merge gc_01_up_gn (in=in1) gc_01_dn_gn (in=in2) gc_1_up_gn (in=in3) gc_1_dn_gn (in=in4);
  by gene_id;
  if in1 then flag_gc_01_up=1; else flag_gc_01_up=0;
  if in2 then flag_gc_01_dn=1; else flag_gc_01_dn=0;
  if in3 then flag_gc_1_up=1; else flag_gc_1_up=0;
  if in4 then flag_gc_1_dn=1; else flag_gc_1_dn=0;
run;

proc freq data=da_gene_tabs ;
 tables flag_gc_01_up*flag_gc_01_dn;
 tables flag_gc_1_up*flag_gc_1_dn;
run;

proc sort data=da_gene_tabs;
  by gene_id;
proc sort data=gene;
  by gene_id;
run;

data gene2;
  merge gene (in=in1) da_gene_tabs (in=in2);
  by gene_id;
  if in1 and in2;
run;




/* prep methylation and accessibility counts */


data meth_data;
  set arabMAP.methylation_data_gc;
  keep site_type chr stop_pos treatment units rep total_C methyl_C perc_methyl_norm;
run;

proc sort data=meth_data;
  by site_type chr stop_pos treatment units  rep ;
proc means data=meth_data noprint;
  by site_type chr stop_pos treatment  units  ;
  var total_C methyl_C perc_methyl_norm;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl_norm)=perc_methyl;
run;


proc transpose data=meth_data2 out=meth_sbys10;
  where total_C >= 10;
  by site_type chr stop_pos;
  id treatment units ;
  var perc_methyl;
run;

data meth_sbys10_2;
  set meth_sbys10;
  if _01Gy100U = . or _01Gy0U = . then _01Gy=.; else _01Gy = _01Gy100U - _01Gy0U;
  if _1Gy100U = . or _1Gy0U = . then _1Gy=.; else _1Gy = _1Gy100U - _1Gy0U;
  if _0Gy100U = . or _0Gy0U = . then _0Gy=.; else _0Gy = _0Gy100U - _0Gy0U;


  if _01Gy ne . and _0Gy ne . then _01Gy_common=_01Gy; else _01Gy_common=.;
  if _1Gy ne . and _0Gy ne . then _1Gy_common=_1Gy; else _1Gy_common=.;
  if (_1Gy ne . or _01Gy ne .) and _0Gy ne . then _0Gy_common=_0Gy; else _0Gy_common=.;

  if _01Gy_common ne . and _1Gy_common ne .  then do;
    _01Gy_common_all = _01Gy_common;
    _1Gy_common_all = _1Gy_common;
    _0Gy_common_all = _0Gy_common;
    end;
  else do;
   _01Gy_common_all = .;
   _1Gy_common_all = .;
   _0Gy_common_all = . ;
    end;
  rename stop_pos=pos;
run;


proc sort data=gene2;
  by chr gene_id;
run;

data gene_2;
  set gene2;
  by chr gene_id;
  gene_length=gene_stop-gene_start;
  do pos = gene_start to gene_stop;
  output;
  end;
run;

proc sort data=gene_2;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data meth_sbys10_gene;
  merge meth_sbys10_2 (in=in1) gene_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data meth_sbys10_gene2;
  set meth_sbys10_gene;
  if strand="+" then pos_bin=int((pos-gene_start)/gene_length * 1000) + 1000 ;
  else if strand="-" then pos_bin=int((gene_stop-pos)/gene_length * 1000) + 1000 ;
run;

data promoter;
  set gene2;
  if strand="+" then do;
    promoter_start=gene_start-1000;
    promoter_stop=gene_start-1;
    end;
  else do;
    promoter_start=gene_stop+1000;
    promoter_stop=gene_stop+1;
    end;
run;

data downstream;
  set gene2;
  if strand="-" then do;
    downstream_start=gene_start-1;
    downstream_stop=gene_start-1000;
    end;
  else do;
    downstream_start=gene_stop+1;
    downstream_stop=gene_stop+1000;
    end;
run;

proc sort data=promoter;
  by chr gene_id ;
proc sort data=downstream;
  by chr gene_id ;
run;


data promoter_2;
  set promoter;
  by chr gene_id;
  if strand="+" then do pos = promoter_start to promoter_stop;
  output;
  end;
  else if strand="-" then do pos = promoter_stop to promoter_start;
  output;
  end;
run;


data downstream_2;
  set downstream;
  by chr gene_id;
  if strand="+" then do pos = downstream_start to downstream_stop;
  output;
  end;
  else if strand="-" then do pos = downstream_stop to downstream_start;
  output;
  end;
run;



proc sort data=promoter_2;
  by chr pos;
proc sort data=downstream_2;
  by chr pos;
run;

data meth_sbys10_promoter;
  merge meth_sbys10_2 (in=in1) promoter_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_downstream;
  merge meth_sbys10_2 (in=in1) downstream_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_promoter2;
  set meth_sbys10_promoter;
  if strand="+" then pos_bin=(pos-promoter_start);
  else if strand="-" then pos_bin=(promoter_start-pos) ;
run;

data meth_sbys10_downstream2;
  set meth_sbys10_downstream;
  if strand="+" then pos_bin=(pos-downstream_start) + 2000;
  else if strand="-" then pos_bin=(downstream_start-pos) + 2000;
run;



data meth_sbys10_genestack;
   set meth_sbys10_gene2 (in=in1) meth_sbys10_promoter2 (in=in2) meth_sbys10_downstream2 (in=in3);
run;

proc sort data=meth_sbys10_genestack;
  by gene_id chr pos_bin;
run;

data meth_data2;
  set meth_sbys10_genestack;
  grouped_pos=int(pos_bin/10)*10;
run;


proc sort data=meth_data2;
  by pos_bin grouped_pos;
run;


proc means data=meth_data2 noprint;
  by  pos_bin grouped_pos  ;
  var _01Gy _1Gy _0Gy   ;
  output out=meth_data3
  mean=;
run;


proc sort data=meth_data3 ;
  by grouped_pos;
run;


proc means data=meth_data3 noprint;
  by grouped_pos  ;
  var  _01Gy _1Gy _0Gy    
       ;
  output out=meth_data4
  mean=;
run;


data meth_data_export_acc;
  set meth_data4;
  drop _TYPE_ _FREQ_ ;
  rename grouped_pos=pos;
run;








/* prep methylation counts */


data meth_data;
  set arabMAP.methylation_data_cg_chg_chh;
  keep site_type chr stop_pos treatment units rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by site_type chr stop_pos treatment   rep ;
proc means data=meth_data noprint;
  by site_type chr stop_pos treatment    ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;


proc transpose data=meth_data2 out=meth_sbys10;
  where total_C >= 10;
  by site_type chr stop_pos;
  id treatment  ;
  var perc_methyl;
run;

data meth_sbys10_2;
  set meth_sbys10;


  if _01Gy ne . and _0Gy ne . then _01Gy_common=_01Gy; else _01Gy_common=.;
  if _1Gy ne . and _0Gy ne . then _1Gy_common=_1Gy; else _1Gy_common=.;
  if (_1Gy ne . or _01Gy ne .) and _0Gy ne . then _0Gy_common=_0Gy; else _0Gy_common=.;

  if _01Gy_common ne . and _1Gy_common ne .  then do;
    _01Gy_common_all = _01Gy_common;
    _1Gy_common_all = _1Gy_common;
    _0Gy_common_all = _0Gy_common;
    end;
  else do;
   _01Gy_common_all = .;
   _1Gy_common_all = .;
   _0Gy_common_all = . ;
    end;
  rename stop_pos=pos;
run;


proc sort data=gene2;
  by chr gene_id;
run;

data gene_2;
  set gene2;
  by chr gene_id;
  gene_length=gene_stop-gene_start;
  do pos = gene_start to gene_stop;
  output;
  end;
run;

proc sort data=gene_2;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data meth_sbys10_gene;
  merge meth_sbys10_2 (in=in1) gene_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data meth_sbys10_gene2;
  set meth_sbys10_gene;
  if strand="+" then pos_bin=int((pos-gene_start)/gene_length * 1000) + 1000 ;
  else if strand="-" then pos_bin=int((gene_stop-pos)/gene_length * 1000) + 1000 ;
run;

data promoter;
  set gene2;
  if strand="+" then do;
    promoter_start=gene_start-1000;
    promoter_stop=gene_start-1;
    end;
  else do;
    promoter_start=gene_stop+1000;
    promoter_stop=gene_stop+1;
    end;
run;

data downstream;
  set gene2;
  if strand="-" then do;
    downstream_start=gene_start-1;
    downstream_stop=gene_start-1000;
    end;
  else do;
    downstream_start=gene_stop+1;
    downstream_stop=gene_stop+1000;
    end;
run;

proc sort data=promoter;
  by chr gene_id ;
proc sort data=downstream;
  by chr gene_id ;
run;


data promoter_2;
  set promoter;
  by chr gene_id;
  if strand="+" then do pos = promoter_start to promoter_stop;
  output;
  end;
  else if strand="-" then do pos = promoter_stop to promoter_start;
  output;
  end;
run;


data downstream_2;
  set downstream;
  by chr gene_id;
  if strand="+" then do pos = downstream_start to downstream_stop;
  output;
  end;
  else if strand="-" then do pos = downstream_stop to downstream_start;
  output;
  end;
run;



proc sort data=promoter_2;
  by chr pos;
proc sort data=downstream_2;
  by chr pos;
run;

data meth_sbys10_promoter;
  merge meth_sbys10_2 (in=in1) promoter_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_downstream;
  merge meth_sbys10_2 (in=in1) downstream_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_promoter2;
  set meth_sbys10_promoter;
  if strand="+" then pos_bin=(pos-promoter_start);
  else if strand="-" then pos_bin=(promoter_start-pos) ;
run;

data meth_sbys10_downstream2;
  set meth_sbys10_downstream;
  if strand="+" then pos_bin=(pos-downstream_start) + 2000;
  else if strand="-" then pos_bin=(downstream_start-pos) + 2000;
run;



data meth_sbys10_genestack;
   set meth_sbys10_gene2 (in=in1) meth_sbys10_promoter2 (in=in2) meth_sbys10_downstream2 (in=in3);
run;

proc sort data=meth_sbys10_genestack;
  by site_type gene_id chr pos_bin;
run;

data meth_data2;
  set meth_sbys10_genestack;
  grouped_pos=int(pos_bin/10)*10;
run;


%macro subsetCalc(flag);


proc sort data=meth_data2;
  by site_type pos_bin grouped_pos;
run;


proc means data=meth_data2 noprint;
  by  site_type pos_bin grouped_pos ;
  where &flag.=1;
  var _01Gy _1Gy _0Gy   ;
  output out=meth_data3
  mean=;
run;


proc sort data=meth_data3 ;
  by site_type grouped_pos;
run;


proc means data=meth_data3 noprint;
  by site_type grouped_pos  ;
  var  _01Gy _1Gy _0Gy    
       ;
  output out=meth_data4
  mean=;
run;


data meth_data_export_meth;
  set meth_data4;
  drop _TYPE_ _FREQ_ ;
  rename grouped_pos=pos;
run;

data meth_data_export_acc2;
  set meth_data_export_acc;
  diff_10cGy_acc = _01Gy - _0Gy;
  diff_100cGy_acc = _1Gy - _0Gy;
run;

data meth_data_export_meth2;
  set meth_data_export_meth;
  diff_10cGy_meth = _01Gy - _0Gy;
  diff_100cGy_meth = _1Gy - _0Gy;
run;

proc sort data=meth_data_export_meth2;
  by pos;
proc sort data=meth_data_export_acc2;
  by pos;
run;

data acc2meth;
  merge meth_data_export_acc2 meth_data_export_meth2 ;
  by pos;
run;

proc sort data=acc2meth;
  by site_type pos;
proc transpose data=acc2meth out=acc2meth_sbys(rename=(_NAME_=measure COL1=methylation_difference));
  by site_type pos;
  var diff_: ;
run;

proc sort data=acc2meth_sbys;
  by site_type measure pos ;
run;


data acc2meth_CG_01_&flag. acc2meth_CHG_01_&flag. acc2meth_CHH_01_&flag. acc2meth_CG_1_&flag. acc2meth_CHG_1_&flag. acc2meth_CHH_1_&flag.;
  set acc2meth_sbys;
  if site_type="CG" and (measure="diff_10cGy_acc" or measure="diff_10cGy_meth") then output acc2meth_CG_01_&flag.;
  if site_type="CG" and (measure="diff_100cGy_acc" or measure="diff_100cGy_meth") then output acc2meth_CG_1_&flag.;

  if site_type="CHG" and (measure="diff_10cGy_acc" or measure="diff_10cGy_meth") then output acc2meth_CHG_01_&flag.;
  if site_type="CHG" and (measure="diff_100cGy_acc" or measure="diff_100cGy_meth") then output acc2meth_CHG_1_&flag.;

  if site_type="CHH" and (measure="diff_10cGy_acc" or measure="diff_10cGy_meth") then output acc2meth_CHH_01_&flag.;
  if site_type="CHH" and (measure="diff_100cGy_acc" or measure="diff_100cGy_meth") then output acc2meth_CHH_1_&flag.;
  drop site_type;
run;

%mend;

%subsetCalc(flag_gc_01_up);
%subsetCalc(flag_gc_01_dn);
%subsetCalc(flag_gc_1_up);
%subsetCalc(flag_gc_1_dn);



proc export data=acc2meth_CG_01_flag_gc_01_up
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hyperDA_genes_10cGy_CG_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CG_1_flag_gc_1_up
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hyperDA_genes_100cGy_CG_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CG_01_flag_gc_01_dn
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hypoDA_genes_10cGy_CG_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CG_1_flag_gc_1_dn
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hypoDA_genes_100cGy_CG_meth.csv"
dbms=csv replace;
run;



proc export data=acc2meth_CHG_01_flag_gc_01_up
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hyperDA_genes_10cGy_CHG_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CHG_1_flag_gc_1_up
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hyperDA_genes_100cGy_CHG_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CHG_01_flag_gc_01_dn
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hypoDA_genes_10cGy_CHG_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CHG_1_flag_gc_1_dn
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hypoDA_genes_100cGy_CHG_meth.csv"
dbms=csv replace;
run;



proc export data=acc2meth_CHH_01_flag_gc_01_up
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hyperDA_genes_10cGy_CHH_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CHH_1_flag_gc_1_up
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hyperDA_genes_100cGy_CHH_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CHH_01_flag_gc_01_dn
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hypoDA_genes_10cGy_CHH_meth.csv"
dbms=csv replace;
run;

proc export data=acc2meth_CHH_1_flag_gc_1_dn
outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/acc2meth_hypoDA_genes_100cGy_CHH_meth.csv"
dbms=csv replace;
run;


















%macro makeAccData(list, interGenic, condit1, condit2, siteType);

data gene;
 set &list.;
 where flag_intergenic=&interGenic.;
run;


proc sort data=gene;
  by chr dar_center;
run;

data gene_2;
  set gene;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=gene_2;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data meth_sbys10_gene;
  merge meth_sbys10_2 (in=in1) gene_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_gene2;
  set meth_sbys10_gene;
  distance_to_center=dar_Center-pos;
run;


data meth_sbys10_gene3;
  set meth_sbys10_gene2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;



proc sort data=meth_sbys10_gene3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=meth_sbys10_gene3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy _0Gy _1Gy _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff
  mean(_01Gy)=mean_01Gy
  mean(_1Gy)=mean_1Gy
  mean( _0Gy)=mean_0Gy
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff ;
  by grouped_pos2;
run;


proc means data=mean_diff noprint;
  by grouped_pos2  ;
  var  mean_01Gy mean_1Gy mean_0Gy mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_2
  mean(mean_01Gy)=mean_01Gy
  mean(mean_1Gy)=mean_1Gy
  mean( mean_0Gy)=mean_0Gy
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_1;
  set mean_diff;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_3;
  set mean_diff_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;

data export;
  retain pos mean&condit1. mean&condit2.;
  set mean_diff_3;
  keep pos mean&condit1. mean&condit2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/DAR_&siteType._&list._&condit1._&condit2._flagIntergenic_&intergenic._accessibility.csv"
     dbms=csv replace;
run;

%mend;

%makeAccData(gc2chh_01_up, 0, _0Gy, _01Gy, CHH);
%makeAccData(gc2chh_01_dn, 0, _0Gy, _01Gy, CHH);
%makeAccData(gc2chh_1_up, 0, _0Gy, _1Gy, CHH);
%makeAccData(gc2chh_1_dn, 0, _0Gy, _1Gy, CHH);

%makeAccData(gc2chg_01_up, 0, _0Gy, _01Gy, CHG);
%makeAccData(gc2chg_01_dn, 0, _0Gy, _01Gy, CHG);
%makeAccData(gc2chg_1_up, 0, _0Gy, _1Gy, CHG);
%makeAccData(gc2chg_1_dn, 0, _0Gy, _1Gy, CHG);

%makeAccData(gc2cg_01_up, 0, _0Gy, _01Gy, CG);
%makeAccData(gc2cg_01_dn, 0, _0Gy, _01Gy, CG);
%makeAccData(gc2cg_1_up, 0, _0Gy, _1Gy, CG);
%makeAccData(gc2cg_1_dn, 0, _0Gy, _1Gy, CG);

%makeAccData(gc2chh_01_up, 1, _0Gy, _01Gy, CHH);
%makeAccData(gc2chh_01_dn, 1, _0Gy, _01Gy, CHH);
%makeAccData(gc2chh_1_up, 1, _0Gy, _1Gy, CHH);
%makeAccData(gc2chh_1_dn, 1, _0Gy, _1Gy, CHH);

%makeAccData(gc2chg_01_up, 1, _0Gy, _01Gy, CHG);
%makeAccData(gc2chg_01_dn, 1, _0Gy, _01Gy, CHG);
%makeAccData(gc2chg_1_up, 1, _0Gy, _1Gy, CHG);
%makeAccData(gc2chg_1_dn, 1, _0Gy, _1Gy, CHG);

%makeAccData(gc2cg_01_up, 1, _0Gy, _01Gy, CG);
%makeAccData(gc2cg_01_dn, 1, _0Gy, _01Gy, CG);
%makeAccData(gc2cg_1_up, 1, _0Gy, _1Gy, CG);
%makeAccData(gc2cg_1_dn, 1, _0Gy, _1Gy, CG);




%makeAccData(chh2gc_01_up, 0, _0Gy, _01Gy, CHH);
%makeAccData(chh2gc_01_dn, 0, _0Gy, _01Gy, CHH);
%makeAccData(chh2gc_1_up, 0, _0Gy, _1Gy, CHH);
%makeAccData(chh2gc_1_dn, 0, _0Gy, _1Gy, CHH);
%makeAccData(chg2gc_01_up, 0, _0Gy, _01Gy, CHG);
%makeAccData(chg2gc_01_dn, 0, _0Gy, _01Gy, CHG);
%makeAccData(chg2gc_1_up, 0, _0Gy, _1Gy, CHG);
%makeAccData(chg2gc_1_dn, 0, _0Gy, _1Gy, CHG);
%makeAccData(cg2gc_01_up, 0, _0Gy, _01Gy, CG);
%makeAccData(cg2gc_01_dn, 0, _0Gy, _01Gy, CG);
%makeAccData(cg2gc_1_up, 0, _0Gy, _1Gy, CG);
%makeAccData(cg2gc_1_dn, 0, _0Gy, _1Gy, CG);
%makeAccData(chh2gc_01_up, 1, _0Gy, _01Gy, CHH);
%makeAccData(chh2gc_01_dn, 1, _0Gy, _01Gy, CHH);
%makeAccData(chh2gc_1_up, 1, _0Gy, _1Gy, CHH);
%makeAccData(chh2gc_1_dn, 1, _0Gy, _1Gy, CHH);
%makeAccData(chg2gc_01_up, 1, _0Gy, _01Gy, CHG);
%makeAccData(chg2gc_01_dn, 1, _0Gy, _01Gy, CHG);
%makeAccData(chg2gc_1_up, 1, _0Gy, _1Gy, CHG);
%makeAccData(chg2gc_1_dn, 1, _0Gy, _1Gy, CHG);
%makeAccData(cg2gc_01_up, 1, _0Gy, _01Gy, CG);
%makeAccData(cg2gc_01_dn, 1, _0Gy, _01Gy, CG);
%makeAccData(cg2gc_1_up, 1, _0Gy, _1Gy, CG);
%makeAccData(cg2gc_1_dn, 1, _0Gy, _1Gy, CG);








/* Methylation */

proc freq data=gc2cg_01_up;  tables flag_intergenic; run;
proc freq data=gc2chg_01_up;  tables flag_intergenic; run;
proc freq data=gc2chh_01_up;  tables flag_intergenic; run;
proc freq data=gc2cg_01_dn;  tables flag_intergenic; run;
proc freq data=gc2chg_01_dn;  tables flag_intergenic; run;
proc freq data=gc2chh_01_dn;  tables flag_intergenic; run;
proc freq data=gc2cg_1_up;  tables flag_intergenic; run;
proc freq data=gc2chg_1_up;  tables flag_intergenic; run;
proc freq data=gc2chh_1_up;  tables flag_intergenic; run;
proc freq data=gc2cg_1_dn;  tables flag_intergenic; run;
proc freq data=gc2chg_1_dn;  tables flag_intergenic; run;
proc freq data=gc2chh_1_dn;  tables flag_intergenic; run;


proc freq data=cg2gc_01_up;  tables flag_intergenic; run;
proc freq data=chg2gc_01_up;  tables flag_intergenic; run;
proc freq data=chh2gc_01_up;  tables flag_intergenic; run;
proc freq data=cg2gc_01_dn;  tables flag_intergenic; run;
proc freq data=chg2gc_01_dn;  tables flag_intergenic; run;
proc freq data=chh2gc_01_dn;  tables flag_intergenic; run;
proc freq data=cg2gc_1_up;  tables flag_intergenic; run;
proc freq data=chg2gc_1_up;  tables flag_intergenic; run;
proc freq data=chh2gc_1_up;  tables flag_intergenic; run;
proc freq data=cg2gc_1_dn;  tables flag_intergenic; run;
proc freq data=chg2gc_1_dn;  tables flag_intergenic; run;
proc freq data=chh2gc_1_dn;  tables flag_intergenic; run;




/* prep methylation and accessibility counts */

data meth_data;
  set arabMAP.methylation_data_CG_CHG_CHH;
  keep site_type chr stop_pos treatment units rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by site_type chr stop_pos treatment  rep ;
proc means data=meth_data noprint;
  by site_type chr stop_pos treatment   ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;


proc transpose data=meth_data2 out=meth_sbys10;
  where total_C >= 10;
  by site_type chr stop_pos;
  id treatment  ;
  var perc_methyl;
run;

data meth_sbys10_2;
  set meth_sbys10;
  if _01Gy ne . and _0Gy ne . then _01Gy_common=_01Gy; else _01Gy_common=.;
  if _1Gy ne . and _0Gy ne . then _1Gy_common=_1Gy; else _1Gy_common=.;
  if (_1Gy ne . or _01Gy ne .) and _0Gy ne . then _0Gy_common=_0Gy; else _0Gy_common=.;

  if _01Gy_common ne . and _1Gy_common ne .  then do;
    _01Gy_common_all = _01Gy_common;
    _1Gy_common_all = _1Gy_common;
    _0Gy_common_all = _0Gy_common;
    end;
  else do;
   _01Gy_common_all = .;
   _1Gy_common_all = .;
   _0Gy_common_all = . ;
    end;
  rename stop_pos=pos;
run;


%macro makeMethData(list, interGenic, condit1, condit2, siteType);

data gene;
 set &list.;
 where flag_intergenic=&interGenic.;
run;


proc sort data=gene;
  by chr dar_center;
run;

data gene_2;
  set gene;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

data meth_sbys10_2a;
  set meth_sbys10_2;
  where site_type="&siteType.";
run;


proc sort data=gene_2;
  by chr pos;
proc sort data=meth_sbys10_2a;
  by chr pos;
run;

data meth_sbys10_gene;
  merge meth_sbys10_2a (in=in1) gene_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_gene2;
  set meth_sbys10_gene;
  distance_to_center=dar_Center-pos;
run;


data meth_sbys10_gene3;
  set meth_sbys10_gene2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;



proc sort data=meth_sbys10_gene3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=meth_sbys10_gene3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy _0Gy _1Gy _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff
  mean(_01Gy)=mean_01Gy
  mean(_1Gy)=mean_1Gy
  mean( _0Gy)=mean_0Gy
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff ;
  by grouped_pos2;
run;


proc means data=mean_diff noprint;
  by grouped_pos2  ;
  var  mean_01Gy mean_1Gy mean_0Gy mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_2
  mean(mean_01Gy)=mean_01Gy
  mean(mean_1Gy)=mean_1Gy
  mean( mean_0Gy)=mean_0Gy
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_1;
  set mean_diff;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_3;
  set mean_diff_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;

data export;
  retain pos mean&condit1. mean&condit2.;
  set mean_diff_3;
  keep pos mean&condit1. mean&condit2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/DAR_&siteType._&list._&condit1._&condit2._flagIntergenic_&intergenic._methylation.csv"
     dbms=csv replace;
run;

%mend;

%makeMethData(gc2chh_01_up, 0, _0Gy, _01Gy, CHH);
%makeMethData(gc2chh_01_dn, 0, _0Gy, _01Gy, CHH);
%makeMethData(gc2chh_1_up, 0, _0Gy, _1Gy, CHH);
%makeMethData(gc2chh_1_dn, 0, _0Gy, _1Gy, CHH);

%makeMethData(gc2chg_01_up, 0, _0Gy, _01Gy, CHG);
%makeMethData(gc2chg_01_dn, 0, _0Gy, _01Gy, CHG);
%makeMethData(gc2chg_1_up, 0, _0Gy, _1Gy, CHG);
%makeMethData(gc2chg_1_dn, 0, _0Gy, _1Gy, CHG);

%makeMethData(gc2cg_01_up, 0, _0Gy, _01Gy, CG);
%makeMethData(gc2cg_01_dn, 0, _0Gy, _01Gy, CG);
%makeMethData(gc2cg_1_up, 0, _0Gy, _1Gy, CG);
%makeMethData(gc2cg_1_dn, 0, _0Gy, _1Gy, CG);

%makeMethData(gc2chh_01_up, 1, _0Gy, _01Gy, CHH);
%makeMethData(gc2chh_01_dn, 1, _0Gy, _01Gy, CHH);
%makeMethData(gc2chh_1_up, 1, _0Gy, _1Gy, CHH);
%makeMethData(gc2chh_1_dn, 1, _0Gy, _1Gy, CHH);

%makeMethData(gc2chg_01_up, 1, _0Gy, _01Gy, CHG);
%makeMethData(gc2chg_01_dn, 1, _0Gy, _01Gy, CHG);
%makeMethData(gc2chg_1_up, 1, _0Gy, _1Gy, CHG);
%makeMethData(gc2chg_1_dn, 1, _0Gy, _1Gy, CHG);

%makeMethData(gc2cg_01_up, 1, _0Gy, _01Gy, CG);
%makeMethData(gc2cg_01_dn, 1, _0Gy, _01Gy, CG);
%makeMethData(gc2cg_1_up, 1, _0Gy, _1Gy, CG);
%makeMethData(gc2cg_1_dn, 1, _0Gy, _1Gy, CG);




%makeMethData(chh2gc_01_up, 0, _0Gy, _01Gy, CHH);
%makeMethData(chh2gc_01_dn, 0, _0Gy, _01Gy, CHH);
%makeMethData(chh2gc_1_up, 0, _0Gy, _1Gy, CHH);
%makeMethData(chh2gc_1_dn, 0, _0Gy, _1Gy, CHH);
%makeMethData(chg2gc_01_up, 0, _0Gy, _01Gy, CHG);
%makeMethData(chg2gc_01_dn, 0, _0Gy, _01Gy, CHG);
%makeMethData(chg2gc_1_up, 0, _0Gy, _1Gy, CHG);
%makeMethData(chg2gc_1_dn, 0, _0Gy, _1Gy, CHG);
%makeMethData(cg2gc_01_up, 0, _0Gy, _01Gy, CG);
%makeMethData(cg2gc_01_dn, 0, _0Gy, _01Gy, CG);
%makeMethData(cg2gc_1_up, 0, _0Gy, _1Gy, CG);
%makeMethData(cg2gc_1_dn, 0, _0Gy, _1Gy, CG);
%makeMethData(chh2gc_01_up, 1, _0Gy, _01Gy, CHH);
%makeMethData(chh2gc_01_dn, 1, _0Gy, _01Gy, CHH);
%makeMethData(chh2gc_1_up, 1, _0Gy, _1Gy, CHH);
%makeMethData(chh2gc_1_dn, 1, _0Gy, _1Gy, CHH);
%makeMethData(chg2gc_01_up, 1, _0Gy, _01Gy, CHG);
%makeMethData(chg2gc_01_dn, 1, _0Gy, _01Gy, CHG);
%makeMethData(chg2gc_1_up, 1, _0Gy, _1Gy, CHG);
%makeMethData(chg2gc_1_dn, 1, _0Gy, _1Gy, CHG);
%makeMethData(cg2gc_01_up, 1, _0Gy, _01Gy, CG);
%makeMethData(cg2gc_01_dn, 1, _0Gy, _01Gy, CG);
%makeMethData(cg2gc_1_up, 1, _0Gy, _1Gy, CG);
%makeMethData(cg2gc_1_dn, 1, _0Gy, _1Gy, CG);


