/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

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
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
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
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
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

proc freq data=all_01_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=all_01_dn;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=all_1_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=all_1_dn;  tables flag_DMC10_ge2*flag_fdr05; run;

proc freq data=cg_01_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=cg_01_dn;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=cg_1_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=cg_1_dn;  tables flag_DMC10_ge2*flag_fdr05; run;

proc freq data=chg_01_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=chg_01_dn;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=chg_1_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=chg_1_dn;  tables flag_DMC10_ge2*flag_fdr05; run;

proc freq data=chh_01_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=chh_01_dn;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=chh_1_up;  tables flag_DMC10_ge2*flag_fdr05; run;
proc freq data=chh_1_dn;  tables flag_DMC10_ge2*flag_fdr05; run;


/* COUNTS:

                ALL     CG      CHG     CHH
10cGy Hyper     10086   1501    1135    7450
10cGy Hypo      61638   3148    8065    50425
100cGy Hyper    28646   88      573     27985
100cGy Hypo     16746   309     423     16014

Dup removal:

                ALL     CG      CHG     CHH
10cGy Hyper     10078   1497    1133    7448
10cGy Hypo      61609   3134    8061    50414
100cGy Hyper    28583   88      573     27922
100cGy Hypo     16746   309     423     16014

Split by test
#diffsites | both | binomial

                ALL         CG          CHG         CHH
10cGy Hyper     0|8|10070   0|4|1493    0|2|1131    0|2|7446
10cGy Hypo      0|29|61578  0|14|3119   0|4|8057    0|11|50402
100cGy Hyper    0|63|28519  0|0|88      0|0|573     0|63|27858
100cGy Hypo     0|0|16746   0|0|309     0|0|423     0|0|16014

*/




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

proc freq data=gc_01_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_01_dn;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_dn;  tables flag_DAC10_ge2*flag_fdr05; run;


/* COUNTS:

                
10cGy Hyper     3688
10cGy Hypo      21
100cGy Hyper    772
100cGy Hypo     208

Dup removal:

                
10cGy Hyper     3571
10cGy Hypo      18
100cGy Hyper    772
100cGy Hypo     208

Split by test
#diffsites | both | binomial


10cGy Hyper     88|117|3361
10cGy Hypo      0|3|21
100cGy Hyper    0|0|772
100cGy Hypo     0|0|208

*/


/* for single site * condition outputs, we can leave these alone 
   for combined DMRs, and for hyper/hypo-DMRs/DARs shared between 0.1Gy and 1Gy, we need to merge these into superregions 

   easiest way is to STACK sites, merge yte region with previous is overlapping ( lag code!! )  and track siteType and comparison
   
   for SHARED between 0.1 and 1 Gy y, count the number of occurances of super-region by unique region. If 2 then shared superregion
   */

%macro mergeDMR(dataIN);

proc sort data=&dataIN.;
   by chr region_Start region_stop;
run;

data &dataIN._make_super;
  retain superregion_num;
  set &dataIN.;
  by chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop);
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

proc sort data=&dataIN._make_super;
   by chr superregion_num region_start region_stop;
proc means data=&dataIN._make_super noprint;
   by chr superregion_num;
   var region_start region_stop;
   output out=&dataIN._merged (drop=_TYPE_ _FREQ_) min(region_start)=region_start max(region_stop)=region_stop;
run;

data &dataIN._merged2;
  set &dataIN._merged;
  length site_type $6.;
  site_type="merged";
run;

%mend;

%mergeDMR(all_01_up);
%mergeDMR(all_01_dn);
%mergeDMR(all_1_up);
%mergeDMR(all_1_dn);


/* Merged regions:
10cGy Hyper     9490
10cGy Hypo      53983
100cGy Hyper    28070
100cGy Hypo     16446
*/


%macro commonSuper(dataA, dataB, outName);

data stack_&outName.;
  set &dataA. &dataB.;
  keep site_type chr  region_start region_stop ;
run;

proc sort data=stack_&outName. nodup;
  by site_type chr  region_start region_stop;
run;


data stack_&outName._make_super;
  retain superregion_num;
  set stack_&outName.;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop);
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;


proc sort data=stack_&outName._make_super;
   by site_type chr superregion_num region_start region_stop;
proc means data=stack_&outName._make_super noprint;
   by site_type chr superregion_num;
   var region_start region_stop;
   output out=&outName._shared min(region_start)=region_start max(region_stop)=region_stop;
run;



data &outName._shared2;
  set &outName._shared;
  where _FREQ_ > 1;
  drop _TYPE_ _FREQ_;
run;

%mend;


%commonSuper(all_01_up_merged2, all_1_up_merged2, up_dmr_72);
%commonSuper(all_01_dn_merged2, all_1_dn_merged2, dn_dmr_72);
%commonSuper(cg_01_up, cg_1_up, up_dmr_72_CG);
%commonSuper(cg_01_dn, cg_1_dn, dn_dmr_72_CG);
%commonSuper(chg_01_up, chg_1_up, up_dmr_72_CHG);
%commonSuper(chg_01_dn, chg_1_dn, dn_dmr_72_CHG);
%commonSuper(chh_01_up, chh_1_up, up_dmr_72_CHH);
%commonSuper(chh_01_dn, chh_1_dn, dn_dmr_72_CHH);
%commonSuper(gc_01_up, gc_1_up, up_dar_72);
%commonSuper(gc_01_dn, gc_1_dn, dn_dar_72);

/* Merged super regions
All hyper   3109
All hypo    11908
CG hyper    28
CG hypo     105
ChG hyper   86
CHG hypo    253
CHH hyper   2618
CHH hypo    10975
GC hyper    317
GC hypo     2


*/


/* Export BED files for DREME analysis */

/* Merge hyper/hypo DMR by site type and dose */


data dmr_01_72_cg; set cg_01_up cg_01_dn; run;
data dmr_01_72_chg; set chg_01_up chg_01_dn; run;
data dmr_01_72_chh; set chg_01_up chh_01_dn; run;

data dmr_1_72_cg; set cg_1_up cg_1_dn; run;
data dmr_1_72_chg; set chg_1_up chg_1_dn; run;
data dmr_1_72_chh; set chg_1_up chh_1_dn; run;

proc sort data=dmr_01_72_cg nodup; by _all_; run;
proc sort data=dmr_01_72_chg nodup; by _all_; run;
proc sort data=dmr_01_72_chh nodup; by _all_; run;
proc sort data=dmr_1_72_cg nodup; by _all_; run;
proc sort data=dmr_1_72_chg nodup; by _all_; run;
proc sort data=dmr_1_72_chh nodup; by _all_; run;


%macro exportBED(inputData);

data format_for_bedfile;
   retain chr region_start2 region_stop2 region_id;
   set &inputData.;
   length region_id $100.;
   region_start2 = region_start - 6;
   region_stop2 = region_stop + 6;
   region_id=compress(catx("_",chr,region_start,region_stop));
   keep chr region_start2 region_stop2 region_id;
run;

proc sort data=format_for_bedfile nodup;
  by chr region_start2 region_stop2 region_id;
run;


proc export data=format_for_bedfile
     outfile="!HOME/concannon/DTRA/at_rad_motif_analysis/input_data/min5_&inputData..bed"
     dbms=tab replace;
     putnames=no;
run;

%mend;


%exportBED(cg_01_up);
%exportBED(cg_1_up);
%exportBED(cg_01_dn);
%exportBED(cg_1_dn);


%exportBED(chg_01_up);
%exportBED(chg_1_up);
%exportBED(chg_01_dn);
%exportBED(chg_1_dn);


%exportBED(chh_01_up);
%exportBED(chh_1_up);
%exportBED(chh_01_dn);
%exportBED(chh_1_dn);

%exportBED(gc_01_up);
%exportBED(gc_1_up);
%exportBED(gc_01_dn);
%exportBED(gc_1_dn);


%exportBED(all_01_up_merged2);
%exportBED(all_1_up_merged2);
%exportBED(all_01_dn_merged2);
%exportBED(all_1_dn_merged2);

%exportBED(up_dmr_72_shared2);
%exportBED(dn_dmr_72_shared2);
%exportBED(up_dmr_72_CG_shared2);
%exportBED(dn_dmr_72_CG_shared2);

%exportBED(up_dmr_72_CHG_shared2);
%exportBED(dn_dmr_72_CHG_shared2);
%exportBED(up_dmr_72_CHH_shared2);
%exportBED(dn_dmr_72_CHH_shared2);

%exportBED(up_dar_72_shared2);
%exportBED(dn_dar_72_shared2);



%exportBED(dmr_01_72_cg);
%exportBED(dmr_1_72_cg);
%exportBED(dmr_01_72_chg);
%exportBED(dmr_1_72_chg);
%exportBED(dmr_01_72_chh);
%exportBED(dmr_1_72_chh);

/* DAR counts:
   Hyper all, exon, intron, promoter, downstream; 0.1G, 1G
   Hypo all, exon, intron, promoter, downstream; 0.1G, 1G

 */


/* counts by genic feature */




proc freq data=cg_01_up;  tables feature; run;
proc freq data=cg_01_dn;  tables feature; run;
proc freq data=cg_1_up;  tables feature; run;
proc freq data=cg_1_dn;  tables feature; run;

proc freq data=chg_01_up;  tables feature; run;
proc freq data=chg_01_dn;  tables feature; run;
proc freq data=chg_1_up;  tables feature; run;
proc freq data=chg_1_dn;  tables feature; run;

proc freq data=chh_01_up;  tables feature; run;
proc freq data=chh_01_dn;  tables feature; run;
proc freq data=chh_1_up;  tables feature; run;
proc freq data=chh_1_dn;  tables feature; run;

proc freq data=gc_01_up;  tables feature; run;
proc freq data=gc_01_dn;  tables feature; run;
proc freq data=gc_1_up;  tables feature; run;
proc freq data=gc_1_dn;  tables feature; run;


/* 

Dose    Site    Direction   Intergenic  TTS     exon    intron  promoter   
10cGy   CG      Hyper       472         244     443     69      269
10cGy   CG      Hypo        1509        453     493     122     557
100cGy  CG      Hyper       36          8       16      2       26
100cGy  CG      Hypo        141         50      48      4       66

10cGy   CHG     Hyper       852         85      98      20      78
10cGy   CHG     Hypo        6298        531     264     152     816
100cGy  CHG     Hyper       427         42      27      13      64
100cGy  CHG     Hypo        337         27      22      6       31

10cGy   CHH     Hyper       4748        665     896     202     937
10cGy   CHH     Hypo        33687       4828    1789    1633    8477
100cGy  CHH     Hyper       19138       2602    962     647     4573
100cGy  CHH     Hypo        10773       1527    597     496     2631

10cGy   GC      Hyper       2382        251     525     161     252
10cGy   GC      Hypo        8           3       6       0       1
100cGy  GC      Hyper       568         39      84      31      50
100cGy  GC      Hypo        144         16      24      10      14


*/


/* count by gene */



%macro geneCount(inData);

data &inData._gene;
  set &inData.;
  length geneID2 $32.;
  if feature="Intergenic" then delete;
  geneID2=scan(scan(scan(annotation,2,"("),1,")"),1,".");
  keep geneID2;
run;

proc sort data=&inData._gene nodup;
  by geneID2;
run;

%mend;

%geneCount(cg_01_up); %geneCount(cg_01_dn); 
%geneCount(cg_1_up); %geneCount(cg_1_dn);
%geneCount(chg_01_up); %geneCount(chg_01_dn);
 %geneCount(chg_1_up); %geneCount(chg_1_dn);
%geneCount(chh_01_up); %geneCount(chh_01_dn); 
%geneCount(chh_1_up); %geneCount(chh_1_dn);
%geneCount(gc_01_up); %geneCount(gc_01_dn);
 %geneCount(gc_1_up); %geneCount(gc_1_dn);


/*
            Hyper   Hypo
10cGy   CG  969     1478
100cGy  CG  52      155
10cGy   CHG 257     1424
100cGy  CHG 138     85
10cGy   CHH 2348    8578
100cGy  CHH 5250    3887
10cGy   GC  1120    10
100cGy  GC  202     60
*/

data gene_list;
  merge cg_01_up_gene (in=in1) cg_01_dn_gene (in=in2) cg_1_up_gene (in=in3) cg_1_dn_gene (in=in4)
        chg_01_up_gene (in=in5) chg_01_dn_gene (in=in6) chg_1_up_gene (in=in7) chg_1_dn_gene (in=in8)
        chh_01_up_gene (in=in9) chh_01_dn_gene (in=in10) chh_1_up_gene (in=in11) chh_1_dn_gene (in=in12)
        gc_01_up_gene (in=in13) gc_01_dn_gene (in=in14) gc_1_up_gene (in=in15) gc_1_dn_gene (in=in16);
  by geneID2;
  if in1 then CG_01_UP=1; else CG_01_UP=0;
  if in2 then CG_01_DN=1; else CG_01_DN=0;
  if in3 then CG_1_UP=1; else CG_1_UP=0;
  if in4 then CG_1_DN=1; else CG_1_DN=0;
  if in5 then CHG_01_UP=1; else CHG_01_UP=0;
  if in6 then CHG_01_DN=1; else CHG_01_DN=0;
  if in7 then CHG_1_UP=1; else CHG_1_UP=0;
  if in8 then CHG_1_DN=1; else CHG_1_DN=0;
  if in9 then CHH_01_UP=1; else CHH_01_UP=0;
  if in10 then CHH_01_DN=1; else CHH_01_DN=0;
  if in11 then CHH_1_UP=1; else CHH_1_UP=0;
  if in12 then CHH_1_DN=1; else CHH_1_DN=0;
  if in13 then GC_01_UP=1; else GC_01_UP=0;
  if in14 then GC_01_DN=1; else GC_01_DN=0;
  if in15 then GC_1_UP=1; else GC_1_UP=0;
  if in16 then GC_1_DN=1; else GC_1_DN=0;
run;


proc freq data=gene_list noprint;
  tables CG_01_UP*CG_01_DN*CG_1_UP*CG_1_DN / out=gene_count_CG ;
  tables CHG_01_UP*CHG_01_DN*CHG_1_UP*CHG_1_DN / out=gene_count_CHG;
  tables CHH_01_UP*CHH_01_DN*CHH_1_UP*CHH_1_DN / out=gene_count_CHH;
  tables GC_01_UP*GC_01_DN*GC_1_UP*GC_1_DN  / out=gene_count_GC;
run;


proc print data=gene_count_CG;
proc print data=gene_count_CHG;
proc print data=gene_count_CHH;
proc print data=gene_count_GC;
run;


/*
CG gene crosstabs:

 CG_01_UP    CG_01_DN    CG_1_UP    CG_1_DN    COUNT

     0           0          0          1          56
     0           0          1          0          18
     0           1          0          0        1288
     0           1          0          1          85
     0           1          1          0           6
     0           1          1          1           1
     1           0          0          0         843
     1           0          0          1           7
     1           0          1          0          21
     1           1          0          0          86
     1           1          0          1           6
     1           1          1          0           6


CHG gene crosstabs:


   CHG_     CHG_
  01_UP    01_DN    CHG_1_UP    CHG_1_DN    COUNT
    0        0          0           1          15
    0        0          1           0          61
    0        0          1           1           1
    0        1          0           0        1261
    0        1          0           1          54
    0        1          1           0          41
    0        1          1           1           8
    1        0          0           0         176
    1        0          0           1           3
    1        0          1           0          18
    1        1          0           0          48
    1        1          0           1           3
    1        1          1           0           8
    1        1          1           1           1

CHH gene crosstabs:


   CHH_     CHH_
  01_UP    01_DN    CHH_1_UP    CHH_1_DN    COUNT
    0        0          0           1         435
    0        0          1           0        1146
    0        0          1           1          68
    0        1          0           0        3104
    0        1          0           1        1274
    0        1          1           0        1548
    0        1          1           1        1372
    1        0          0           0         811
    1        0          0           1          48
    1        0          1           0         190
    1        0          1           1          19
    1        1          0           0         218
    1        1          0           1         155
    1        1          1           0         391
    1        1          1           1         516


GC gene crosstabs:


 GC_01_UP    GC_01_DN    GC_1_UP    GC_1_DN    COUNT
     0           0          0          1          48
     0           0          1          0          91
     0           1          0          0           6
     0           1          0          1           1
     0           1          1          0           1
     1           0          0          0         999
     1           0          0          1          10
     1           0          1          0         109
     1           0          1          1           1
     1           1          0          0           1



*/



/* 4 way venn of regions */


%macro fourWayCtab(inData1,inData2,inData3,inData4);

data stack_all;
  set &inData1. (in=in1) &inData2. (in=in2) &inData3. (in=in3) &inData4. (in=in4);
  length comp $20.;
  if in1 then comp="&inData1.";
  if in2 then comp="&inData2.";
  if in3 then comp="&inData3.";
  if in4 then comp="&inData4.";
  keep comp site_type chr region_start region_stop ;
run;

proc sort data=stack_all nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_all_make_super;
  retain superregion_num;
  set stack_all;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop > region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_all_super1;
   set stack_all_make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_all_super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_all_super1 out=stack_all_super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;


data stack_all_super_sbys2;
  set stack_all_super_sbys;
  if &inData1.=. then &inData1.=0;
  if &inData2.=. then &inData2.=0;
  if &inData3.=. then &inData3.=0;
  if &inData4.=. then &inData4.=0;
run;

proc freq data=stack_all_super_sbys2 noprint;
  tables &inData1.*&inData2.*&inData3.*&inData4. / out=region_compare;
run;
proc print data=region_compare;
run;

%mend;

%fourWayCtab(cg_01_up,cg_01_dn,cg_1_up,cg_1_dn);
%fourWayCtab(chg_01_up,chg_01_dn,chg_1_up,chg_1_dn);
%fourWayCtab(chh_01_up,chh_01_dn,chh_1_up,chh_1_dn);
%fourWayCtab(gc_01_up,gc_01_dn,gc_1_up,gc_1_dn);



/* Used the program "intervene" to do this instead... */





/* Export data for making CDF plots of length in python */


data cg_01_up_length; set cg_01_up; keep chr region_num region_length; run;
data cg_01_dn_length; set cg_01_dn; keep chr region_num region_length; run;
data cg_01_length; set cg_01_up cg_01_dn; keep chr region_num region_length; run;
data cg_1_up_length; set cg_1_up; keep chr region_num region_length; run;
data cg_1_dn_length; set cg_1_dn; keep chr region_num region_length; run;
data cg_1_length; set cg_1_up cg_1_dn; keep chr region_num region_length; run;

data chg_01_up_length; set chg_01_up; keep chr region_num region_length; run;
data chg_01_dn_length; set chg_01_dn; keep chr region_num region_length; run;
data chg_01_length; set chg_01_up chg_01_dn; keep chr region_num region_length; run;
data chg_1_up_length; set chg_1_up; keep chr region_num region_length; run;
data chg_1_dn_length; set chg_1_dn; keep chr region_num region_length; run;
data chg_1_length; set chg_1_up chg_1_dn; keep chr region_num region_length; run;

data chh_01_up_length; set chh_01_up; keep chr region_num region_length; run;
data chh_01_dn_length; set chh_01_dn; keep chr region_num region_length; run;
data chh_01_length; set chh_01_up chh_01_dn; keep chr region_num region_length; run;
data chh_1_up_length; set chh_1_up; keep chr region_num region_length; run;
data chh_1_dn_length; set chh_1_dn; keep chr region_num region_length; run;
data chh_1_length; set chh_1_up chh_1_dn; keep chr region_num region_length; run;

data gc_01_up_length; set gc_01_up; keep chr region_num region_length; run;
data gc_01_dn_length; set gc_01_dn; keep chr region_num region_length; run;
data gc_01_length; set gc_01_up gc_01_dn; keep chr region_num region_length; run;
data gc_1_up_length; set gc_1_up; keep chr region_num region_length; run;
data gc_1_dn_length; set gc_1_dn; keep chr region_num region_length; run;
data gc_1_length; set gc_1_up gc_1_dn; keep chr region_num region_length; run;


data cg_all;
  set cg_01_up (in=in1) cg_01_dn (in=in2) cg_1_up (in=in3) cg_1_dn (in=in4);
  length dmr_set $20.;
  if in1 then dmr_set="CG_01";
  if in2 then dmr_set="CG_01";
  if in3 then dmr_set="CG_1";
  if in4 then dmr_set="CG_1";
run;
data chg_all;
  set chg_01_up (in=in1) chg_01_dn (in=in2) chg_1_up (in=in3) chg_1_dn (in=in4);
  length dmr_set $20.;
  if in1 then dmr_set="CHG_01";
  if in2 then dmr_set="CHG_01";
  if in3 then dmr_set="CHG_1";
  if in4 then dmr_set="CHG_1";
run;
data chh_all;
  set chh_01_up (in=in1) chh_01_dn (in=in2) chh_1_up (in=in3) chh_1_dn (in=in4);
  length dmr_set $20.;
  if in1 then dmr_set="CHH_01";
  if in2 then dmr_set="CHH_01";
  if in3 then dmr_set="CHH_1";
  if in4 then dmr_set="CHH_1";
run;
data gc_all;
  set gc_01_up (in=in1) gc_01_dn (in=in2) gc_1_up (in=in3) gc_1_dn (in=in4);
  length dmr_set $20.;
  if in1 then dmr_set="GC_01";
  if in2 then dmr_set="GC_01";
  if in3 then dmr_set="GC_1";
  if in4 then dmr_set="GC_1";
run;

proc sort data=cg_all;
 by dmr_set;
proc sort data=chg_all;
 by dmr_set;
proc sort data=chh_all;
 by dmr_set;
proc sort data=gc_all;
 by dmr_set;
proc means data=cg_all noprint;
  by dmr_set ;
  var region_length;
  output out=cg_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;
proc means data=chg_all noprint;
  by dmr_set ;
  var region_length;
  output out=chg_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;
proc means data=chh_all noprint;
  by dmr_set ;
  var region_length;
  output out=chh_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;
proc means data=gc_all noprint;
  by dmr_set ;
  var region_length;
  output out=gc_distrib mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;

proc anova data=cg_all;
 class dmr_set;
 model region_length  = dmr_set;
 ods output modelanova=cg_len_test;
run;

proc anova data=chg_all;
 class dmr_set;
 model region_length  = dmr_set;
 ods output modelanova=chg_len_test;
run;

proc anova data=chh_all;
 class dmr_set;
 model region_length  = dmr_set;
 ods output modelanova=chh_len_test;
run;

proc anova data=gc_all;
 class dmr_set;
 model region_length  = dmr_set;
 ods output modelanova=gc_len_test;
run;



proc export data=cg_01_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_10cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=cg_01_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_10cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=cg_01_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_10cGy_all_length.txt" dbms=tab replace; run;
proc export data=cg_1_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_100cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=cg_1_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_100cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=cg_1_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_100cGy_all_length.txt" dbms=tab replace; run;

proc export data=chg_01_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_10cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=chg_01_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_10cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=chg_01_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_10cGy_all_length.txt" dbms=tab replace; run;
proc export data=chg_1_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_100cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=chg_1_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_100cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=chg_1_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_100cGy_all_length.txt" dbms=tab replace; run;

proc export data=chh_01_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_10cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=chh_01_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_10cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=chh_01_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_10cGy_all_length.txt" dbms=tab replace; run;
proc export data=chh_1_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_100cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=chh_1_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_100cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=chh_1_length outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_100cGy_all_length.txt" dbms=tab replace; run;

proc export data=gc_01_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_10cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=gc_01_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_10cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=gc_01_length outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_10cGy_all_length.txt" dbms=tab replace; run;
proc export data=gc_1_up_length outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_100cGy_hyper_length.txt" dbms=tab replace; run;
proc export data=gc_1_dn_length outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_100cGy_hypo_length.txt" dbms=tab replace; run;
proc export data=gc_1_length outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_100cGy_all_length.txt" dbms=tab replace; run;




/* CDF plot data for all regions */

data cg_length_01 chg_length_01 chh_length_01
     cg_length_1  chg_length_1  chh_length_1;
  set dmr_w_annot2;
  if site_type = "CG" and comparison = "0Gy_vs_01G" then output cg_length_01;
  if site_type = "CHG" and comparison = "0Gy_vs_01G" then output chg_length_01;
  if site_type = "CHH" and comparison = "0Gy_vs_01G" then output chh_length_01;
  if site_type = "CG" and comparison = "0Gy_vs_1Gy" then output cg_length_1;
  if site_type = "CHG" and comparison = "0Gy_vs_1Gy" then output chg_length_1;
  if site_type = "CHH" and comparison = "0Gy_vs_1Gy" then output chh_length_1;
  keep chr region_num region_length;
run;

data gc_length_01 gc_length_1 ;
  set dar_w_annot2;
  if comparison = "01Gy_0Gy" then output gc_length_01;
  if comparison = "1Gy_0Gy" then output gc_length_1;
  keep chr region_num region_length;
run;

proc export data=cg_length_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_10cGy_all_regions_length.txt"
     dbms=tab replace;
     run;

proc export data=chg_length_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_10cGy_all_regions_length.txt"
     dbms=tab replace;
     run;

proc export data=chh_length_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_10cGy_all_regions_length.txt"
     dbms=tab replace;
     run;

proc export data=cg_length_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_100cGy_all_regions_length.txt"
     dbms=tab replace;
     run;

proc export data=chg_length_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_100cGy_all_regions_length.txt"
     dbms=tab replace;
     run;

proc export data=chh_length_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_100cGy_all_regions_length.txt"
     dbms=tab replace;
     run;


proc export data=gc_length_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_10cGy_all_regions_length.txt"
     dbms=tab replace;
     run;

proc export data=gc_length_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/GC_100cGy_all_regions_length.txt"
     dbms=tab replace;
     run;




/* DAGs vs DEGs */


data de_results;
  set arabRNA.arab_results_by_gene;
  keep gene_id flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
  flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
  flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
  mean_cpm_: 
  fdr_: 
  flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_3h_fdr05
  flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_72h_fdr05
  flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_3h_fdr05
  flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_72h_fdr05 ;
run;


data de_results2;
  set de_results;
  log2fc_01gy_Mock_1h = log2(mean_cpm_01gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_01gy_Mock_3h = log2(mean_cpm_01gy_3h ) - log2(mean_cpm_Mock_3h );
  log2fc_01gy_Mock_24h = log2(mean_cpm_01gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_01gy_Mock_72h = log2(mean_cpm_01gy_72h ) - log2(mean_cpm_Mock_72h );
  log2fc_1gy_Mock_1h = log2(mean_cpm_1gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_1gy_Mock_3h = log2(mean_cpm_1gy_3h) - log2(mean_cpm_Mock_3h );
  log2fc_1gy_Mock_24h = log2(mean_cpm_1gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_1gy_Mock_72h = log2(mean_cpm_1gy_72h ) - log2(mean_cpm_Mock_72h );
run;



data up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
     dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1
     up_1_1_fc1  up_1_3_fc1  up_1_24_fc1  up_1_72_fc1
     dn_1_1_fc1  dn_1_3_fc1  dn_1_24_fc1  dn_1_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;
keep gene_id;
run;

proc sort data=up_01_1_fc1; by gene_id;
 proc sort data=up_01_3_fc1; by gene_id;
  proc sort data=up_01_24_fc1; by gene_id;
    proc sort data=up_01_72_fc1; by gene_id;

 proc sort data=dn_01_1_fc1; by gene_id;
  proc sort data=dn_01_3_fc1; by gene_id;
   proc sort data=dn_01_24_fc1; by gene_id;
    proc sort data=dn_01_72_fc1; by gene_id;

  proc sort data=up_1_1_fc1; by gene_id;
   proc sort data=up_1_3_fc1; by gene_id;
   proc sort data=up_1_24_fc1; by gene_id;
  proc sort data=up_1_72_fc1; by gene_id;

  proc sort data=dn_1_1_fc1; by gene_id;
   proc sort data=dn_1_3_fc1; by gene_id;
   proc sort data=dn_1_24_fc1; by gene_id;
    proc sort data=dn_1_72_fc1; by gene_id;
run;


%macro countDAGDEG(inDEG, inDAG);

data stack_DEG;
   set dn_&inDEG._fc1 up_&inDEG._fc1;
run;

data stack_DAG;
   set gc_&inDAG._dn_gene gc_&inDAG._up_gene;
   length gene_id $32.;
   gene_id = upcase(geneID2);
run;

proc sort data=stack_DEG nodup;
by gene_id;
run;

proc sort data=stack_DAG nodup;
by gene_id;
run;


data deg_vs_dag;
  merge stack_DAG (in=in1) stack_DEG (in=in2);
  by gene_id;
  if in1 then flag_DAG=1; else flag_DAG=0;
  if in2 then flag_DEG=1; else flag_DEG=0;
run;

proc freq data=deg_vs_dag;
  tables flag_DAG*flag_DEG;
run;

%mend;


%countDAGDEG(01_1, 01);
%countDAGDEG(01_3, 01);
%countDAGDEG(01_24, 01);
%countDAGDEG(01_72, 01);

data stack_DEG;
  set up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
      dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1;
run;

data stack_DAG;
   set gc_01_dn_gene gc_01_up_gene;
   length gene_id $32.;
   gene_id = upcase(geneID2);
run;

proc sort data=stack_DEG nodup;
  by gene_id;
proc sort data=stack_DAG;
  by gene_id;
run;

data deg_vs_dag;
  merge stack_DAG (in=in1) stack_DEG (in=in2);
  by gene_id;
  if in1 then flag_DAG=1; else flag_DAG=0;
  if in2 then flag_DEG=1; else flag_DEG=0;
run;

proc freq data=deg_vs_dag;
  tables flag_DAG*flag_DEG;
run;



%countDAGDEG(1_1, 1);
%countDAGDEG(1_3, 1);
%countDAGDEG(1_24, 1);
%countDAGDEG(1_72, 1);



data stack_DEG;
  set up_1_1_fc1 up_1_3_fc1 up_1_24_fc1 up_1_72_fc1
      dn_1_1_fc1 dn_1_3_fc1 dn_1_24_fc1 dn_1_72_fc1;
run;

data stack_DAG;
   set gc_1_dn_gene gc_1_up_gene;
   length gene_id $32.;
   gene_id = upcase(geneID2);
run;

proc sort data=stack_DEG nodup;
  by gene_id;
proc sort data=stack_DAG;
  by gene_id;
run;

data deg_vs_dag;
  merge stack_DAG (in=in1) stack_DEG (in=in2);
  by gene_id;
  if in1 then flag_DAG=1; else flag_DAG=0;
  if in2 then flag_DEG=1; else flag_DEG=0;
run;

proc freq data=deg_vs_dag;
  tables flag_DAG*flag_DEG;
run;


/* 

Dose    Time    DAG     Both    DE
10cGy   1h      1119    9       489
10cGy   3h      1127    1       32
10cGy   24h     1126    2       30
10cGy   72h     1127    1       65
10cGy   Any     1117    12      589
100cGy  1h      258     3       340
100cGy  3h      259     2       164
100cGy  24h     261     0       5
100cGy  72h     261     0       44
100cGy  Any     258     4       532

*/

/* DMR vs DAR */



%macro countBySite(siteType);
%macro countDMGDEG(inDEG, inDAG);

data stack_DEG;
   set dn_&inDEG._fc1 up_&inDEG._fc1;
run;

data stack_DAG;
   set &siteType._&inDAG._dn_gene &siteType._&inDAG._up_gene;
   length gene_id $32.;
   gene_id = upcase(geneID2);
run;

proc sort data=stack_DEG nodup;
by gene_id;
run;

proc sort data=stack_DAG nodup;
by gene_id;
run;


data deg_vs_dag;
  merge stack_DAG (in=in1) stack_DEG (in=in2);
  by gene_id;
  if in1 then flag_DAG=1; else flag_DAG=0;
  if in2 then flag_DEG=1; else flag_DEG=0;
run;

proc freq data=deg_vs_dag;
  tables flag_DAG*flag_DEG;
run;

%mend;


%countDMGDEG(01_1, 01);
%countDMGDEG(01_3, 01);
%countDMGDEG(01_24, 01);
%countDMGDEG(01_72, 01);

data stack_DEG;
  set up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
      dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1;
run;

data stack_DAG;
   set &siteType._01_dn_gene &siteType._01_up_gene;
   length gene_id $32.;
   gene_id = upcase(geneID2);
run;

proc sort data=stack_DEG nodup;
  by gene_id;
proc sort data=stack_DAG nodup;
  by gene_id;
run;

data deg_vs_dag;
  merge stack_DAG (in=in1) stack_DEG (in=in2);
  by gene_id;
  if in1 then flag_DAG=1; else flag_DAG=0;
  if in2 then flag_DEG=1; else flag_DEG=0;
run;

proc freq data=deg_vs_dag;
  tables flag_DAG*flag_DEG;
run;



%countDMGDEG(1_1, 1);
%countDMGDEG(1_3, 1);
%countDMGDEG(1_24, 1);
%countDMGDEG(1_72, 1);



data stack_DEG;
  set up_1_1_fc1 up_1_3_fc1 up_1_24_fc1 up_1_72_fc1
      dn_1_1_fc1 dn_1_3_fc1 dn_1_24_fc1 dn_1_72_fc1;
run;

data stack_DAG;
   set &siteType._1_dn_gene &siteType._1_up_gene;
   length gene_id $32.;
   gene_id = upcase(geneID2);
run;

proc sort data=stack_DEG nodup;
  by gene_id;
proc sort data=stack_DAG nodup;
  by gene_id;
run;

data deg_vs_dag;
  merge stack_DAG (in=in1) stack_DEG (in=in2);
  by gene_id;
  if in1 then flag_DAG=1; else flag_DAG=0;
  if in2 then flag_DEG=1; else flag_DEG=0;
run;

proc freq data=deg_vs_dag;
  tables flag_DAG*flag_DEG;
run;

%mend;


%countBySite(cg);
%countBySite(chg);
%countBySite(chh);



/*

Type    Dose    Time    DMG     Both    DE
CG      10cGy   1h      2320    29      469
CG      10cGy   3h      2347    2       31
CG      10cGy   24h     2348    1       31
CG      10cGy   72h     2348    1       65
CG      10cGy   Any     2316    33      568
CG      100cGy  1h      205     1       342
CG      100cGy  3h      204     2       164
CG      100cGy  24h     206     0       5
CG      100cGy  72h     206     0       44
CG      100cGy  Any     203     3       533
CHG     10cGy   1h      1601    20      478
CHG     10cGy   3h      1619    2       31
CHG     10cGy   24h     1619    2       30
CHG     10cGy   72h     1620    1       65
CHG     10cGy   Any     1598    23      578
CHG     100cGy  1h      211     2       341
CHG     100cGy  3h      213     0       166
CHG     100cGy  24h     213     0       5
CHG     100cGy  72h     213     0       44
CHG     100cGy  Any     211     2       534
CHH     10cGy   1h      9520    126     372
CHH     10cGy   3h      9632    14      19
CHH     10cGy   24h     9639    7       25
CHH     10cGy   72h     9629    17      49
CHH     10cGy   Any     9492    154     447
CHH     100cGy  1h      7089    73      270
CHH     100cGy  3h      7128    34      132
CHH     100cGy  24h     7162    9       5
CHH     100cGy  72h     7153    9       35
CHH     100cGy  Any     7051    111     425

*/


/*   TE analysis */


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_DMR_to_TE_min_5_sites.txt"
out=dmr_te_annot dbms=tab replace;
guessingrows=all;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_DAR_to_TE_min_5_sites.txt"
out=dar_te_annot dbms=tab replace;
guessingrows=all;
run;




data dmr_te_annot2;
  set dmr_te_annot;
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

data dar_te_annot2;
  set dar_te_annot;
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
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;


data results_by_dmr;
   set wgbsA.results_by_dmr_5sites;
run;

data results_by_dar;
   set wgbsA.results_by_dar_5sites;
run;


proc sort data=dmr_te_annot2;
  by comparison site_type chr  region_num;
proc sort data=dar_te_annot2;
  by comparison site_type chr  region_num;
proc sort data=results_by_dmr;
  by comparison site_type chr  region_num;
proc sort data=results_by_dar;
  by comparison site_type chr  region_num;
run;

data dmr_w_te_annot;
  merge dmr_te_annot2 (in=in1) results_by_dmr (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

data dar_w_te_annot;
  merge dar_te_annot2 (in=in1) results_by_dar (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


data dmr_w_te_annot2;
  set dmr_w_te_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;


data dar_w_te_annot2;
  set dar_w_te_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;





data all_01_up all_01_dn all_1_up all_1_dn
     cg_01_up cg_01_dn  cg_1_up cg_1_dn
     chg_01_up chg_01_dn  chg_1_up chg_1_dn
     chh_01_up chh_01_dn  chh_1_up chh_1_dn;
     set dmr_w_te_annot2;
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


proc freq data=cg_01_up;  tables feature; run;
proc freq data=cg_01_dn;  tables feature; run;
proc freq data=cg_1_up;  tables feature; run;
proc freq data=cg_1_dn;   tables feature; run;

proc freq data=chg_01_up;  tables feature; run;
proc freq data=chg_01_dn;  tables feature; run;
proc freq data=chg_1_up;  tables feature; run;
proc freq data=chg_1_dn;   tables feature; run;

proc freq data=chh_01_up;  tables feature; run;
proc freq data=chh_01_dn;  tables feature; run;
proc freq data=chh_1_up;  tables feature; run;
proc freq data=chh_1_dn;  tables feature; run;






/* Count DARs */

data gc_01_up gc_01_dn gc_1_up gc_1_dn;
     set dar_w_te_annot2;
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

proc freq data=gc_01_up;  tables feature; run;
proc freq data=gc_01_dn;  tables feature; run;
proc freq data=gc_1_up; tables feature; run;
proc freq data=gc_1_dn;  tables feature; run;



/* COUNTS:

Dose    Sign    Type    Intergenic  TTS     Gene    Promoter
10cGy   Hyper   CG      804         192     208     297
10cGy   Hypo    CG      1080        464     669     935
100cGy  Hyper   CG      35          20      18      15
100cGy  Hypo    CG      111         67      64      67
10cGy   Hyper   CHG     171         168     455     341
10cGy   Hypo    CHG     641         1326    3286    2812
100cGy  Hyper   CHG     51          95      252     175
100cGy  Hypo    CHG     33          65      185     140
10cGy   Hyper   CHH     1654        1132    2292    2372
10cGy   Hypo    CHH     7166        8678    14350   20231
100cGy  Hyper   CHH     3207        5002    8642    11134
100cGy  Hypo    CHH     2310        2691    4581    6432
10cGy   Hyper   GC      891         516     1166    998
10cGy   Hypo    GC      8           2       5       3
100cGy  Hyper   GC      134         129     267     242
100cGy  Hypo    GC      39          29      84      56

*/


data dmr_w_te_annot3;
  set dmr_w_te_annot2;
  keep site_type comparison chr region_num feature;
  rename feature=TE_feature;
run;

data dar_w_te_annot3;
  set dar_w_te_annot2;
  keep site_type comparison chr region_num feature;
  rename feature=TE_feature;
run;

data dmr_w_annot3;
  set dmr_w_annot2;
  *keep site_type comparison chr region_num feature;
  rename feature=genic_feature;
run;

data dar_w_annot3;
  set dar_w_annot2;
  *keep site_type comparison chr region_num feature;
  rename feature=genic_feature;
run;

proc sort data=dmr_w_annot3;
  by comparison site_Type chr region_num;
proc sort data=dar_w_annot3;
  by comparison site_Type chr region_num;
proc sort data=dmr_w_te_annot3;
  by comparison site_Type chr region_num;
proc sort data=dar_w_te_annot3;
  by comparison site_Type chr region_num;
run;

data dmr_w_annot4;
  merge dmr_w_annot3 (in=in1) dmr_w_te_annot3 (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

data dar_w_annot4;
  merge dar_w_annot3 (in=in1) dar_w_te_annot3 (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;





data all_01_up all_01_dn all_1_up all_1_dn
     cg_01_up cg_01_dn  cg_1_up cg_1_dn
     chg_01_up chg_01_dn  chg_1_up chg_1_dn
     chh_01_up chh_01_dn  chh_1_up chh_1_dn;
     set dmr_w_annot4;
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
     set dar_w_annot4;
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





proc freq data=cg_01_up;  tables genic_feature*TE_feature; run;
proc freq data=cg_01_dn;  tables genic_feature*TE_feature; run;
proc freq data=cg_1_up;  tables genic_feature*TE_feature; run;
proc freq data=cg_1_dn;   tables genic_feature*TE_feature; run;

proc freq data=chg_01_up;  tables genic_feature*TE_feature; run;
proc freq data=chg_01_dn;  tables genic_feature*TE_feature; run;
proc freq data=chg_1_up;  tables genic_feature*TE_feature; run;
proc freq data=chg_1_dn;   tables genic_feature*TE_feature; run;

proc freq data=chh_01_up;  tables genic_feature*TE_feature; run;
proc freq data=chh_01_dn;  tables genic_feature*TE_feature; run;
proc freq data=chh_1_up;  tables genic_feature*TE_feature; run;
proc freq data=chh_1_dn;  tables genic_feature*TE_feature; run;


proc freq data=gc_01_up;  tables genic_feature*TE_feature; run;
proc freq data=gc_01_dn;  tables genic_feature*TE_feature; run;
proc freq data=gc_1_up; tables genic_feature*TE_feature; run;
proc freq data=gc_1_dn;  tables genic_feature*TE_feature; run;


/*

10cGy CG Hyper:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     34 |    105 |    180 |    153 |    472
             |   2.27 |   7.01 |  12.02 |  10.22 |  31.53
             |   7.20 |  22.25 |  38.14 |  32.42 |
             |   4.24 |  54.69 |  86.54 |  51.69 |
-------------+--------+--------+--------+--------+
TTS          |    176 |     24 |     11 |     33 |    244
             |  11.76 |   1.60 |   0.73 |   2.20 |  16.30
             |  72.13 |   9.84 |   4.51 |  13.52 |
             |  21.97 |  12.50 |   5.29 |  11.15 |
-------------+--------+--------+--------+--------+
exon         |    376 |     30 |      0 |     37 |    443
             |  25.12 |   2.00 |   0.00 |   2.47 |  29.59
             |  84.88 |   6.77 |   0.00 |   8.35 |
             |  46.94 |  15.63 |   0.00 |  12.50 |
-------------+--------+--------+--------+--------+
intron       |     53 |      4 |      4 |      8 |     69
             |   3.54 |   0.27 |   0.27 |   0.53 |   4.61
             |  76.81 |   5.80 |   5.80 |  11.59 |
             |   6.62 |   2.08 |   1.92 |   2.70 |
-------------+--------+--------+--------+--------+
promoter-TSS |    162 |     29 |     13 |     65 |    269
             |  10.82 |   1.94 |   0.87 |   4.34 |  17.97
             |  60.22 |  10.78 |   4.83 |  24.16 |
             |  20.22 |  15.10 |   6.25 |  21.96 |
-------------+--------+--------+--------+--------+
Total             801      192      208      296     1497
                53.51    12.83    13.89    19.77   100.00

10cGy CG Hypo:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     80 |    262 |    583 |    584 |   1509
             |   2.55 |   8.36 |  18.60 |  18.63 |  48.15
             |   5.30 |  17.36 |  38.63 |  38.70 |
             |   7.50 |  56.59 |  87.14 |  62.46 |
-------------+--------+--------+--------+--------+
TTS          |    243 |     69 |     38 |    103 |    453
             |   7.75 |   2.20 |   1.21 |   3.29 |  14.45
             |  53.64 |  15.23 |   8.39 |  22.74 |
             |  22.77 |  14.90 |   5.68 |  11.02 |
-------------+--------+--------+--------+--------+
exon         |    394 |     41 |      3 |     55 |    493
             |  12.57 |   1.31 |   0.10 |   1.75 |  15.73
             |  79.92 |   8.32 |   0.61 |  11.16 |
             |  36.93 |   8.86 |   0.45 |   5.88 |
-------------+--------+--------+--------+--------+
intron       |     91 |      8 |      6 |     17 |    122
             |   2.90 |   0.26 |   0.19 |   0.54 |   3.89
             |  74.59 |   6.56 |   4.92 |  13.93 |
             |   8.53 |   1.73 |   0.90 |   1.82 |
-------------+--------+--------+--------+--------+
promoter-TSS |    259 |     83 |     39 |    176 |    557
             |   8.26 |   2.65 |   1.24 |   5.62 |  17.77
             |  46.50 |  14.90 |   7.00 |  31.60 |
             |  24.27 |  17.93 |   5.83 |  18.82 |
-------------+--------+--------+--------+--------+
Total            1067      463      669      935     3134
                34.05    14.77    21.35    29.83   100.00




100cGy CG Hyper:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |      1 |     10 |     16 |      9 |     36
             |   1.14 |  11.36 |  18.18 |  10.23 |  40.91
             |   2.78 |  27.78 |  44.44 |  25.00 |
             |   2.86 |  50.00 |  88.89 |  60.00 |
-------------+--------+--------+--------+--------+
TTS          |      3 |      2 |      1 |      2 |      8
             |   3.41 |   2.27 |   1.14 |   2.27 |   9.09
             |  37.50 |  25.00 |  12.50 |  25.00 |
             |   8.57 |  10.00 |   5.56 |  13.33 |
-------------+--------+--------+--------+--------+
exon         |     12 |      3 |      0 |      1 |     16
             |  13.64 |   3.41 |   0.00 |   1.14 |  18.18
             |  75.00 |  18.75 |   0.00 |   6.25 |
             |  34.29 |  15.00 |   0.00 |   6.67 |
-------------+--------+--------+--------+--------+
intron       |      1 |      1 |      0 |      0 |      2
             |   1.14 |   1.14 |   0.00 |   0.00 |   2.27
             |  50.00 |  50.00 |   0.00 |   0.00 |
             |   2.86 |   5.00 |   0.00 |   0.00 |
-------------+--------+--------+--------+--------+
promoter-TSS |     18 |      4 |      1 |      3 |     26
             |  20.45 |   4.55 |   1.14 |   3.41 |  29.55
             |  69.23 |  15.38 |   3.85 |  11.54 |
             |  51.43 |  20.00 |   5.56 |  20.00 |
-------------+--------+--------+--------+--------+
Total              35       20       18       15       88
                39.77    22.73    20.45    17.05   100.00






100cGy CG Hypo:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     11 |     36 |     56 |     38 |    141
             |   3.56 |  11.65 |  18.12 |  12.30 |  45.63
             |   7.80 |  25.53 |  39.72 |  26.95 |
             |   9.91 |  53.73 |  87.50 |  56.72 |
-------------+--------+--------+--------+--------+
TTS          |     27 |     10 |      4 |      9 |     50
             |   8.74 |   3.24 |   1.29 |   2.91 |  16.18
             |  54.00 |  20.00 |   8.00 |  18.00 |
             |  24.32 |  14.93 |   6.25 |  13.43 |
-------------+--------+--------+--------+--------+
exon         |     31 |      9 |      1 |      7 |     48
             |  10.03 |   2.91 |   0.32 |   2.27 |  15.53
             |  64.58 |  18.75 |   2.08 |  14.58 |
             |  27.93 |  13.43 |   1.56 |  10.45 |
-------------+--------+--------+--------+--------+
intron       |      3 |      0 |      1 |      0 |      4
             |   0.97 |   0.00 |   0.32 |   0.00 |   1.29
             |  75.00 |   0.00 |  25.00 |   0.00 |
             |   2.70 |   0.00 |   1.56 |   0.00 |
-------------+--------+--------+--------+--------+
promoter-TSS |     39 |     12 |      2 |     13 |     66
             |  12.62 |   3.88 |   0.65 |   4.21 |  21.36
             |  59.09 |  18.18 |   3.03 |  19.70 |
             |  35.14 |  17.91 |   3.13 |  19.40 |
-------------+--------+--------+--------+--------+
Total             111       67       64       67      309
                35.92    21.68    20.71    21.68   100.00


10cGy CHG Hyper:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     17 |    134 |    425 |    276 |    852
             |   1.50 |  11.83 |  37.51 |  24.36 |  75.20
             |   2.00 |  15.73 |  49.88 |  32.39 |
             |  10.06 |  79.76 |  93.41 |  80.94 |
-------------+--------+--------+--------+--------+
TTS          |     39 |      9 |     16 |     21 |     85
             |   3.44 |   0.79 |   1.41 |   1.85 |   7.50
             |  45.88 |  10.59 |  18.82 |  24.71 |
             |  23.08 |   5.36 |   3.52 |   6.16 |
-------------+--------+--------+--------+--------+
exon         |     74 |      9 |      0 |     15 |     98
             |   6.53 |   0.79 |   0.00 |   1.32 |   8.65
             |  75.51 |   9.18 |   0.00 |  15.31 |
             |  43.79 |   5.36 |   0.00 |   4.40 |
-------------+--------+--------+--------+--------+
intron       |      9 |      4 |      3 |      4 |     20
             |   0.79 |   0.35 |   0.26 |   0.35 |   1.77
             |  45.00 |  20.00 |  15.00 |  20.00 |
             |   5.33 |   2.38 |   0.66 |   1.17 |
-------------+--------+--------+--------+--------+
promoter-TSS |     30 |     12 |     11 |     25 |     78
             |   2.65 |   1.06 |   0.97 |   2.21 |   6.88
             |  38.46 |  15.38 |  14.10 |  32.05 |
             |  17.75 |   7.14 |   2.42 |   7.33 |
-------------+--------+--------+--------+--------+
Total             169      168      455      341     1133
                14.92    14.83    40.16    30.10   100.00




10cGy CHG Hypo:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |    158 |   1003 |   3083 |   2054 |   6298
             |   1.96 |  12.44 |  38.25 |  25.48 |  78.13
             |   2.51 |  15.93 |  48.95 |  32.61 |
             |  24.76 |  75.64 |  93.82 |  73.07 |
-------------+--------+--------+--------+--------+
TTS          |    112 |    117 |     79 |    223 |    531
             |   1.39 |   1.45 |   0.98 |   2.77 |   6.59
             |  21.09 |  22.03 |  14.88 |  42.00 |
             |  17.55 |   8.82 |   2.40 |   7.93 |
-------------+--------+--------+--------+--------+
exon         |    119 |     54 |      6 |     85 |    264
             |   1.48 |   0.67 |   0.07 |   1.05 |   3.28
             |  45.08 |  20.45 |   2.27 |  32.20 |
             |  18.65 |   4.07 |   0.18 |   3.02 |
-------------+--------+--------+--------+--------+
intron       |     48 |     18 |     23 |     63 |    152
             |   0.60 |   0.22 |   0.29 |   0.78 |   1.89
             |  31.58 |  11.84 |  15.13 |  41.45 |
             |   7.52 |   1.36 |   0.70 |   2.24 |
-------------+--------+--------+--------+--------+
promoter-TSS |    201 |    134 |     95 |    386 |    816
             |   2.49 |   1.66 |   1.18 |   4.79 |  10.12
             |  24.63 |  16.42 |  11.64 |  47.30 |
             |  31.50 |  10.11 |   2.89 |  13.73 |
-------------+--------+--------+--------+--------+
Total             638     1326     3286     2811     8061
                 7.91    16.45    40.76    34.87   100.00



100cGy CHG Hyper:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     11 |     68 |    224 |    124 |    427
             |   1.92 |  11.87 |  39.09 |  21.64 |  74.52
             |   2.58 |  15.93 |  52.46 |  29.04 |
             |  21.57 |  71.58 |  88.89 |  70.86 |
-------------+--------+--------+--------+--------+
TTS          |      9 |      7 |     11 |     15 |     42
             |   1.57 |   1.22 |   1.92 |   2.62 |   7.33
             |  21.43 |  16.67 |  26.19 |  35.71 |
             |  17.65 |   7.37 |   4.37 |   8.57 |
-------------+--------+--------+--------+--------+
exon         |     11 |      5 |      1 |     10 |     27
             |   1.92 |   0.87 |   0.17 |   1.75 |   4.71
             |  40.74 |  18.52 |   3.70 |  37.04 |
             |  21.57 |   5.26 |   0.40 |   5.71 |
-------------+--------+--------+--------+--------+
intron       |      4 |      3 |      2 |      4 |     13
             |   0.70 |   0.52 |   0.35 |   0.70 |   2.27
             |  30.77 |  23.08 |  15.38 |  30.77 |
             |   7.84 |   3.16 |   0.79 |   2.29 |
-------------+--------+--------+--------+--------+
promoter-TSS |     16 |     12 |     14 |     22 |     64
             |   2.79 |   2.09 |   2.44 |   3.84 |  11.17
             |  25.00 |  18.75 |  21.88 |  34.38 |
             |  31.37 |  12.63 |   5.56 |  12.57 |
-------------+--------+--------+--------+--------+
Total              51       95      252      175      573
                 8.90    16.58    43.98    30.54   100.00



100cGy CHG Hypo:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |      8 |     54 |    175 |    100 |    337
             |   1.89 |  12.77 |  41.37 |  23.64 |  79.67
             |   2.37 |  16.02 |  51.93 |  29.67 |
             |  24.24 |  83.08 |  94.59 |  71.43 |
-------------+--------+--------+--------+--------+
TTS          |      9 |      2 |      3 |     13 |     27
             |   2.13 |   0.47 |   0.71 |   3.07 |   6.38
             |  33.33 |   7.41 |  11.11 |  48.15 |
             |  27.27 |   3.08 |   1.62 |   9.29 |
-------------+--------+--------+--------+--------+
exon         |      7 |      5 |      1 |      9 |     22
             |   1.65 |   1.18 |   0.24 |   2.13 |   5.20
             |  31.82 |  22.73 |   4.55 |  40.91 |
             |  21.21 |   7.69 |   0.54 |   6.43 |
-------------+--------+--------+--------+--------+
intron       |      2 |      1 |      0 |      3 |      6
             |   0.47 |   0.24 |   0.00 |   0.71 |   1.42
             |  33.33 |  16.67 |   0.00 |  50.00 |
             |   6.06 |   1.54 |   0.00 |   2.14 |
-------------+--------+--------+--------+--------+
promoter-TSS |      7 |      3 |      6 |     15 |     31
             |   1.65 |   0.71 |   1.42 |   3.55 |   7.33
             |  22.58 |   9.68 |  19.35 |  48.39 |
             |  21.21 |   4.62 |   3.24 |  10.71 |
-------------+--------+--------+--------+--------+
Total              33       65      185      140      423
                 7.80    15.37    43.74    33.10   100.00



10cGy CHH Hyper:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |    155 |    784 |   2144 |   1665 |   4748
             |   2.08 |  10.53 |  28.79 |  22.35 |  63.75
             |   3.26 |  16.51 |  45.16 |  35.07 |
             |   9.38 |  69.26 |  93.54 |  70.19 |
-------------+--------+--------+--------+--------+
TTS          |    321 |    102 |     46 |    196 |    665
             |   4.31 |   1.37 |   0.62 |   2.63 |   8.93
             |  48.27 |  15.34 |   6.92 |  29.47 |
             |  19.43 |   9.01 |   2.01 |   8.26 |
-------------+--------+--------+--------+--------+
exon         |    694 |     82 |      7 |    113 |    896
             |   9.32 |   1.10 |   0.09 |   1.52 |  12.03
             |  77.46 |   9.15 |   0.78 |  12.61 |
             |  42.01 |   7.24 |   0.31 |   4.76 |
-------------+--------+--------+--------+--------+
intron       |    121 |     25 |     16 |     40 |    202
             |   1.62 |   0.34 |   0.21 |   0.54 |   2.71
             |  59.90 |  12.38 |   7.92 |  19.80 |
             |   7.32 |   2.21 |   0.70 |   1.69 |
-------------+--------+--------+--------+--------+
promoter-TSS |    361 |    139 |     79 |    358 |    937
             |   4.85 |   1.87 |   1.06 |   4.81 |  12.58
             |  38.53 |  14.83 |   8.43 |  38.21 |
             |  21.85 |  12.28 |   3.45 |  15.09 |
-------------+--------+--------+--------+--------+
Total            1652     1132     2292     2372     7448
                22.18    15.20    30.77    31.85   100.00



10cGy CHH Hypo:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |   1681 |   5770 |  12856 |  13380 |  33687
             |   3.33 |  11.45 |  25.50 |  26.54 |  66.82
             |   4.99 |  17.13 |  38.16 |  39.72 |
             |  23.48 |  66.49 |  89.60 |  66.15 |
-------------+--------+--------+--------+--------+
TTS          |   1218 |    941 |    525 |   2144 |   4828
             |   2.42 |   1.87 |   1.04 |   4.25 |   9.58
             |  25.23 |  19.49 |  10.87 |  44.41 |
             |  17.01 |  10.84 |   3.66 |  10.60 |
-------------+--------+--------+--------+--------+
exon         |   1077 |    226 |     32 |    454 |   1789
             |   2.14 |   0.45 |   0.06 |   0.90 |   3.55
             |  60.20 |  12.63 |   1.79 |  25.38 |
             |  15.04 |   2.60 |   0.22 |   2.24 |
-------------+--------+--------+--------+--------+
intron       |    831 |    245 |    126 |    431 |   1633
             |   1.65 |   0.49 |   0.25 |   0.85 |   3.24
             |  50.89 |  15.00 |   7.72 |  26.39 |
             |  11.61 |   2.82 |   0.88 |   2.13 |
-------------+--------+--------+--------+--------+
promoter-TSS |   2352 |   1496 |    810 |   3819 |   8477
             |   4.67 |   2.97 |   1.61 |   7.58 |  16.81
             |  27.75 |  17.65 |   9.56 |  45.05 |
             |  32.85 |  17.24 |   5.64 |  18.88 |
-------------+--------+--------+--------+--------+
Total            7159     8678    14349    20228    50414
                14.20    17.21    28.46    40.12   100.00


100cGy CHH Hyper:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |    738 |   3342 |   7696 |   7362 |  19138
             |   2.64 |  11.97 |  27.56 |  26.37 |  68.54
             |   3.86 |  17.46 |  40.21 |  38.47 |
             |  23.13 |  66.97 |  89.23 |  66.23 |
-------------+--------+--------+--------+--------+
TTS          |    551 |    544 |    319 |   1188 |   2602
             |   1.97 |   1.95 |   1.14 |   4.25 |   9.32
             |  21.18 |  20.91 |  12.26 |  45.66 |
             |  17.27 |  10.90 |   3.70 |  10.69 |
-------------+--------+--------+--------+--------+
exon         |    541 |    159 |     28 |    234 |    962
             |   1.94 |   0.57 |   0.10 |   0.84 |   3.45
             |  56.24 |  16.53 |   2.91 |  24.32 |
             |  16.95 |   3.19 |   0.32 |   2.11 |
-------------+--------+--------+--------+--------+
intron       |    235 |    113 |     76 |    223 |    647
             |   0.84 |   0.40 |   0.27 |   0.80 |   2.32
             |  36.32 |  17.47 |  11.75 |  34.47 |
             |   7.36 |   2.26 |   0.88 |   2.01 |
-------------+--------+--------+--------+--------+
promoter-TSS |   1126 |    832 |    506 |   2109 |   4573
             |   4.03 |   2.98 |   1.81 |   7.55 |  16.38
             |  24.62 |  18.19 |  11.06 |  46.12 |
             |  35.29 |  16.67 |   5.87 |  18.97 |
-------------+--------+--------+--------+--------+
Total            3191     4990     8625    11116    27922
                11.43    17.87    30.89    39.81   100.00



100cGy CHH Hypo:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |    592 |   1820 |   4136 |   4225 |  10773
             |   3.70 |  11.37 |  25.83 |  26.38 |  67.27
             |   5.50 |  16.89 |  38.39 |  39.22 |
             |  25.63 |  67.63 |  90.29 |  65.69 |
-------------+--------+--------+--------+--------+
TTS          |    385 |    298 |    162 |    682 |   1527
             |   2.40 |   1.86 |   1.01 |   4.26 |   9.54
             |  25.21 |  19.52 |  10.61 |  44.66 |
             |  16.67 |  11.07 |   3.54 |  10.60 |
-------------+--------+--------+--------+--------+
exon         |    361 |     77 |     12 |    147 |    597
             |   2.25 |   0.48 |   0.07 |   0.92 |   3.73
             |  60.47 |  12.90 |   2.01 |  24.62 |
             |  15.63 |   2.86 |   0.26 |   2.29 |
-------------+--------+--------+--------+--------+
intron       |    255 |     67 |     35 |    129 |    486
             |   1.59 |   0.42 |   0.22 |   0.81 |   3.03
             |  52.47 |  13.79 |   7.20 |  26.54 |
             |  11.04 |   2.49 |   0.76 |   2.01 |
-------------+--------+--------+--------+--------+
promoter-TSS |    717 |    429 |    236 |   1249 |   2631
             |   4.48 |   2.68 |   1.47 |   7.80 |  16.43
             |  27.25 |  16.31 |   8.97 |  47.47 |
             |  31.04 |  15.94 |   5.15 |  19.42 |
-------------+--------+--------+--------+--------+
Total            2310     2691     4581     6432    16014
                14.42    16.80    28.61    40.16   100.00


10cGy GC Hyper:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     82 |    432 |   1117 |    751 |   2382
             |   2.30 |  12.10 |  31.28 |  21.03 |  66.70
             |   3.44 |  18.14 |  46.89 |  31.53 |
             |   9.20 |  83.72 |  95.80 |  75.25 |
-------------+--------+--------+--------+--------+
TTS          |    142 |     28 |     18 |     63 |    251
             |   3.98 |   0.78 |   0.50 |   1.76 |   7.03
             |  56.57 |  11.16 |   7.17 |  25.10 |
             |  15.94 |   5.43 |   1.54 |   6.31 |
-------------+--------+--------+--------+--------+
exon         |    445 |     21 |      7 |     52 |    525
             |  12.46 |   0.59 |   0.20 |   1.46 |  14.70
             |  84.76 |   4.00 |   1.33 |   9.90 |
             |  49.94 |   4.07 |   0.60 |   5.21 |
-------------+--------+--------+--------+--------+
intron       |    129 |     11 |      4 |     17 |    161
             |   3.61 |   0.31 |   0.11 |   0.48 |   4.51
             |  80.12 |   6.83 |   2.48 |  10.56 |
             |  14.48 |   2.13 |   0.34 |   1.70 |
-------------+--------+--------+--------+--------+
promoter-TSS |     93 |     24 |     20 |    115 |    252
             |   2.60 |   0.67 |   0.56 |   3.22 |   7.06
             |  36.90 |   9.52 |   7.94 |  45.63 |
             |  10.44 |   4.65 |   1.72 |  11.52 |
-------------+--------+--------+--------+--------+
Total             891      516     1166      998     3571
                24.95    14.45    32.65    27.95   100.00

                      The SAS System


10cGy GC Hypo:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |      0 |      1 |      5 |      2 |      8
             |   0.00 |   5.56 |  27.78 |  11.11 |  44.44
             |   0.00 |  12.50 |  62.50 |  25.00 |
             |   0.00 |  50.00 | 100.00 |  66.67 |
-------------+--------+--------+--------+--------+
TTS          |      3 |      0 |      0 |      0 |      3
             |  16.67 |   0.00 |   0.00 |   0.00 |  16.67
             | 100.00 |   0.00 |   0.00 |   0.00 |
             |  37.50 |   0.00 |   0.00 |   0.00 |
-------------+--------+--------+--------+--------+
exon         |      5 |      1 |      0 |      0 |      6
             |  27.78 |   5.56 |   0.00 |   0.00 |  33.33
             |  83.33 |  16.67 |   0.00 |   0.00 |
             |  62.50 |  50.00 |   0.00 |   0.00 |
-------------+--------+--------+--------+--------+
promoter-TSS |      0 |      0 |      0 |      1 |      1
             |   0.00 |   0.00 |   0.00 |   5.56 |   5.56
             |   0.00 |   0.00 |   0.00 | 100.00 |
             |   0.00 |   0.00 |   0.00 |  33.33 |
-------------+--------+--------+--------+--------+
Total               8        2        5        3       18
                44.44    11.11    27.78    16.67   100.00


100cGy GC Hyper:

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |     19 |    103 |    256 |    190 |    568
             |   2.46 |  13.34 |  33.16 |  24.61 |  73.58
             |   3.35 |  18.13 |  45.07 |  33.45 |
             |  14.18 |  79.84 |  95.88 |  78.51 |
-------------+--------+--------+--------+--------+
TTS          |     16 |      5 |      7 |     11 |     39
             |   2.07 |   0.65 |   0.91 |   1.42 |   5.05
             |  41.03 |  12.82 |  17.95 |  28.21 |
             |  11.94 |   3.88 |   2.62 |   4.55 |
-------------+--------+--------+--------+--------+
exon         |     54 |     11 |      0 |     19 |     84
             |   6.99 |   1.42 |   0.00 |   2.46 |  10.88
             |  64.29 |  13.10 |   0.00 |  22.62 |
             |  40.30 |   8.53 |   0.00 |   7.85 |
-------------+--------+--------+--------+--------+
intron       |     25 |      2 |      0 |      4 |     31
             |   3.24 |   0.26 |   0.00 |   0.52 |   4.02
             |  80.65 |   6.45 |   0.00 |  12.90 |
             |  18.66 |   1.55 |   0.00 |   1.65 |
-------------+--------+--------+--------+--------+
promoter-TSS |     20 |      8 |      4 |     18 |     50
             |   2.59 |   1.04 |   0.52 |   2.33 |   6.48
             |  40.00 |  16.00 |   8.00 |  36.00 |
             |  14.93 |   6.20 |   1.50 |   7.44 |
-------------+--------+--------+--------+--------+
Total             134      129      267      242      772
                17.36    16.71    34.59    31.35   100.00


100cGy GC Hypo:


Frequency    |
Percent      |
Row Pct      |
Col Pct      |Intergen|TTS     |exon    |promoter|  Total
             |ic      |        |        |-TSS    |
-------------+--------+--------+--------+--------+
Intergenic   |      3 |     18 |     82 |     41 |    144
             |   1.44 |   8.65 |  39.42 |  19.71 |  69.23
             |   2.08 |  12.50 |  56.94 |  28.47 |
             |   7.69 |  62.07 |  97.62 |  73.21 |
-------------+--------+--------+--------+--------+
TTS          |      7 |      3 |      1 |      5 |     16
             |   3.37 |   1.44 |   0.48 |   2.40 |   7.69
             |  43.75 |  18.75 |   6.25 |  31.25 |
             |  17.95 |  10.34 |   1.19 |   8.93 |
-------------+--------+--------+--------+--------+
exon         |     18 |      3 |      0 |      3 |     24
             |   8.65 |   1.44 |   0.00 |   1.44 |  11.54
             |  75.00 |  12.50 |   0.00 |  12.50 |
             |  46.15 |  10.34 |   0.00 |   5.36 |
-------------+--------+--------+--------+--------+
intron       |      6 |      3 |      0 |      1 |     10
             |   2.88 |   1.44 |   0.00 |   0.48 |   4.81
             |  60.00 |  30.00 |   0.00 |  10.00 |
             |  15.38 |  10.34 |   0.00 |   1.79 |
-------------+--------+--------+--------+--------+
promoter-TSS |      5 |      2 |      1 |      6 |     14
             |   2.40 |   0.96 |   0.48 |   2.88 |   6.73
             |  35.71 |  14.29 |   7.14 |  42.86 |
             |  12.82 |   6.90 |   1.19 |  10.71 |
-------------+--------+--------+--------+--------+
Total              39       29       84       56      208
                18.75    13.94    40.38    26.92   100.00





*/



/* DMR vs DAR vs TE annotation */










/* PRep and export data for the following LINE PLOTS:

(1) Average GC accessibility in hyper/hypo DARs
(2) TSS accessibility plots

*/





* GC methylation data;
data gc_data;
  set arabMAP.methylation_data_gc;
  where flag_normalized=1;
  keep chr stop_pos treatment units rep total_C_norm methyl_C_norm perc_methyl_norm;
run;

proc sort data=gc_data;
  by chr stop_pos treatment units rep ;
proc means data=gc_data noprint;
  by chr stop_pos treatment units  ;
  var total_C_norm methyl_C_norm perc_methyl_norm;
  output out=gc_data2 sum(total_C_norm)=total_C sum(methyl_C_norm)=methyl_C mean(perc_methyl_norm)=perc_methyl;
run;

data gc_data3;
  set gc_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=gc_data3 out=gc_sbys10;
  where total_C >= 10;
  by chr stop_pos;
  id treatment units;
  var perc_methyl2;
run;

data gc_sbys10_2;
  set gc_sbys10;
  if _01Gy0U ne  . and _01Gy100U ne . then _01Gy_100U_0U=_01Gy100U - _01Gy0U;
  else _01Gy_100U_0U=.;

  if _1Gy0U ne . and _1Gy100U ne . then _1Gy_100U_0U=_1Gy100U - _1Gy0U;
  else _1Gy_100U_0U=.;

  if _0Gy0U ne . and _0Gy100U ne . then _0Gy_100U_0U=_0Gy100U - _0Gy0U;
  else _0Gy_100U_0U=.;

  if _01Gy_100U_0U ne . and _0Gy_100U_0U ne . then _01Gy_100U_0U_common=_01Gy_100U_0U; else _01Gy_100U_0U_common=.;
  if _1Gy_100U_0U ne . and _0Gy_100U_0U ne . then _1Gy_100U_0U_common=_1Gy_100U_0U; else _1Gy_100U_0U_common=.;
  if (_1Gy_100U_0U ne . or _01Gy_100U_0U ne .) and _0Gy_100U_0U ne . then _0Gy_100U_0U_common=_0Gy_100U_0U; else _0Gy_100U_0U_common=.;

  if _01Gy_100U_0U_common ne . and _1Gy_100U_0U_common ne .  then do;
    _01Gy_100U_0U_common_all = _01Gy_100U_0U_common;
    _1Gy_100U_0U_common_all = _1Gy_100U_0U_common;
    _0Gy_100U_0U_common_all = _0Gy_100U_0U_common;
    end;
  else do;
   _01Gy_100U_0U_common_all = .;
   _1Gy_100U_0U_common_all = .;
   _0Gy_100U_0U_common_all = . ;
    end;
  rename stop_pos=pos;
run;


/* Get DARs */


data up_dar_01_1kb up_dar_1_1kb dn_dar_01_1kb dn_dar_1_1kb;
     set arabMAP.results_by_dar_annot;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";


     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;

     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;


     keep comparison chr dar_Center plot_start plot_stop ;
run;

proc sort data=up_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;





%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=gc_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) gc_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dar_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_100U_0U_common _1Gy_100U_0U_common  _0Gy_100U_0U_common
       _01Gy_100U_0U_common_all _1Gy_100U_0U_common_all  _0Gy_100U_0U_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean( _0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean( _0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy_100U_0U_common mean_1Gy_100U_0U_common  mean_0Gy_100U_0U_common
       mean_01Gy_100U_0U_common_all mean_1Gy_100U_0U_common_all  mean_0Gy_100U_0U_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(mean_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean(mean_0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(mean_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(mean_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean(mean_0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dar_01_1kb);
%mergeMETH(dn_dar_01_1kb);
%mergeMETH(up_dar_1_1kb);
%mergeMETH(dn_dar_1_1kb);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/GC_accessibility_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hyper_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hyper_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hyper_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hyper_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hyper_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hyper_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hyper_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hyper_DAR_1kb_1Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hypo_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, hypo_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hypo_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common_all, hypo_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hypo_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common, hypo_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hypo_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_100U_0U_common, mean_1Gy_100U_0U_common_all, hypo_DAR_1kb_1Gy_0Gy_common_binned);

/* As above, but now for TSSs */


data site2promoter;
  set arabMAP.results_by_dac_annot;
  length gene_id $20.;
  if abs(distance_to_tss) > 999 then delete;
  if count(nearest_promoterID, "-T1") > 0 then gene_ID=compress(upcase(tranwrd(Nearest_PromoterID,"-T1","")));
  else gene_ID=compress(upcase(scan(Nearest_PromoterID,1,".")));
  keep gene_ID chr start_pos stop_pos strand distance_to_tss nearest_promoterID;
  rename stop_pos=pos nearest_promoterID=transcript_id;
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
  keep gene_id VAR7; 
  rename VAR7=strand;
run;

proc sort data=gene nodup;
  by gene_id;
proc sort data=site2promoter nodup;
    by gene_id;
run;

data site2promoter2 no_strand no_site;
  merge site2promoter (in=in1) gene (in=in2);
  by gene_id;
  if in1 and in2 then output site2promoter2;
  else if in1 then output no_strand;
  else output no_site;
run;

/* distance to TSS is relative to strand of transcript, so I don't need to flip anything!!! */

proc sort data=site2promoter2 nodup;
  by chr pos;
proc sort data=gc_sbys10_2;
  by chr pos;
run;

data site2promoter2_w_meth;
  merge site2promoter2 (in=in1) gc_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data site2promoter2_w_meth2;
  set site2promoter2_w_meth;
  grouped_pos=int(distance_to_TSS/10) * 10;
  if distance_to_TSS < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=site2promoter2_w_meth2;
  by distance_to_TSS grouped_pos2 ;
run;


proc means data=site2promoter2_w_meth2 noprint;
  by distance_to_TSS grouped_pos2  ;
  var  _01Gy_100U_0U_common _1Gy_100U_0U_common  _0Gy_100U_0U_common
       _01Gy_100U_0U_common_all _1Gy_100U_0U_common_all  _0Gy_100U_0U_common_all ;
  output out=mean_diff_tss
  mean(_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean( _0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean( _0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;

proc sort data=mean_diff_tss ;
  by grouped_pos2;
run;


proc means data=mean_diff_tss noprint;
  by grouped_pos2  ;
  var  mean_01Gy_100U_0U_common mean_1Gy_100U_0U_common  mean_0Gy_100U_0U_common
       mean_01Gy_100U_0U_common_all mean_1Gy_100U_0U_common_all  mean_0Gy_100U_0U_common_all ;
  output out=mean_diff_tss_2
  mean(mean_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(mean_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean(mean_0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(mean_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(mean_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean(mean_0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;


data mean_diff_tss_1;
  set mean_diff_tss;
  drop _TYPE_ _FREQ_;
  keep distance_To_tss mean_: ;
  rename distance_to_tss=pos;
run;

data mean_diff_tss_3;
  set mean_diff_tss_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;


/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, var3, outName);

data export;
  retain pos &var1. &var2. &var3.;
  set &inData.;
  keep pos &var1. &var2. &var3.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/GC_accessibility_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_tss_1, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, mean_1Gy_100U_0U_common, GC_acc_TSS_1kb);
%exportLine( mean_diff_tss_3, mean_0Gy_100U_0U_common, mean_01Gy_100U_0U_common, mean_1Gy_100U_0U_common, GC_acc_TSS_1kb_binned);

%exportLine( mean_diff_tss_1, mean_0Gy_100U_0U_common_all, mean_01Gy_100U_0U_common_all, mean_1Gy_100U_0U_common_all, GC_acc_TSS_1kb_common);
%exportLine( mean_diff_tss_3, mean_0Gy_100U_0U_common_all, mean_01Gy_100U_0U_common_all, mean_1Gy_100U_0U_common_all, GC_acc_TSS_1kb_common_binned);



/*************************************************************************************/


/* Methylation plots:
   same as GC plots (se we can reuse the same code!)
   but also need to add 
   methylation in DMRs by site type (i.e. CG for CG DMRs, CHG for CHG DMRs, CHH for CHH DMRs) not just DARs (but code is good)

   Make a big-ass macro to do everything for ONE site type at a time, then kick off CG, CHG, CHH

*/


%macro methDataGen(siteType);


/* DMR counts:
   Hyper all, exon, intron, promoter, downstream; 0.1G, 1G
   Hypo all, exon, intron, promoter, downstream; 0.1G, 1G

 */

data up_dmr_01_72 up_dmr_1_72 dn_dmr_01_72 dn_dmr_1_72;
     set arabMAP.results_by_dmr_annot;
     length feature $32.;
     where site_type="&siteType.";
     feature = scan(annotation, 1, " ");
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72;

     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
run;

proc sort data=up_dmr_01_72 nodup; by _all_; run;
proc sort data=dn_dmr_01_72 nodup; by _all_; run;
proc sort data=up_dmr_1_72 nodup; by _all_; run;
proc sort data=dn_dmr_1_72 nodup; by _all_; run;


proc freq data=up_dmr_01_72; tables feature; run;
proc freq data=dn_dmr_01_72; tables feature; run;
proc freq data=up_dmr_1_72; tables feature; run;
proc freq data=dn_dmr_1_72; tables feature; run;

/* 

CG sites:
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   5           9       10      0       3           27
0.1Gy   Hypo    21          20      11      6       54          112
1Gy     Hyper   0           1       1       0       0           2
1Gy     Hypo    2           4       1       1       12          20



CHG sites:
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   9           12      13      0       3           37
0.1Gy   Hypo    258         29      4       9       25          325
1Gy     Hyper   4           2       1       0       0           7
1Gy     Hypo    2           0       0       0       1           3


CHH sites:
Dose    up/down intergenic  TTS     exon    intron  promoter    Total
0.1Gy   Hyper   737         119     126     31      139         1152
0.1Gy   Hypo    16996       2510    531     603     4705        25345
1Gy     Hyper   7203        978     271     229     1762        10443
1Gy     Hypo    2372        386     96      93      617         3564


*/

/* count by gene */

data up_dmr_01_72_gn up_dmr_1_72_gn dn_dmr_01_72_gn dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="&siteType.";
     if geneID = "" then delete;
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_72_gn;
     keep geneID;
run;

proc sort  data=up_dmr_01_72_gn nodup; by geneID; run;
proc sort data=dn_dmr_01_72_gn nodup; by geneID; run;
proc sort  data=up_dmr_1_72_gn nodup; by geneID; run;
proc sort  data=dn_dmr_1_72_gn nodup; by geneID; run;

data up_01_1_gn;
  merge up_dmr_01_72_gn (in=in1) up_dmr_1_72_gn (in=in2);
  by geneID;
  if in1 then dmr_up_01=1; else dmr_up_01=0;
  if in2 then dmr_up_1=1; else dmr_up_1=0;
run;

data dn_01_1_gn;
  merge dn_dmr_01_72_gn (in=in1) dn_dmr_1_72_gn (in=in2);
  by geneID;
  if in1 then dmr_dn_01=1; else dmr_dn_01=0;
  if in2 then dmr_dn_1=1; else dmr_dn_1=0;
run;


data up_dn_01_gn;
  merge up_dmr_01_72_gn (in=in1) dn_dmr_01_72_gn (in=in2);
  by geneID;
  if in1 then dmr_up_01=1; else dmr_up_01=0;
  if in2 then dmr_dn_01=1; else dmr_dn_01=0;
run;

data up_dn_1_gn;
  merge up_dmr_1_72_gn (in=in1) dn_dmr_1_72_gn (in=in2);
  by geneID;
  if in1 then dmr_up_1=1; else dmr_up_1=0;
  if in2 then dmr_dn_1=1; else dmr_dn_1=0;
run;

proc freq data=up_01_1_gn; tables dmr_up_01*dmr_up_1; run;
proc freq data=dn_01_1_gn; tables dmr_dn_01*dmr_dn_1; run;
proc freq data=up_dn_01_gn; tables dmr_up_01*dmr_dn_01; run;
proc freq data=up_dn_1_gn; tables dmr_up_1*dmr_dn_1; run;


/* 
CG:
0.1Gy hyper-DMG 25
0.1Gy hypo-DMG  92
1Gy hyper-DMG   2
1Gy hypo-DMG    20


       Table of dmr_up_01 by dmr_up_1

    dmr_up_01     dmr_up_1

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |      0 |      1 |      1
             |   0.00 |   3.85 |   3.85
             |   0.00 | 100.00 |
             |   0.00 |  50.00 |
    ---------+--------+--------+
           1 |     24 |      1 |     25
             |  92.31 |   3.85 |  96.15
             |  96.00 |   4.00 |
             | 100.00 |  50.00 |
    ---------+--------+--------+
    Total          24        2       26
                92.31     7.69   100.00


    Table of dmr_dn_01 by dmr_dn_1

 dmr_dn_01     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |      6 |      6
          |   0.00 |   6.12 |   6.12
          |   0.00 | 100.00 |
          |   0.00 |  30.00 |
 ---------+--------+--------+
        1 |     78 |     14 |     92
          |  79.59 |  14.29 |  93.88
          |  84.78 |  15.22 |
          | 100.00 |  70.00 |
 ---------+--------+--------+
 Total          78       20       98
             79.59    20.41   100.00

            The SAS System


   dmr_up_01     dmr_dn_01

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |     91 |     91
            |   0.00 |  78.45 |  78.45
            |   0.00 | 100.00 |
            |   0.00 |  98.91 |
   ---------+--------+--------+
          1 |     24 |      1 |     25
            |  20.69 |   0.86 |  21.55
            |  96.00 |   4.00 |
            | 100.00 |   1.09 |
   ---------+--------+--------+
   Total          24       92      116
               20.69    79.31   100.00



 dmr_up_1     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |     20 |     20
          |   0.00 |  90.91 |  90.91
          |   0.00 | 100.00 |
          |   0.00 | 100.00 |
 ---------+--------+--------+
        1 |      2 |      0 |      2
          |   9.09 |   0.00 |   9.09
          | 100.00 |   0.00 |
          | 100.00 |   0.00 |
 ---------+--------+--------+
 Total           2       20       22
              9.09    90.91   100.00




CHG:
0.1Gy hyper-DMG 
0.1Gy hypo-DMG  
1Gy hyper-DMG   
1Gy hypo-DMG    


     Table of dmr_up_01 by dmr_up_1

  dmr_up_01     dmr_up_1

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      6 |      6
           |   0.00 |  14.29 |  14.29
           |   0.00 | 100.00 |
           |   0.00 |  85.71 |
  ---------+--------+--------+
         1 |     35 |      1 |     36
           |  83.33 |   2.38 |  85.71
           |  97.22 |   2.78 |
           | 100.00 |  14.29 |
  ---------+--------+--------+
  Total          35        7       42
              83.33    16.67   100.00

             The SAS System


   Table of dmr_dn_01 by dmr_dn_1

dmr_dn_01     dmr_dn_1

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |      0 |      1 |      1
         |   0.00 |   0.37 |   0.37
         |   0.00 | 100.00 |
         |   0.00 |  33.33 |
---------+--------+--------+
       1 |    269 |      2 |    271
         |  98.90 |   0.74 |  99.63
         |  99.26 |   0.74 |
         | 100.00 |  66.67 |
---------+--------+--------+
Total         269        3      272
            98.90     1.10   100.00

   Table of dmr_up_01 by dmr_dn_01

 dmr_up_01     dmr_dn_01

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    266 |    266
          |   0.00 |  88.08 |  88.08
          |   0.00 | 100.00 |
          |   0.00 |  98.15 |
 ---------+--------+--------+
        1 |     31 |      5 |     36
          |  10.26 |   1.66 |  11.92
          |  86.11 |  13.89 |
          | 100.00 |   1.85 |
 ---------+--------+--------+
 Total          31      271      302
             10.26    89.74   100.00

     Table of dmr_up_1 by dmr_dn_1

  dmr_up_1     dmr_dn_1

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      3 |      3
           |   0.00 |  30.00 |  30.00
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |      7 |      0 |      7
           |  70.00 |   0.00 |  70.00
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total           7        3       10
              70.00    30.00   100.00

             The SAS System



CHH:
0.1Gy hyper-DMG 
0.1Gy hypo-DMG  
1Gy hyper-DMG   
1Gy hypo-DMG    

      Table of dmr_up_01 by dmr_up_1

   dmr_up_01     dmr_up_1

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |   3134 |   3134
            |   0.00 |  79.16 |  79.16
            |   0.00 | 100.00 |
            |   0.00 |  83.71 |
   ---------+--------+--------+
          1 |    215 |    610 |    825
            |   5.43 |  15.41 |  20.84
            |  26.06 |  73.94 |
            | 100.00 |  16.29 |
   ---------+--------+--------+
   Total         215     3744     3959
                5.43    94.57   100.00

              The SAS System

            The FREQ Procedure

    Table of dmr_dn_01 by dmr_dn_1

 dmr_dn_01     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    119 |    119
          |   0.00 |   1.84 |   1.84
          |   0.00 | 100.00 |
          |   0.00 |   5.62 |
 ---------+--------+--------+
        1 |   4366 |   1997 |   6363
          |  67.36 |  30.81 |  98.16
          |  68.62 |  31.38 |
          | 100.00 |  94.38 |
 ---------+--------+--------+
 Total        4366     2116     6482
             67.36    32.64   100.00

    Table of dmr_up_01 by dmr_dn_01

  dmr_up_01     dmr_dn_01

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   5735 |   5735
           |   0.00 |  87.42 |  87.42
           |   0.00 | 100.00 |
           |   0.00 |  90.13 |
  ---------+--------+--------+
         1 |    197 |    628 |    825
           |   3.00 |   9.57 |  12.58
           |  23.88 |  76.12 |
           | 100.00 |   9.87 |
  ---------+--------+--------+
  Total         197     6363     6560
               3.00    97.00   100.00

    Table of dmr_up_1 by dmr_dn_1

 dmr_up_1     dmr_dn_1

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    794 |    794
          |   0.00 |  17.50 |  17.50
          |   0.00 | 100.00 |
          |   0.00 |  37.52 |
 ---------+--------+--------+
        1 |   2422 |   1322 |   3744
          |  53.37 |  29.13 |  82.50
          |  64.69 |  35.31 |
          | 100.00 |  62.48 |
 ---------+--------+--------+
 Total        2422     2116     4538
             53.37    46.63   100.00

*/



data dmr_up_dn_all_gn;
  merge  up_dmr_01_72_gn (in=in1) up_dmr_1_72_gn (in=in2) dn_dmr_01_72_gn (in=in3) dn_dmr_1_72_gn (in=in4);
  by geneID;
  if in1 then dmr_up_01=1; else dmr_up_01=0;
  if in1 then dmr_up_1=1; else dmr_up_1=0;
  if in2 then dmr_dn_01=1; else dmr_dn_01=0;
  if in2 then dmr_dn_1=1; else dmr_dn_1=0;
run;

proc freq data=dmr_up_dn_all_gn noprint;
  tables dmr_up_01*dmr_up_1*dmr_dn_01*dmr_dn_1 / out=dag_compare;
proc print data=dag_compare;
run;

/*
CG:

    dmr_                 dmr_
   up_01    dmr_up_1    dn_01    dmr_dn_1    COUNT

     0          0         0          0         97
     0          0         1          1          1
     1          1         0          0         24
     1          1         1          1          1





CHG:

   dmr_                 dmr_
  up_01    dmr_up_1    dn_01    dmr_dn_1    COUNT

    0          0         0          0        267
    0          0         1          1          6
    1          1         0          0         35
    1          1         1          1          1



CHH:

     dmr_                 dmr_
    up_01    dmr_up_1    dn_01    dmr_dn_1    COUNT

      0          0         0          0        3268
      0          0         1          1        3134
      1          1         0          0         215
      1          1         1          1         610



*/



/* overlapping DMRs now

need to merge region here!!

0.1Gy vs 1Gy hyper (DMR)
0.1Gy vs 1Gy hypo (DMR)
0.1Gy hyper vs hypo (DMR)
0.1Gy hyper vs hypo (DMR)

*/



%macro commonDMR(dataA, dataB, outName);

data stack_&outName.;
  set &dataA. (in=in1) &dataB. (in=in2);
  length comp $20.;
  if in1 then comp="&dataA.";
  if in2 then comp="&dataB.";
  keep comp site_type chr  region_start region_stop ;
run;

proc sort data=stack_&outName. nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_&outName._make_super;
  retain superregion_num;
  set stack_&outName.;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_&outName._super1;
   set stack_&outName._make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_&outName._super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_&outName._super1 out=stack_&outName._super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;

data stack_&outName._super_sbys2;
  set stack_&outName._super_sbys;
  if &dataA.=. then &dataA.=0;
  if &dataB.=. then &dataB.=0;
run;

proc freq data=stack_&outName._super_sbys2;
  tables &dataA.*&dataB. ;
run;
%mend;

%commonDMR(up_dmr_01_72, up_dmr_1_72, up_dmr_72);
%commonDMR(dn_dmr_01_72, dn_dmr_1_72, dn_dmr_72);
%commonDMR(up_dmr_01_72, dn_dmr_01_72, updn_01_dmr_72);
%commonDMR(up_dmr_1_72, dn_dmr_1_72, updn_1_dmr_72);


/*


CG:
   Table of up_dmr_01_72 by up_dmr_1_72

   up_dmr_01_72     up_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |      1 |      1
            |   0.00 |   3.57 |   3.57
            |   0.00 | 100.00 |
            |   0.00 |  50.00 |
   ---------+--------+--------+
          1 |     26 |      1 |     27
            |  92.86 |   3.57 |  96.43
            |  96.30 |   3.70 |
            | 100.00 |  50.00 |
   ---------+--------+--------+
   Total          26        2       28
               92.86     7.14   100.00

              The SAS System


  dn_dmr_01_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      6 |      6
           |   0.00 |   5.17 |   5.17
           |   0.00 | 100.00 |
           |   0.00 |  30.00 |
  ---------+--------+--------+
         1 |     96 |     14 |    110
           |  82.76 |  12.07 |  94.83
           |  87.27 |  12.73 |
           | 100.00 |  70.00 |
  ---------+--------+--------+
  Total          96       20      116
              82.76    17.24   100.00


Table of up_dmr_01_72 by dn_dmr_01_72

 up_dmr_01_72     dn_dmr_01_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    112 |    112
          |   0.00 |  80.58 |  80.58
          |   0.00 | 100.00 |
          |   0.00 | 100.00 |
 ---------+--------+--------+
        1 |     27 |      0 |     27
          |  19.42 |   0.00 |  19.42
          | 100.00 |   0.00 |
          | 100.00 |   0.00 |
 ---------+--------+--------+
 Total          27      112      139
             19.42    80.58   100.00


  Table of up_dmr_1_72 by dn_dmr_1_72

  up_dmr_1_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |     20 |     20
           |   0.00 |  90.91 |  90.91
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |      2 |      0 |      2
           |   9.09 |   0.00 |   9.09
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total           2       20       22
               9.09    90.91   100.00



CHG:
   Table of up_dmr_01_72 by up_dmr_1_72

   up_dmr_01_72     up_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |      6 |      6
            |   0.00 |  13.95 |  13.95
            |   0.00 | 100.00 |
            |   0.00 |  85.71 |
   ---------+--------+--------+
          1 |     36 |      1 |     37
            |  83.72 |   2.33 |  86.05
            |  97.30 |   2.70 |
            | 100.00 |  14.29 |
   ---------+--------+--------+
   Total          36        7       43
               83.72    16.28   100.00

 Table of dn_dmr_01_72 by dn_dmr_1_72

 dn_dmr_01_72     dn_dmr_1_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |      2 |      2
          |   0.00 |   0.61 |   0.61
          |   0.00 | 100.00 |
          |   0.00 |  66.67 |
 ---------+--------+--------+
        1 |    324 |      1 |    325
          |  99.08 |   0.31 |  99.39
          |  99.69 |   0.31 |
          | 100.00 |  33.33 |
 ---------+--------+--------+
 Total         324        3      327
             99.08     0.92   100.00

            The SAS System

 Table of up_dmr_01_72 by dn_dmr_01_72

  up_dmr_01_72     dn_dmr_01_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    325 |    325
           |   0.00 |  89.78 |  89.78
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |     37 |      0 |     37
           |  10.22 |   0.00 |  10.22
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total          37      325      362
              10.22    89.78   100.00

  Table of up_dmr_1_72 by dn_dmr_1_72

  up_dmr_1_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |      3 |      3
           |   0.00 |  30.00 |  30.00
           |   0.00 | 100.00 |
           |   0.00 | 100.00 |
  ---------+--------+--------+
         1 |      7 |      0 |      7
           |  70.00 |   0.00 |  70.00
           | 100.00 |   0.00 |
           | 100.00 |   0.00 |
  ---------+--------+--------+
  Total           7        3       10
              70.00    30.00   100.00

CHH:



   Table of up_dmr_01_72 by up_dmr_1_72

   up_dmr_01_72     up_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |  10118 |  10118
            |   0.00 |  89.80 |  89.80
            |   0.00 | 100.00 |
            |   0.00 |  96.96 |
   ---------+--------+--------+
          1 |    832 |    317 |   1149
            |   7.38 |   2.81 |  10.20
            |  72.41 |  27.59 |
            | 100.00 |   3.04 |
   ---------+--------+--------+
   Total         832    10435    11267
                7.38    92.62   100.00




  dn_dmr_01_72     dn_dmr_1_72

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    876 |    876
           |   0.00 |   3.37 |   3.37
           |   0.00 | 100.00 |
           |   0.00 |  24.85 |
  ---------+--------+--------+
         1 |  22439 |   2649 |  25088
           |  86.42 |  10.20 |  96.63
           |  89.44 |  10.56 |
           | 100.00 |  75.15 |
  ---------+--------+--------+
  Total       22439     3525    25964
              86.42    13.58   100.00


Table of up_dmr_01_72 by dn_dmr_01_72

 up_dmr_01_72     dn_dmr_01_72

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |  25345 |  25345
          |   0.00 |  95.65 |  95.65
          |   0.00 | 100.00 |
          |   0.00 | 100.00 |
 ---------+--------+--------+
        1 |   1152 |      0 |   1152
          |   4.35 |   0.00 |   4.35
          | 100.00 |   0.00 |
          | 100.00 |   0.00 |
 ---------+--------+--------+
 Total        1152    25345    26497
              4.35    95.65   100.00

 Table of up_dmr_1_72 by dn_dmr_1_72

    up_dmr_1_72     dn_dmr_1_72

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |      0 |   3564 |   3564
            |   0.00 |  25.44 |  25.44
            |   0.00 | 100.00 |
            |   0.00 | 100.00 |
   ---------+--------+--------+
          1 |  10443 |      0 |  10443
            |  74.56 |   0.00 |  74.56
            | 100.00 |   0.00 |
            | 100.00 |   0.00 |
   ---------+--------+--------+
   Total       10443     3564    14007
               74.56    25.44   100.00

*/

/* 4 way venn of regions */


data stack_all_dmr;
  set up_dmr_01_72 (in=in1) dn_dmr_01_72 (in=in2) up_dmr_1_72 (in=in3) dn_dmr_1_72 (in=in4);
  length comp $20.;
  if in1 then comp="up_dmr_01_72";
  if in2 then comp="dn_dmr_01_72";
  if in3 then comp="up_dmr_1_72";
  if in4 then comp="dn_dmr_1_72";
  keep comp site_type chr  region_start region_stop ;
run;

proc sort data=stack_all_dmr nodup;
  by  site_type chr  region_start region_stop comp;
run;


data stack_all_dmr_make_super;
  retain superregion_num;
  set stack_all_dmr;
  by site_type chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_all_dmr_super1;
   set stack_all_dmr_make_super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_all_dmr_super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_all_dmr_super1 out=stack_all_dmr_super_sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;


data stack_all_dmr_super_sbys2;
  set stack_all_dmr_super_sbys;
  if up_dmr_01_72=. then up_dmr_01_72=0;
  if dn_dmr_01_72=. then dn_dmr_01_72=0;
  if up_dmr_1_72=. then up_dmr_1_72=0;
  if dn_dmr_1_72=. then dn_dmr_1_72=0;
run;

proc freq data=stack_all_dmr_super_sbys2 noprint;
  tables up_dmr_01_72*up_dmr_1_72*dn_dmr_01_72*dn_dmr_1_72 / out=dmr_compare;
run;
proc print data=dmr_compare;
run;


/*
CG:

up_dmr_    up_dmr_    dn_dmr_    dn_dmr_
 01_72       1_72      01_72       1_72     COUNT

   0          0          0          1          6
   0          0          1          0         96
   0          0          1          1         14 <- hypo
   0          1          0          0          1
   1          0          0          0         26
   1          1          0          0          1 <- hyper

                    The SAS System

CHG:


 up_dmr_    up_dmr_    dn_dmr_    dn_dmr_
  01_72       1_72      01_72       1_72     COUNT

    0          0          0          1          2
    0          0          1          0        324
    0          0          1          1          1 <- hypo
    0          1          0          0          6
    1          0          0          0         36
    1          1          0          0          1 <- hyper


CHH:


 up_dmr_    up_dmr_    dn_dmr_    dn_dmr_
  01_72       1_72      01_72       1_72     COUNT

    0          0          0          1         876
    0          0          1          0       19628
    0          0          1          1        2532   <- hypomethylated-CHH
    0          1          0          0        7300
    0          1          1          0        2585
    0          1          1          1          96
    1          0          0          0         814
    1          0          0          1          14
    1          0          1          1           5
    1          1          0          0         295   <- hypermethylated-CHH
    1          1          1          0          18
    1          1          1          1           2



 */


/* PRep and export data for the following LINE PLOTS:

(1) Average accessibility in hyper/hypo DARs
(1) Average accessibility in hyper/hypo DMRs
(2) TSS accessibility plots

*/


/* methylation in DARs */


* GC methylation data;
data meth_data;
  set arabMAP.methylation_data_cg_chg_chh;
  where site_type="&siteType.";
  keep chr stop_pos treatment rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by chr stop_pos treatment  rep ;
proc means data=meth_data noprint;
  by chr stop_pos treatment   ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;

data meth_data3;
  set meth_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=meth_data3 out=meth_sbys10;
  where total_C >= 10;
  by chr stop_pos;
  id treatment ;
  var perc_methyl2;
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


/* Get DARs */


data up_dar_01_1kb up_dar_1_1kb dn_dar_01_1kb dn_dar_1_1kb;
     set arabMAP.results_by_dar_annot;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";

     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;

     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_1kb;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_1kb;


     keep comparison chr dar_Center plot_start plot_stop ;
run;

proc sort data=up_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_01_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_1_1kb nodup;  by chr dar_center plot_start plot_stop; run;





%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=meth_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dar_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dar_01_1kb);
%mergeMETH(dn_dar_01_1kb);
%mergeMETH(up_dar_1_1kb);
%mergeMETH(dn_dar_1_1kb);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hyper_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hyper_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hyper_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_up_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hyper_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hyper_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hyper_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_up_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hyper_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_up_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hyper_DAR_1kb_1Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hypo_DAR_1kb_01Gy_0Gy);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hypo_DAR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hypo_DAR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hypo_DAR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hypo_DAR_1kb_1Gy_0Gy);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hypo_DAR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_dn_dar_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hypo_DAR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_dn_dar_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hypo_DAR_1kb_1Gy_0Gy_common_binned);



/* Do the same but for DMRs */

/* Get DMRs */




data up_dmr_01_1kb up_dmr_1_1kb dn_dmr_01_1kb dn_dmr_1_1kb;
     set arabMAP.results_by_dmr_annot;
     length feature $32.;
     where site_type="&siteType.";
     feature = scan(annotation, 1, " ");

     dmr_center=int((region_start + region_stop) / 2);
     plot_start=dmr_center - 999;
     plot_stop=dmr_center + 999;

     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output up_dmr_01_1kb;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dn_dmr_01_1kb;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output up_dmr_1_1kb;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_1kb;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output up_dmr_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output up_dmr_1_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dn_dmr_01_1kb;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dn_dmr_1_1kb;

     keep comparison chr dmr_Center plot_start plot_stop ;
run;

proc sort data=up_dmr_01_1kb nodup;  by chr dmr_center plot_start plot_stop; run;
proc sort data=dn_dmr_01_1kb nodup;  by chr dmr_center plot_start plot_stop; run;
proc sort data=up_dmr_1_1kb nodup;  by chr dmr_center plot_start plot_stop; run;
proc sort data=dn_dmr_1_1kb nodup;  by chr dmr_center plot_start plot_stop; run;





%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dmr_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=meth_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dmr_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dmr_01_1kb);
%mergeMETH(dn_dmr_01_1kb);
%mergeMETH(up_dmr_1_1kb);
%mergeMETH(dn_dmr_1_1kb);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hyper_DMR_1kb_01Gy_0Gy);
%exportLine( mean_diff_up_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hyper_DMR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_up_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hyper_DMR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_up_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hyper_DMR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_up_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hyper_DMR_1kb_1Gy_0Gy);
%exportLine( mean_diff_up_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hyper_DMR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_up_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hyper_DMR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_up_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hyper_DMR_1kb_1Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common, hypo_DMR_1kb_01Gy_0Gy);
%exportLine( mean_diff_dn_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common, hypo_DMR_1kb_01Gy_0Gy_binned);
%exportLine( mean_diff_dn_dmr_01_1kb_1, mean_0Gy_common, mean_01Gy_common_all, hypo_DMR_1kb_01Gy_0Gy_common);
%exportLine( mean_diff_dn_dmr_01_1kb_3, mean_0Gy_common, mean_01Gy_common_all, hypo_DMR_1kb_01Gy_0Gy_common_binned);

%exportLine( mean_diff_dn_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common, hypo_DMR_1kb_1Gy_0Gy);
%exportLine( mean_diff_dn_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common, hypo_DMR_1kb_1Gy_0Gy_binned);
%exportLine( mean_diff_dn_dmr_1_1kb_1, mean_0Gy_common, mean_1Gy_common_all, hypo_DMR_1kb_1Gy_0Gy_common);
%exportLine( mean_diff_dn_dmr_1_1kb_3, mean_0Gy_common, mean_1Gy_common_all, hypo_DMR_1kb_1Gy_0Gy_common_binned);


/* As above, but now for TSSs */


data site2promoter;
  set arabMAP.results_by_dmc_annot;
  where site_type="&siteType.";
  length gene_id $20.;
  if abs(distance_to_tss) > 999 then delete;
  if count(nearest_promoterID, "-T1") > 0 then gene_ID=compress(upcase(tranwrd(Nearest_PromoterID,"-T1","")));
  else gene_ID=compress(upcase(scan(Nearest_PromoterID,1,".")));
  keep gene_ID chr start_pos stop_pos strand distance_to_tss nearest_promoterID;
  rename stop_pos=pos nearest_promoterID=transcript_id;
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
  keep gene_id VAR7; 
  rename VAR7=strand;
run;

proc sort data=gene nodup;
  by gene_id;
proc sort data=site2promoter nodup;
    by gene_id;
run;

data site2promoter2 no_strand no_site;
  merge site2promoter (in=in1) gene (in=in2);
  by gene_id;
  if in1 and in2 then output site2promoter2;
  else if in1 then output no_strand;
  else output no_site;
run;

/* distance to TSS is relative to strand of transcript, so I don't need to flip anything!!! */

proc sort data=site2promoter2 nodup;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data site2promoter2_w_meth;
  merge site2promoter2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data site2promoter2_w_meth2;
  set site2promoter2_w_meth;
  grouped_pos=int(distance_to_TSS/10) * 10;
  if distance_to_TSS < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=site2promoter2_w_meth2;
  by distance_to_TSS grouped_pos2 ;
run;


proc means data=site2promoter2_w_meth2 noprint;
  by distance_to_TSS grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_tss
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_tss ;
  by grouped_pos2;
run;


proc means data=mean_diff_tss noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_tss_2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_tss_1;
  set mean_diff_tss;
  drop _TYPE_ _FREQ_;
  keep distance_To_tss mean_: ;
  rename distance_to_tss=pos;
run;

data mean_diff_tss_3;
  set mean_diff_tss_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;


/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, var3, outName);

data export;
  retain pos &var1. &var2. &var3.;
  set &inData.;
  keep pos &var1. &var2. &var3.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_tss_1, mean_0Gy_common, mean_01Gy_common, mean_1Gy_common, TSS_1kb);
%exportLine( mean_diff_tss_3, mean_0Gy_common, mean_01Gy_common, mean_1Gy_common, TSS_1kb_binned);

%exportLine( mean_diff_tss_1, mean_0Gy_common_all, mean_01Gy_common_all, mean_1Gy_common_all, TSS_1kb_common);
%exportLine( mean_diff_tss_3, mean_0Gy_common_all, mean_01Gy_common_all, mean_1Gy_common_all, TSS_1kb_common_binned);


%mend;

%methDataGen(CG);
%methDataGen(CHG);
%methDataGen(CHH);



/* scatterplots :
  (1) rep to rep concordance count shared/unique, pearson correlation on perc methyl (all, shared)
  (2) HCG (mean) shared only, pearson correlation
  */

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
  else do; flag_rep2_&treatment._&units.=0; perc_methyl_&treatment._&units._2=0; total_c_1=0; end;
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

data export_Data;
   set meth_data_compare2;
   keep chr pos perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_&siteType._&treatment._&units..csv"
     dbms=csv replace;
run;

%mend;


%concordance(CG, 0Gy, 0U);
/*
      flag_rep1_0Gy_0U
             flag_rep2_0Gy_0U

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       1|  Total
   ---------+--------+
          0 |   1128 |   1128
            |   0.03 |   0.03
            | 100.00 |
            |   0.03 |
   ---------+--------+
          1 |3924221 |3924221
            |  99.97 |  99.97
            | 100.00 |
            |  99.97 |
   ---------+--------+
   Total     3925349  3925349
              100.00   100.00

         The SAS System
  
Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

perc_methyl_0Gy_0U_1     5154586       0.30002       0.40726       1546473             0       1.00000
perc_methyl_0Gy_0U_2     5154586       0.30033       0.40696       1548099             0       1.00000


                            Pearson Correlation Coefficients, N = 5154586
                                      Prob > |r| under H0: Rho=0

                                                 perc_methyl_      perc_methyl_
                                                     0Gy_0U_1          0Gy_0U_2

                       perc_methyl_0Gy_0U_1           1.00000           0.94429
                                                                         <.0001

                       perc_methyl_0Gy_0U_2           0.94429           1.00000
                                                       <.0001                            

   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_0Gy_0U_1     3924221       0.33743       0.41563       1324132             0       1.00000
   perc_methyl_0Gy_0U_2     3924221       0.33727       0.41543       1323525             0       1.00000


                               Pearson Correlation Coefficients, N = 3924221
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        0Gy_0U_1          0Gy_0U_2

                          perc_methyl_0Gy_0U_1           1.00000           0.96658
                                                                            <.0001

                          perc_methyl_0Gy_0U_2           0.96658           1.00000
                                                          <.0001


*/
%concordance(CG, 01Gy, 0U);
/*
 flag_rep1_01Gy_0U
           flag_rep2_01Gy_0U

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        0 |    331 |    331
          |   0.01 |   0.01
          | 100.00 |
          |   0.01 |
 ---------+--------+
        1 |4489853 |4489853
          |  99.99 |  99.99
          | 100.00 |
          |  99.99 |
 ---------+--------+
 Total     4490184  4490184
            100.00   100.00



    Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_01Gy_0U_1     5316165       0.29626       0.39809       1574957             0       1.00000
    perc_methyl_01Gy_0U_2     5316165       0.29402       0.39736       1563068             0       1.00000


                                Pearson Correlation Coefficients, N = 5316165
                                          Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                         01Gy_0U_1         01Gy_0U_2

                           perc_methyl_01Gy_0U_1           1.00000           0.94806
                                                                              <.0001

                           perc_methyl_01Gy_0U_2           0.94806           1.00000
                                                            <.0001



   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1     4489853       0.31312       0.40114       1405852             0       1.00000
   perc_methyl_01Gy_0U_2     4489853       0.31086       0.40033       1395720             0       1.00000


                               Pearson Correlation Coefficients, N = 4489853
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.96354
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.96354           1.00000
                                                           <.0001
*/
%concordance(CG, 1Gy, 0U);
/*

    flag_rep1_1Gy_0U
              flag_rep2_1Gy_0U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           0 |   7961 |   7961
             |   0.22 |   0.22
             | 100.00 |
             |   0.22 |
    ---------+--------+
           1 |3560334 |3560334
             |  99.78 |  99.78
             | 100.00 |
             |  99.78 |
    ---------+--------+
    Total     3568295  3568295
               100.00   100.00


     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_1Gy_0U_1     5081067       0.29735       0.41020       1510831             0       1.00000
     perc_methyl_1Gy_0U_2     5081067       0.30495       0.40880       1549447             0       1.00000


                                 Pearson Correlation Coefficients, N = 5081067
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          1Gy_0U_1          1Gy_0U_2

                            perc_methyl_1Gy_0U_1           1.00000           0.92965
                                                                              <.0001

                            perc_methyl_1Gy_0U_2           0.92965           1.00000
                                                            <.0001

                                                 The SAS System                                          19

       Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

       perc_methyl_1Gy_0U_1     3560334       0.34248       0.41945       1219338             0       1.00000
       perc_methyl_1Gy_0U_2     3560334       0.34610       0.41779       1232244             0       1.00000


                                   Pearson Correlation Coefficients, N = 3560334
                                             Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                            1Gy_0U_1          1Gy_0U_2

                              perc_methyl_1Gy_0U_1           1.00000           0.96308
                                                                                <.0001

                              perc_methyl_1Gy_0U_2           0.96308           1.00000
                                                              <.0001



*/

%concordance(CHG, 0Gy, 0U);
/*

    flag_rep1_0Gy_0U
              flag_rep2_0Gy_0U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           0 |   1230 |   1230
             |   0.03 |   0.03
             | 100.00 |
             |   0.03 |
    ---------+--------+
           1 |4273540 |4273540
             |  99.97 |  99.97
             | 100.00 |
             |  99.97 |
    ---------+--------+
    Total     4274770  4274770
               100.00   100.00
      Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

      perc_methyl_0Gy_0U_1     5650661       0.12442       0.24155        703079             0       1.00000
      perc_methyl_0Gy_0U_2     5650661       0.12518       0.24309        707376             0       1.00000


                                  Pearson Correlation Coefficients, N = 5650661
                                            Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                           0Gy_0U_1          0Gy_0U_2

                             perc_methyl_0Gy_0U_1           1.00000           0.84020
                                                                               <.0001

                             perc_methyl_0Gy_0U_2           0.84020           1.00000
                                                             <.0001

                                                  The SAS System                                          19:4


                                                   Simple Statistics

         Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

         perc_methyl_0Gy_0U_1     4273540       0.14149       0.24878        604670             0       1.00000
         perc_methyl_0Gy_0U_2     4273540       0.14179       0.25037        605957             0       1.00000


                                     Pearson Correlation Coefficients, N = 4273540
                                               Prob > |r| under H0: Rho=0

                                                          perc_methyl_      perc_methyl_
                                                              0Gy_0U_1          0Gy_0U_2

                                perc_methyl_0Gy_0U_1           1.00000           0.90438
                                                                                  <.0001

                                perc_methyl_0Gy_0U_2           0.90438           1.00000
                                                                <.0001

                                                     The SAS System                                          19:40

*/
%concordance(CHG, 01Gy, 0U);
/*
      flag_rep1_01Gy_0U
                flag_rep2_01Gy_0U

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       1|  Total
      ---------+--------+
             0 |    331 |    331
               |   0.01 |   0.01
               | 100.00 |
               |   0.01 |
      ---------+--------+
             1 |4946525 |4946525
               |  99.99 |  99.99
               | 100.00 |
               |  99.99 |
      ---------+--------+
      Total     4946856  4946856
                 100.00   100.00



                                              Simple Statistics

   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1     5837991       0.12251       0.22717        715213             0       1.00000
   perc_methyl_01Gy_0U_2     5837991       0.12044       0.22753        703117             0       1.00000


                               Pearson Correlation Coefficients, N = 5837991
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.84203
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.84203           1.00000
                                                           <.0001

                                               Simple Statistics

    Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_01Gy_0U_1     4946525       0.12825       0.22627        634400             0       1.00000
    perc_methyl_01Gy_0U_2     4946525       0.12620       0.22650        624268             0       1.00000


                                Pearson Correlation Coefficients, N = 4946525
                                          Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                         01Gy_0U_1         01Gy_0U_2

                           perc_methyl_01Gy_0U_1           1.00000           0.88876
                                                                              <.0001

                           perc_methyl_01Gy_0U_2           0.88876           1.00000
                                                            <.0001

                                                 The SAS System                                          19:40


*/
%concordance(CHG, 1Gy, 0U);
/*
  flag_rep1_1Gy_0U
            flag_rep2_1Gy_0U

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       1|  Total
  ---------+--------+
         0 |   7888 |   7888
           |   0.20 |   0.20
           | 100.00 |
           |   0.20 |
  ---------+--------+
         1 |3874859 |3874859
           |  99.80 |  99.80
           | 100.00 |
           |  99.80 |
  ---------+--------+
  Total     3882747  3882747
             100.00   100.00

        The SAS System


                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_1Gy_0U_1     5575382       0.12117       0.24744        675590             0       1.00000
  perc_methyl_1Gy_0U_2     5575382       0.13336       0.25278        743518             0       1.00000


                              Pearson Correlation Coefficients, N = 5575382
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       1Gy_0U_1          1Gy_0U_2

                         perc_methyl_1Gy_0U_1           1.00000           0.79521
                                                                           <.0001

                         perc_methyl_1Gy_0U_2           0.79521           1.00000
                                                         <.0001

                                                Simple Statistics

   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_1Gy_0U_1     3874859       0.14311       0.25683        554515             0       1.00000
   perc_methyl_1Gy_0U_2     3874859       0.15000       0.25710        581233             0       1.00000


                               Pearson Correlation Coefficients, N = 3874859
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        1Gy_0U_1          1Gy_0U_2

                          perc_methyl_1Gy_0U_1           1.00000           0.88591
                                                                            <.0001

                          perc_methyl_1Gy_0U_2           0.88591           1.00000
                                                          <.0001                                           The SAS System                                          19:40 S

*/

%concordance(CHH, 0Gy, 0U);
/*

  flag_rep1_0Gy_0U
            flag_rep2_0Gy_0U

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       1|  Total
  ---------+--------+
         0 |  10311 |  10311
           |   0.06 |   0.06
           | 100.00 |
           |   0.06 |
  ---------+--------+
         1 |1.829E7 |1.829E7
           |  99.94 |  99.94
           | 100.00 |
           |  99.94 |
  ---------+--------+
  Total     1.831E7  1.831E7
             100.00   100.00

                                              Simple Statistics

    Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_0Gy_0U_1    27234750       0.06308       0.14888       1717888             0       1.00000
    perc_methyl_0Gy_0U_2    27234750       0.06332       0.14614       1724404             0       1.00000


                               Pearson Correlation Coefficients, N = 27234750
                                          Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                         0Gy_0U_1          0Gy_0U_2

                           perc_methyl_0Gy_0U_1           1.00000           0.53149
                                                                             <.0001

                           perc_methyl_0Gy_0U_2           0.53149           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 S



                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_0Gy_0U_1    18294791       0.07021       0.14067       1284486             0       1.00000
  perc_methyl_0Gy_0U_2    18294791       0.06908       0.13745       1263805             0       1.00000


                             Pearson Correlation Coefficients, N = 18294791
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       0Gy_0U_1          0Gy_0U_2

                         perc_methyl_0Gy_0U_1           1.00000           0.72832
                                                                           <.0001

                         perc_methyl_0Gy_0U_2           0.72832           1.00000
                                                         <.0001

*/
%concordance(CHH, 01Gy, 0U);
/*

       flag_rep1_01Gy_0U
                 flag_rep2_01Gy_0U

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       1|  Total
       ---------+--------+
              0 |   3249 |   3249
                |   0.01 |   0.01
                | 100.00 |
                |   0.01 |
       ---------+--------+
              1 | 2.37E7 | 2.37E7
                |  99.99 |  99.99
                | 100.00 |
                |  99.99 |
       ---------+--------+
       Total      2.37E7   2.37E7
                  100.00   100.00

             The SAS System


                                                   Simple Statistics

        Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

        perc_methyl_01Gy_0U_1    29288053       0.06168       0.12649       1806543             0       1.00000
        perc_methyl_01Gy_0U_2    29288053       0.05897       0.12353       1727175             0       1.00000


                                    Pearson Correlation Coefficients, N = 29288053
                                              Prob > |r| under H0: Rho=0

                                                          perc_methyl_      perc_methyl_
                                                             01Gy_0U_1         01Gy_0U_2

                               perc_methyl_01Gy_0U_1           1.00000           0.51221
                                                                                  <.0001

                               perc_methyl_01Gy_0U_2           0.51221           1.00000
                                                                <.0001

                                                     The SAS System                                          19:

                                                   The CORR Procedure

                                                 Simple Statistics

      Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

      perc_methyl_01Gy_0U_1    23696673       0.06363       0.11653       1507771             0       1.00000
      perc_methyl_01Gy_0U_2    23696673       0.06093       0.11436       1443834             0       1.00000


                                  Pearson Correlation Coefficients, N = 23696673
                                            Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                           01Gy_0U_1         01Gy_0U_2

                             perc_methyl_01Gy_0U_1           1.00000           0.64156
                                                                                <.0001

                             perc_methyl_01Gy_0U_2           0.64156           1.00000
                                                              <.0001

                                                   The SAS System                                          19:40 S

*/
%concordance(CHH, 1Gy, 0U);
/*

      flag_rep1_1Gy_0U
                flag_rep2_1Gy_0U

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       1|  Total
      ---------+--------+
             0 |  76558 |  76558
               |   0.44 |   0.44
               | 100.00 |
               |   0.44 |
      ---------+--------+
             1 | 1.73E7 | 1.73E7
               |  99.56 |  99.56
               | 100.00 |
               |  99.56 |
      ---------+--------+
      Total     1.738E7  1.738E7
                 100.00   100.00

            The SAS System



                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_1Gy_0U_1    27286581       0.05769       0.14796       1574213             0       1.00000
  perc_methyl_1Gy_0U_2    27286581       0.09236       0.19675       2520142             0       1.00000


                             Pearson Correlation Coefficients, N = 27286581
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       1Gy_0U_1          1Gy_0U_2

                         perc_methyl_1Gy_0U_1           1.00000           0.41093
                                                                           <.0001

                         perc_methyl_1Gy_0U_2           0.41093           1.00000
                                                         <.0001

                                              Simple Statistics

    Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_0U_1    17299403       0.06694       0.14378       1158074             0       1.00000
    perc_methyl_1Gy_0U_2    17299403       0.08553       0.16447       1479686             0       1.00000


                               Pearson Correlation Coefficients, N = 17299403
                                          Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                         1Gy_0U_1          1Gy_0U_2

                           perc_methyl_1Gy_0U_1           1.00000           0.64831
                                                                             <.0001

                           perc_methyl_1Gy_0U_2           0.64831           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 S

*/

%concordance(GC, 0Gy, 0U);
/*

 flag_rep1_0Gy_0U
           flag_rep2_0Gy_0U

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        1 |4996963 |4996963
          | 100.00 | 100.00
          | 100.00 |
          | 100.00 |
 ---------+--------+
 Total     4996963  4996963
            100.00   100.00

       The SAS System



          Pearson Correlation Coefficients, N = 4996963
                    Prob > |r| under H0: Rho=0

                               perc_methyl_      perc_methyl_
                                   0Gy_0U_1          0Gy_0U_2

     perc_methyl_0Gy_0U_1           1.00000           0.90848
                                                       <.0001

     perc_methyl_0Gy_0U_2           0.90848           1.00000
                                     <.0001

                          The SAS System                                          19:4
*/
%concordance(GC, 01Gy, 0U);
/*

      flag_rep1_01Gy_0U
                flag_rep2_01Gy_0U

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       1|  Total
      ---------+--------+
             1 |5803543 |5803543
               | 100.00 | 100.00
               | 100.00 |
               | 100.00 |
      ---------+--------+
      Total     5803543  5803543
                 100.00   100.00

                                                Simple Statistics

     Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_01Gy_0U_1     5803543       0.09914       0.20332        575387             0       0.96874
     perc_methyl_01Gy_0U_2     5803543       0.10323       0.21691        599089             0       1.03335


                                 Pearson Correlation Coefficients, N = 5803543
                                           Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                          01Gy_0U_1         01Gy_0U_2

                            perc_methyl_01Gy_0U_1           1.00000           0.89538
                                                                               <.0001

                            perc_methyl_01Gy_0U_2           0.89538           1.00000
                                                             <.0001

                                                  The SAS System                                          19:40

*/
%concordance(GC, 1Gy, 0U);
/*

    flag_rep1_1Gy_0U
              flag_rep2_1Gy_0U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           1 |4539526 |4539526
             | 100.00 | 100.00
             | 100.00 |
             | 100.00 |
    ---------+--------+
    Total     4539526  4539526
               100.00   100.00


   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_1Gy_0U_1     4539526       0.11949       0.25791        542407             0       1.10532
   perc_methyl_1Gy_0U_2     4539526       0.10815       0.21407        490956             0       0.91300


                               Pearson Correlation Coefficients, N = 4539526
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        1Gy_0U_1          1Gy_0U_2

                          perc_methyl_1Gy_0U_1           1.00000           0.88470
                                                                            <.0001

                          perc_methyl_1Gy_0U_2           0.88470           1.00000
                                                          <.0001



                                           Simple Statistics

 Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

 perc_methyl_1Gy_0U_1     4539526       0.11949       0.25791        542407             0       1.10532
 perc_methyl_1Gy_0U_2     4539526       0.10815       0.21407        490956             0       0.91300


                             Pearson Correlation Coefficients, N = 4539526
                                       Prob > |r| under H0: Rho=0

                                                  perc_methyl_      perc_methyl_
                                                      1Gy_0U_1          1Gy_0U_2

                        perc_methyl_1Gy_0U_1           1.00000           0.88470
                                                                          <.0001

                        perc_methyl_1Gy_0U_2           0.88470           1.00000
                                                        <.0001

                                             The SAS System                                          19:40 Sunda
*/

%concordance(GC, 0Gy, 100U);
/*

         flag_rep1_0Gy_100U
                   flag_rep2_0Gy_100U

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       1|  Total
         ---------+--------+
                1 |4712204 |4712204
                  | 100.00 | 100.00
                  | 100.00 |
                  | 100.00 |
         ---------+--------+
         Total     4712204  4712204
                    100.00   100.00



                                                Simple Statistics

     Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_100U_1     4712204       0.29215       0.21867       1376651             0       0.99418
     perc_methyl_0Gy_100U_2     4712204       0.29735       0.23239       1401161             0       1.00588


                                  Pearson Correlation Coefficients, N = 4712204
                                            Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                          0Gy_100U_1        0Gy_100U_2

                            perc_methyl_0Gy_100U_1           1.00000           0.60706
                                                                                <.0001

                            perc_methyl_0Gy_100U_2           0.60706           1.00000
                                                              <.0001

*/
%concordance(GC, 01Gy, 100U);
/*

    flag_rep1_01Gy_100U
              flag_rep2_01Gy_100U

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       1|  Total
    ---------+--------+
           1 |5988381 |5988381
             | 100.00 | 100.00
             | 100.00 |
             | 100.00 |
    ---------+--------+
    Total     5988381  5988381
               100.00   100.00

          The SAS System

                                               Simple Statistics

   Variable                          N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_100U_1     5988381       0.44072       0.20608       2639227             0       1.06212
   perc_methyl_01Gy_100U_2     5988381       0.44296       0.18222       2652589             0       0.97123


                                Pearson Correlation Coefficients, N = 5988381
                                          Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                        01Gy_100U_1       01Gy_100U_2

                          perc_methyl_01Gy_100U_1           1.00000           0.46823
                                                                               <.0001

                          perc_methyl_01Gy_100U_2           0.46823           1.00000
                                                             <.0001

*/
%concordance(GC, 1Gy, 100U);
/*
          flag_rep1_1Gy_100U
                    flag_rep2_1Gy_100U

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       1|  Total
          ---------+--------+
                 1 |4703751 |4703751
                   | 100.00 | 100.00
                   | 100.00 |
                   | 100.00 |
          ---------+--------+
          Total     4703751  4703751
                     100.00   100.00

                                               Simple Statistics

    Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_100U_1     4703751       0.34516       0.21316       1623528             0       0.97086
    perc_methyl_1Gy_100U_2     4703751       0.34052       0.22484       1601728             0       1.03094


                                 Pearson Correlation Coefficients, N = 4703751
                                           Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                         1Gy_100U_1        1Gy_100U_2

                           perc_methyl_1Gy_100U_1           1.00000           0.55288
                                                                               <.0001

                           perc_methyl_1Gy_100U_2           0.55288           1.00000
                                                             <.0001



*/




/* 100X sites only and bin into 100bp windows */



%macro concordance(siteType,treatment,units);

%if &siteType.=GC %then %do;

data meth_data_rep1;
   set arabMAP.methylation_data_gc;
   where rep=1 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads; 
   rename  perc_methyl2=perc_methyl_&treatment._&units._1;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_gc;
   where rep=2 and site_type="&siteType." and treatment="&treatment." and units="&units." and  flag_normalized=1;;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl_norm * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads;
   rename perc_methyl2=perc_methyl_&treatment._&units._2;
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
   keep chr chr_bin perc_methyl2 flag_lt100_reads;
   rename perc_methyl2=perc_methyl_&treatment._&units._1 ;
run;


data meth_data_rep2;
   set arabMAP.methylation_data_cg_chg_chh;
   where rep=2 and site_type="&siteType." and treatment="&treatment."  ;
   chr_bin=int(stop_pos/100) + 1;
   perc_methyl2 = perc_methyl * 100;
   if total_C < 100 then flag_lt100_reads=1;
   else flag_lt100_reads=0;
   keep chr chr_bin perc_methyl2 flag_lt100_reads;
   rename perc_methyl2=perc_methyl_&treatment._&units._2 ;
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
  var perc_methyl_&treatment._&units._1;
  output out=mean_methyl_by_bin_rep1 (drop=_TYPE_ _FREQ_) mean=;
run;

proc means data=meth_data_rep2 noprint;
 by chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl_&treatment._&units._2;
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
  else do; flag_rep2_&treatment._&units.=0; perc_methyl_&treatment._&units._2=0; total_c_1=0; end;
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

data export_Data;
   set meth_data_compare;
   keep chr pos chr_bin perc_methyl_: ;
run;

proc export data=export_data
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/rep_concordance_&siteType._&treatment._&units._100X_100bp_binned.csv"
     dbms=csv replace;
run;

%mend;


%concordance(CG, 0Gy, 0U);
/*


                                               Simple Statistics

     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_0U_1       13009      50.14782      39.86885        652373             0     100.00000
     perc_methyl_0Gy_0U_2       13009      50.15524      39.67361        652469             0     100.00000


                                  Pearson Correlation Coefficients, N = 13009
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          0Gy_0U_1          0Gy_0U_2

                            perc_methyl_0Gy_0U_1           1.00000           0.99253
                                                                              <.0001

                            perc_methyl_0Gy_0U_2           0.99253           1.00000
                                                            <.0001


*/


%concordance(CG, 01Gy, 0U);
/*

                                               Simple Statistics

    Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_01Gy_0U_1        2636      34.14883      31.56339         90016             0      95.40990
    perc_methyl_01Gy_0U_2        2636      33.51767      31.61295         88353             0      95.80665


                                  Pearson Correlation Coefficients, N = 2636
                                          Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                         01Gy_0U_1         01Gy_0U_2

                           perc_methyl_01Gy_0U_1           1.00000           0.99318
                                                                              <.0001

                           perc_methyl_01Gy_0U_2           0.99318           1.00000
                                                            <.0001




*/


%concordance(CG, 1Gy, 0U);
/*


   Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_1Gy_0U_1        4195      40.39137      37.47797        169442             0     100.00000
   perc_methyl_1Gy_0U_2        4195      41.21943      37.27794        172915             0      97.34594


                                Pearson Correlation Coefficients, N = 4195
                                         Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                        1Gy_0U_1          1Gy_0U_2

                          perc_methyl_1Gy_0U_1           1.00000           0.99080
                                                                            <.0001

                          perc_methyl_1Gy_0U_2           0.99080           1.00000
                                                          <.0001




*/


%concordance(CHG, 0Gy, 0U);
/*



                                               Simple Statistics

     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_0U_1       12966      27.31133      28.11818        354119             0     100.00000
     perc_methyl_0Gy_0U_2       12966      27.37780      27.99949        354981             0      96.26168


                                  Pearson Correlation Coefficients, N = 12966
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          0Gy_0U_1          0Gy_0U_2

                            perc_methyl_0Gy_0U_1           1.00000           0.98143
                                                                              <.0001

                            perc_methyl_0Gy_0U_2           0.98143           1.00000
                                                            <.0001

                                                 The SAS System                                          19:40

*/


%concordance(CHG, 01Gy, 0U);
/*

                                             Simple Statistics

  Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_01Gy_0U_1        2427      16.73066      17.18876         40605             0      95.12805
  perc_methyl_01Gy_0U_2        2427      16.41853      17.29947         39848             0      98.71795


                                Pearson Correlation Coefficients, N = 2427
                                        Prob > |r| under H0: Rho=0

                                                    perc_methyl_      perc_methyl_
                                                       01Gy_0U_1         01Gy_0U_2

                         perc_methyl_01Gy_0U_1           1.00000           0.97850
                                                                            <.0001

                         perc_methyl_01Gy_0U_2           0.97850           1.00000
                                                          <.0001

                                               The SAS System                                          19:40 Su



*/


%concordance(CHG, 1Gy, 0U);
/*


                                               Simple Statistics

     Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_1Gy_0U_1        3872      18.80897      22.31672         72828             0      93.16239
     perc_methyl_1Gy_0U_2        3872      19.33267      21.92126         74856             0      93.03798


                                  Pearson Correlation Coefficients, N = 3872
                                           Prob > |r| under H0: Rho=0

                                                      perc_methyl_      perc_methyl_
                                                          1Gy_0U_1          1Gy_0U_2

                            perc_methyl_1Gy_0U_1           1.00000           0.96834
                                                                              <.0001

                            perc_methyl_1Gy_0U_2           0.96834           1.00000
                                                            <.0001

                                                 The SAS System                                          19:40 Su


*/


%concordance(CHH, 0Gy, 0U);
/*


                                          Simple Statistics

Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

perc_methyl_0Gy_0U_1       16905      11.81273      12.42799        199694             0      86.13861
perc_methyl_0Gy_0U_2       16905      11.73605      12.32780        198398             0      88.11881


                             Pearson Correlation Coefficients, N = 16905
                                      Prob > |r| under H0: Rho=0

                                                 perc_methyl_      perc_methyl_
                                                     0Gy_0U_1          0Gy_0U_2

                       perc_methyl_0Gy_0U_1           1.00000           0.91475
                                                                         <.0001

                       perc_methyl_0Gy_0U_2           0.91475           1.00000
                                                       <.0001



*/


%concordance(CHH, 01Gy, 0U);
/*

                                              Simple Statistics

   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1        3125       8.53391       8.06098         26668       0.88496      96.10853
   perc_methyl_01Gy_0U_2        3125       8.38016       8.32528         26188             0      96.25674


                                 Pearson Correlation Coefficients, N = 3125
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.96061
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.96061           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 S



*/


%concordance(CHH, 1Gy, 0U);
/*


                                              Simple Statistics

    Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_0U_1        5329       9.39021      10.68359         50040             0      84.61539
    perc_methyl_1Gy_0U_2        5329       9.83599      10.07930         52416             0      79.41176


                                 Pearson Correlation Coefficients, N = 5329
                                          Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                         1Gy_0U_1          1Gy_0U_2

                           perc_methyl_1Gy_0U_1           1.00000           0.85068
                                                                             <.0001

                           perc_methyl_1Gy_0U_2           0.85068           1.00000
                                                           <.0001

                                                The SAS System                                          19:40 Sun


*/


%concordance(GC, 0Gy, 0U);
/*


                                            Simple Statistics

  Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

  perc_methyl_0Gy_0U_1       13935      20.65310      23.55993        287801             0      93.50084
  perc_methyl_0Gy_0U_2       13935      22.98440      27.00505        320288             0     104.75641


                               Pearson Correlation Coefficients, N = 13935
                                        Prob > |r| under H0: Rho=0

                                                   perc_methyl_      perc_methyl_
                                                       0Gy_0U_1          0Gy_0U_2

                         perc_methyl_0Gy_0U_1           1.00000           0.94378
                                                                           <.0001

                         perc_methyl_0Gy_0U_2           0.94378           1.00000
                                                         <.0001

                                              The SAS System                                          19:40


*/


%concordance(GC, 01Gy, 0U);
/*


                                              Simple Statistics

   Variable                        N          Mean       Std Dev           Sum       Minimum       Maximum

   perc_methyl_01Gy_0U_1        2797      12.56083      13.82068         35133             0      92.79118
   perc_methyl_01Gy_0U_2        2797      13.15119      15.19590         36784             0     100.11650


                                 Pearson Correlation Coefficients, N = 2797
                                         Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                        01Gy_0U_1         01Gy_0U_2

                          perc_methyl_01Gy_0U_1           1.00000           0.96533
                                                                             <.0001

                          perc_methyl_01Gy_0U_2           0.96533           1.00000
                                                           <.0001


*/


%concordance(GC, 1Gy, 0U);
/*


                                          Simple Statistics

Variable                       N          Mean       Std Dev           Sum       Minimum       Maximum

perc_methyl_1Gy_0U_1        4416      17.78988      24.49459         78560             0     109.49917
perc_methyl_1Gy_0U_2        4416      14.69169      18.10247         64879             0      89.10914


                             Pearson Correlation Coefficients, N = 4416
                                      Prob > |r| under H0: Rho=0

                                                 perc_methyl_      perc_methyl_
                                                     1Gy_0U_1          1Gy_0U_2

                       perc_methyl_1Gy_0U_1           1.00000           0.92891
                                                                         <.0001

                       perc_methyl_1Gy_0U_2           0.92891           1.00000
                                                       <.0001


*/


%concordance(GC, 0Gy, 100U);
/*

                                                Simple Statistics

     Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

     perc_methyl_0Gy_100U_1        7702      37.89239      21.87415        291847             0      97.44328
     perc_methyl_0Gy_100U_2        7702      38.89964      22.09051        299605             0      98.89780


                                   Pearson Correlation Coefficients, N = 7702
                                            Prob > |r| under H0: Rho=0

                                                        perc_methyl_      perc_methyl_
                                                          0Gy_100U_1        0Gy_100U_2

                            perc_methyl_0Gy_100U_1           1.00000           0.93661
                                                                                <.0001

                            perc_methyl_0Gy_100U_2           0.93661           1.00000
                                                              <.0001



*/


%concordance(GC, 01Gy, 100U);
/*


                                             Simple Statistics

 Variable                          N          Mean       Std Dev           Sum       Minimum       Maximum

 perc_methyl_01Gy_100U_1        3406      63.06014      20.89049        214783             0     101.01176
 perc_methyl_01Gy_100U_2        3406      62.37673      19.94840        212455       0.39725      92.14254


                                Pearson Correlation Coefficients, N = 3406
                                        Prob > |r| under H0: Rho=0

                                                     perc_methyl_      perc_methyl_
                                                      01Gy_100U_1       01Gy_100U_2

                        perc_methyl_01Gy_100U_1           1.00000           0.97556
                                                                             <.0001

                        perc_methyl_01Gy_100U_2           0.97556           1.00000
                                                           <.0001

*/


%concordance(GC, 1Gy, 100U);

/*

                                               Simple Statistics

    Variable                         N          Mean       Std Dev           Sum       Minimum       Maximum

    perc_methyl_1Gy_100U_1        8295      40.27541      18.96855        334085             0      95.28848
    perc_methyl_1Gy_100U_2        8295      41.41053      20.28028        343500             0     102.18158


                                  Pearson Correlation Coefficients, N = 8295
                                           Prob > |r| under H0: Rho=0

                                                       perc_methyl_      perc_methyl_
                                                         1Gy_100U_1        1Gy_100U_2

                           perc_methyl_1Gy_100U_1           1.00000           0.93477
                                                                               <.0001

                           perc_methyl_1Gy_100U_2           0.93477           1.00000
                                                             <.0001





*/





