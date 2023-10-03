libname wgbs "!PATCON/arabidopsis_wgbs_cold/sas_data";
libname wgbsloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";
ods html close;
ods listing;

/* Create DMAP2 summaries for DMCs/DACs:
    -> full list of tested sites
    -> DMCs/DACs only (FDR05)

    Pull in annotations where possible (i.e. promoter region, gene etc)

    Count DMCs by dose, region, site
        (1) P<0.05, FDR P<0.05
        (2) Up/down, diffs any, 10%, 20%, 30%, 50%

    Count overlap of:
        (1) Tested sites (CG, CHG, CHH, GC)
        (2) Significant DMCs (by dose, region, and site type)

        -> Up/down: DMCs, DMRs and DMGs (# with up, # with down, # both, average direction)
        -> Diff 10%, 20%, 30%
        -> For DMCs/DMRs/DMGs: significance threshold (+/- minimum difference)
        -> For DACs/DARs/DAGs: either 100Uv0U test to be signficant (+/- minimum diff of diffs)
        -> DMAP2 overlaps (P05, FDR05; w/wo logical methylation diff; ignore splitting on direction for now)
        by site:
        -> sites in common between site types
        -> DMCs in common between dosages


    Enrichment of orphan DMCs vs non-orphan DMCs
    Enrichment of promoter DMCs vs non-promoter DMCs
    Enrichment of gene body DMCs vs non-gene body DMCs
*/
proc datasets lib=work kill noprint;
run;
quit;


/* (1) Create summary of all DMCs : both dosage, all three native methylation sites */

data dmc_results;
  length site_type $3.;
  set wgbsloc.results_by_dmc_cg_4C_22C (in=in1)
      wgbsloc.results_by_dmc_chg_4C_22C (in=in2)
      wgbsloc.results_by_dmc_chh_4C_22C (in=in3);
  if in1 then site_type="CG";
  else if in2 then site_type="CHG";
  else if in3 then site_type="CHH";
      drop call pval;
      rename p_value=P_4C_22C
             fdr_p=FDR_P_4C_22C
             flag_p05=flag_p05_4C_22C
             flag_fdr05=flag_fdr05_4C_22C
             flag_diff10=flag_diff_4C_22C_10prc
             flag_diff20=flag_diff_4C_22C_20prc
             flag_diff30=flag_diff_4C_22C_30prc
             flag_diff50=flag_diff_4C_22C_50prc;
run;

/* Want to add coordinates too */

data region_coord;
   set wgbsloc.meth_pos2region;
   length gene_id $15.;
   length region_type $10.;
   length strand $10.;
   format region_start best12.;
   format region_stop best12.;
   do i=1 by 1 while(scan(gene_id_cat,i,"|") ^= "");
       gene_id=compress(scan(gene_id_cat,i,"|"));
       region_type=compress(scan(region_type_cat,i,"|"));
       strand=compress(scan(strand_cat,i,"|"));
       region_start=scan(region_start_cat,i,"|") * 1;
       region_stop=scan(region_stop_cat,i,"|") * 1;
       output;
       end;
    keep chr pos gene_id region_type strand region_start region_stop;
run;

data region_coord2;
   set region_coord;
   format region_id $50.;
   region_id=catx(":",gene_id,region_type,chr,region_start,region_stop,strand);
   keep chr pos gene_id region_type region_id;
run;

proc sort data=region_coord2 nodup;
  by chr pos gene_id region_type;
proc sort data=dmc_results;
  by chr pos gene_id region_type;
run;

data dmc_results2;
  merge dmc_results (in=in1) region_coord2;
  by chr pos gene_id region_type;
  if in1;
run;


* collapse geneIDs and regions;

proc sort data=dmc_results2;
  by site_type chr pos;
proc freq data=dmc_results2 noprint;
  by site_type chr;
  tables pos / out=obs_per_pos;
proc sort data=obs_per_pos;
  by descending count;
run; *14 max;


data _null_;
  set obs_per_pos (firstobs=1 obs=1);
  call symputx('regionCnt',count);
  stop;
run;

%put &regionCnt.;


data dmc_results_cat;
   array gene[&regionCnt.] $15.;
   array region[&regionCnt.] $9.;
   array regionID[&regionCnt.] $50.;
   retain gene1-gene&regionCnt. region1-region&regionCnt. regionID1-regionID&regionCnt. ;
   set dmc_results2;
   by site_type chr pos;
   if first.pos then do;
      call missing(of gene1-gene&regionCnt.);
      call missing(of region1-region&regionCnt.);
      records=0;
      end;
   records + 1;
   gene[records]=gene_id;
   region[records]=region_type;
   regionID[records]=region_id;
   if last.pos then output;
run;

data dmc_results_cat2;
   set dmc_results_cat;
   length gene_id_cat $224.;
   length region_type_cat $140.;
   length region_id_cat $700.;
   gene_id_cat = catx("|", OF gene1-gene&regionCnt.);
   region_type_cat = catx("|", OF region1-region&regionCnt.);
   region_id_cat = catx("|", OF  regionID1-regionID&regionCnt.);
   rename records=region_count;
   drop gene1-gene&regionCnt. region1-region&regionCnt.  regionID1-regionID&regionCnt. gene_id region_type  region_id;
run;

proc datasets lib=work memtype=data;
   modify dmc_results_cat2;
     attrib _all_ label=' ';
run;
quit;


proc contents data=dmc_results_cat2;
run;

data dmc_results_cat3;
   retain chr pos site_type
          meth_22C_0U meth_4C_0U
          P_4C_22C flag_p05_4C_22C
          FDR_P_4C_22C flag_fdr05_4C_22C
          diff_4C_22C flag_diff_4C_22C_10prc flag_diff_4C_22C_20prc
                      flag_diff_4C_22C_30prc flag_diff_4C_22C_50prc
          region_count gene_id_cat flag_orphan_site region_type_cat region_id_cat;
   format meth_22C_0U f8.3;
   format meth_4C_0U f8.3;
   format P_4C_22C best12.;
   format FDR_P_4C_22C best12.;
   format diff_4C_22C f8.3 ;
   set dmc_results_cat2;
   if flag_orphan_site=1 then region_count=0;
run;

/* (2) Subset significant DMCs only */

data sig_dmcs;
  set dmc_results_cat3;
  where FDR_P_4C_22C < 0.05 ;
run;


/* (3) Create summary of all DACs : both dosage, GC sites */


data dac_results;
  length site_type $3.;
  set wgbsloc.results_by_dmc_gc_4C_22C;
  site_type="GC";
  drop methdiff_4C_100U_0U methdiff_22C_100U_0U;
      rename flag_diff10=flag_diff_4C_22C_10prc
             flag_diff20=flag_diff_4C_22C_20prc
             flag_diff30=flag_diff_4C_22C_30prc
             flag_diff50=flag_diff_4C_22C_50prc;
run;

proc sort data=region_coord2 nodup;
  by chr pos gene_id region_type;
proc sort data=dac_results;
  by chr pos gene_id region_type;
run;

data dac_results2;
  merge dac_results (in=in1) region_coord2;
  by chr pos gene_id region_type;
  if in1;
run;


* collapse geneIDs and regions;

proc sort data=dac_results2;
  by site_type chr pos;
proc freq data=dac_results2 noprint;
  by site_type chr;
  tables pos / out=obs_per_pos;
proc sort data=obs_per_pos;
  by descending count;
run; *14 max;


data _null_;
  set obs_per_pos (firstobs=1 obs=1);
  call symputx('regionCnt',count);
  stop;
run;

%put &regionCnt.;


data dac_results_cat;
   array gene[&regionCnt.] $15.;
   array region[&regionCnt.] $9.;
   array regionID[&regionCnt.] $50.;
   retain gene1-gene&regionCnt. region1-region&regionCnt. regionID1-regionID&regionCnt. ;
   set dac_results2;
   by site_type chr pos;
   if first.pos then do;
      call missing(of gene1-gene&regionCnt.);
      call missing(of region1-region&regionCnt.);
      records=0;
      end;
   records + 1;
   gene[records]=gene_id;
   region[records]=region_type;
   regionID[records]=region_id;
   if last.pos then output;
run;

data dac_results_cat2;
   set dac_results_cat;
   length gene_id_cat $224.;
   length region_type_cat $140.;
   length region_id_cat $700.;
   gene_id_cat = catx("|", OF gene1-gene&regionCnt.);
   region_type_cat = catx("|", OF region1-region&regionCnt.);
   region_id_cat = catx("|", OF  regionID1-regionID&regionCnt.);
   rename records=region_count;
   drop gene1-gene&regionCnt. region1-region&regionCnt.  regionID1-regionID&regionCnt. gene_id region_type  region_id;
run;

proc datasets lib=work memtype=data;
   modify dac_results_cat2;
     attrib _all_ label=' ';
run;
quit;


data dac_results_cat3;
   retain chr pos site_type meth_22C_0U meth_22C_100U meth_4C_0U meth_4C_100U 
          P_22C_100U_0U flag_P05_22C_100U_0U fdr_22C_100U_0U flag_FDR05_22C_100U_0U
          P_4C_100U_0U flag_P05_4C_100U_0U fdr_4C_100U_0U flag_FDR05_4C_100U_0U
          diff_22C_100U_0U diff_4C_100U_0U flag_22C_0U_gt_100U  flag_4C_0U_gt_100U
          diff_4C_22C flag_22C_diff10 flag_22C_diff20 flag_22C_diff30 flag_22C_diff50
          flag_4C_diff10 flag_4C_diff20 flag_4C_diff30 flag_4C_diff50
          flag_diff_4C_22C_10prc flag_diff_4C_22C_20prc  flag_diff_4C_22C_30prc flag_diff_4C_22C_50prc
          region_count flag_orphan_site gene_id_cat region_type_cat region_id_cat;
   format meth_22C_0U f8.3;
   format meth_22C_100U f8.3;
   format meth_4C_0U f8.3;
   format meth_4C_100U f8.3;
   format P_22C_100U_0U best12.;
   format fdr_22C_100U_0U best12.;
   format P_4C_100U_0U best12.;
   format fdr_4C_100U_0U best12.;
   format diff_22C_100U_0U f8.3 diff_4C_100U_0U f8.3  diff_4C_22C f8.3 ;
   set dac_results_cat2;
   if flag_orphan_site=1 then region_count=0;
   drop flag_22C_100U_lt_0U flag_4C_100U_lt_0U;
   rename fdr_22C_100U_0U=FDR_P_22C_100U_0U
          fdr_4C_100U_0U=FDR_P_4C_100U_0U;
run;

/* (4) Subset significant DACs only */


data sig_dacs;
  set dac_results_cat3;
  if FDR_P_22C_100U_0U < 0.05 or  FDR_P_4C_100U_0U < 0.05 ;
run;

data dac_results_cat3_noGT;
  set dac_results_cat3;
  if flag_22C_0U_gt_100U=1 or flag_4C_0U_gt_100U=1 then delete;
run;


data sig_dacs_noGT;
  set dac_results_cat3_noGT;
  if FDR_P_22C_100U_0U < 0.05 or  FDR_P_4C_100U_0U < 0.05;
run;

/* COUNTS */


proc freq data=dmc_results_cat3 noprint;
  tables site_type*flag_fdr05_4C_22C / out=dmc_counts;
  tables site_type*flag_fdr05_4C_22C*flag_diff_4C_22C_10prc / out=dmc_counts_10perc;
  tables site_type*flag_fdr05_4C_22C*flag_diff_4C_22C_20prc / out=dmc_counts_20perc;
  tables site_type*flag_fdr05_4C_22C*flag_diff_4C_22C_30prc / out=dmc_counts_30perc;
  tables site_type*flag_fdr05_4C_22C*flag_diff_4C_22C_50prc / out=dmc_counts_50perc;
run;

proc print data=dmc_counts; run;
proc print data=dmc_counts_10perc; run;
proc print data=dmc_counts_20perc; run;
proc print data=dmc_counts_30perc; run;
proc print data=dmc_counts_50perc; run;


proc sort data=dmc_results_cat3;
  by site_type;
proc freq data=dmc_results_cat3 noprint;
  by site_type;
  tables flag_fdr05_4C_22C*flag_orphan_site /out=orphan;
  tables flag_fdr05_4C_22C*flag_diff_4C_22C_10prc*flag_orphan_site /out=orphan_10;
  tables flag_fdr05_4C_22C*flag_diff_4C_22C_20prc*flag_orphan_site /out=orphan_20;
  tables flag_fdr05_4C_22C*flag_diff_4C_22C_30prc*flag_orphan_site /out=orphan_30;
  tables flag_fdr05_4C_22C*flag_diff_4C_22C_50prc*flag_orphan_site /out=orphan_50;
run;

proc print data=orphan;
proc print data=orphan_10;
proc print data=orphan_20;
proc print data=orphan_30;
proc print data=orphan_50;
run;                

proc freq data=dac_results_cat3 noprint;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U / out=dac;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_10prc / out=dac_10;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff10*flag_4C_diff10*flag_diff_4C_22C_10prc / out=dac_10_all;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_20prc / out=dac_20;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff20*flag_4C_diff20*flag_diff_4C_22C_20prc / out=dac_20_all;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_30prc / out=dac_30;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff30*flag_4C_diff30*flag_diff_4C_22C_30prc / out=dac_30_all;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_50prc / out=dac_50;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff50*flag_4C_diff50*flag_diff_4C_22C_50prc / out=dac_50_all;



  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U*flag_orphan_site / out=dac_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_10prc*flag_orphan_site / out=dac_10_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff10*flag_4C_diff10*flag_diff_4C_22C_10prc*flag_orphan_site / out=dac_10_all_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_20prc*flag_orphan_site / out=dac_20_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff20*flag_4C_diff20*flag_diff_4C_22C_20prc*flag_orphan_site / out=dac_20_all_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_30prc*flag_orphan_site / out=dac_30_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff30*flag_4C_diff30*flag_diff_4C_22C_30prc*flag_orphan_site / out=dac_30_all_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_diff_4C_22C_50prc*flag_orphan_site / out=dac_50_orphan;
  tables flag_FDR05_22C_100U_0U*flag_FDR05_4C_100U_0U*
         flag_22C_0U_gt_100U*flag_4C_0U_gt_100U
         *flag_22C_diff50*flag_4C_diff50*flag_diff_4C_22C_50prc*flag_orphan_site / out=dac_50_all_orphan;
run;


proc print data=dac;
proc print data=dac_10;
proc print data=dac_10_all;
proc print data=dac_20;
proc print data=dac_20_all;
proc print data=dac_30;
proc print data=dac_30_all;
proc print data=dac_50;
proc print data=dac_50_all;
proc print data=dac_orphan;
proc print data=dac_10_orphan;
proc print data=dac_10_all_orphan;
proc print data=dac_20_orphan;
proc print data=dac_20_all_orphan;
proc print data=dac_30_orphan;
proc print data=dac_30_all_orphan;
proc print data=dac_50_orphan;
proc print data=dac_50_all_orphan;
run;

/* Make perm */

data wgbsloc.dmc_results_combined_annot;
   set dmc_results_cat3;
run;

data wgbsloc.dac_results_combined_annot;
   set dac_results_cat3;
run;

data wgbsloc.dac_results_combined_annot_noGT;
   set dac_results_cat3_noGT;
run;



/* EXPORT */

proc export data=dmc_results_cat3
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/DMAP2_native_methylation_sites_CG_CHG_CHH_min_10_reads.csv"
     dbms=csv replace;
run;

proc export data=sig_dmcs
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/DMAP2_native_methylation_sites_CG_CHG_CHH_min_10_reads_DMCs_FDR05.csv"
     dbms=csv replace;
run;

proc export data=dac_results_cat3
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/DMAP2_accessibility_sites_GC_min_10_reads.csv"
     dbms=csv replace;
run;

proc export data=sig_dacs
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/DMAP2_accessibility_sites_GC_min_10_reads_DACs_FDR05.csv"
     dbms=csv replace;
run;


proc export data=dac_results_cat3_noGT
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/DMAP2_accessibility_sites_GC_min_10_reads_no_0U_gt_100U.csv"
     dbms=csv replace;
run;

proc export data=sig_dacs_noGT
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/DMAP2_accessibility_sites_GC_min_10_reads_no_0U_gt_100U_DACs_FDR05.csv"
     dbms=csv replace;
run;




/* Count: Overlap between sites and DMCs/DACs by site_Type */

data cg chg chh;
set dmc_results_cat3;
if site_type="CG" then output cg;
if site_type="CHG" then output chg;
if site_type="CHH" then output chh;
keep chr pos;
run;


data gc;
  set dac_results_cat3;
  keep chr pos;
run;

proc sort data=cg nodup;
  by chr pos;
proc sort data=chg nodup;
  by chr pos;
proc sort data=chh nodup;
  by chr pos;
proc sort data=gc nodup;
  by chr pos;
run;

data site_compare;
  merge cg (in=in1) chg (in=in2) chh (in=in3) gc (in=in4);
  by chr pos;
  if in1 then flag_cg=1; else flag_cg=0;
  if in2 then flag_chg=1; else flag_chg=0;
  if in3 then flag_chh=1; else flag_chh=0;
  if in4 then flag_gc=1; else flag_gc=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg*flag_chg*flag_chh*flag_gc / out=site_count;
proc print data=site_count;
run;

data cg cg_10 cg_20 cg_30 cg_50
     chg chg_10 chg_20 chg_30 chg_50
     chh chh_10 chh_20 chh_30 chh_50;
     set dmc_results_cat3;
     if site_type="CG" and flag_fdr05_4C_22C=1 then output cg;
     if site_type="CG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_10prc=1 then output cg_10;
     if site_type="CG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_20prc=1 then output cg_20;
     if site_type="CG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_30prc=1 then output cg_30;
     if site_type="CG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_50prc=1 then output cg_50;

     if site_type="CHG" and flag_fdr05_4C_22C=1 then output chg;
     if site_type="CHG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_10prc=1 then output chg_10;
     if site_type="CHG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_20prc=1 then output chg_20;
     if site_type="CHG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_30prc=1 then output chg_30;
     if site_type="CHG" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_50prc=1 then output chg_50;

     if site_type="CHH" and flag_fdr05_4C_22C=1 then output chh;
     if site_type="CHH" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_10prc=1 then output chh_10;
     if site_type="CHH" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_20prc=1 then output chh_20;
     if site_type="CHH" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_30prc=1 then output chh_30;
     if site_type="CHH" and flag_fdr05_4C_22C=1 and flag_diff_4C_22C_50prc=1 then output chh_50;

keep chr pos;
run;



data gc gc_10 gc_10_all gc_20 gc_20_all gc_30 gc_30_all gc_50 gc_50_all;
  set dac_results_cat3;
  if flag_4C_0U_gt_100U=1 or flag_22C_0U_gt_100U=1 then delete;

  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) then output gc;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_10prc=1 then output gc_10;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_20prc=1 then output gc_20;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_30prc=1 then output gc_30;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_50prc=1 then output gc_50;

  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_10prc=1 
     and  flag_22C_diff10=1 and flag_4C_diff10=1 then output gc_10_all;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_20prc=1 
     and  flag_22C_diff20=1 and flag_4C_diff20=1 then output gc_20_all;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_30prc=1 
     and  flag_22C_diff30=1 and flag_4C_diff30=1 then output gc_30_all;
  if (flag_FDR05_22C_100U_0U=1 or flag_FDR05_4C_100U_0U=1) and flag_diff_4C_22C_50prc=1 
     and  flag_22C_diff50=1 and flag_4C_diff50=1 then output gc_50_all;

  keep chr pos;
run;

proc sort data=cg nodup; by chr pos;
proc sort data=cg_10 nodup; by chr pos;
proc sort data=cg_20 nodup; by chr pos;
proc sort data=cg_30 nodup; by chr pos;
proc sort data=cg_50 nodup; by chr pos;
proc sort data=chg nodup; by chr pos;
proc sort data=chg_10 nodup; by chr pos;
proc sort data=chg_20 nodup; by chr pos;
proc sort data=chg_30 nodup; by chr pos;
proc sort data=chg_50 nodup; by chr pos;
proc sort data=chh nodup; by chr pos;
proc sort data=chh_10 nodup; by chr pos;
proc sort data=chh_20 nodup; by chr pos;
proc sort data=chh_30 nodup; by chr pos;
proc sort data=chh_50 nodup; by chr pos;
proc sort data=gc nodup; by chr pos;
proc sort data=gc_10 nodup; by chr pos;
proc sort data=gc_10_all nodup; by chr pos;
proc sort data=gc_20 nodup; by chr pos;
proc sort data=gc_20_all nodup; by chr pos;
proc sort data=gc_30 nodup; by chr pos;
proc sort data=gc_30_all nodup; by chr pos;
proc sort data=gc_50 nodup; by chr pos;
proc sort data=gc_50_all nodup; by chr pos;
run;

data site_compare;
  merge cg (in=in1) chg (in=in2) chh (in=in3) gc (in=in4);
  by chr pos;
  if in1 then flag_cg=1; else flag_cg=0;
  if in2 then flag_chg=1; else flag_chg=0;
  if in3 then flag_chh=1; else flag_chh=0;
  if in4 then flag_gc=1; else flag_gc=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg*flag_chg*flag_chh*flag_gc / out=site_count;
proc print data=site_count;
run;




data site_compare;
  merge cg_10 (in=in1) chg_10 (in=in2) chh_10 (in=in3) gc_10 (in=in4);
  by chr pos;
  if in1 then flag_cg_10=1; else flag_cg_10=0;
  if in2 then flag_chg_10=1; else flag_chg_10=0;
  if in3 then flag_chh_10=1; else flag_chh_10=0;
  if in4 then flag_gc_10=1; else flag_gc_10=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_10*flag_chg_10*flag_chh_10*flag_gc_10 / out=site_count;
proc print data=site_count;
run;


data site_compare;
  merge cg_10 (in=in1) chg_10 (in=in2) chh_10 (in=in3) gc_10_all (in=in4);
  by chr pos;
  if in1 then flag_cg_10=1; else flag_cg_10=0;
  if in2 then flag_chg_10=1; else flag_chg_10=0;
  if in3 then flag_chh_10=1; else flag_chh_10=0;
  if in4 then flag_gc_10_all=1; else flag_gc_10_all=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_10*flag_chg_10*flag_chh_10*flag_gc_10_all / out=site_count;
proc print data=site_count;
run;


data site_compare;
  merge cg_20 (in=in1) chg_20 (in=in2) chh_20 (in=in3) gc_20 (in=in4);
  by chr pos;
  if in1 then flag_cg_20=1; else flag_cg_20=0;
  if in2 then flag_chg_20=1; else flag_chg_20=0;
  if in3 then flag_chh_20=1; else flag_chh_20=0;
  if in4 then flag_gc_20=1; else flag_gc_20=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_20*flag_chg_20*flag_chh_20*flag_gc_20 / out=site_count;
proc print data=site_count;
run;


data site_compare;
  merge cg_20 (in=in1) chg_20 (in=in2) chh_20 (in=in3) gc_20_all (in=in4);
  by chr pos;
  if in1 then flag_cg_20=1; else flag_cg_20=0;
  if in2 then flag_chg_20=1; else flag_chg_20=0;
  if in3 then flag_chh_20=1; else flag_chh_20=0;
  if in4 then flag_gc_20_all=1; else flag_gc_20_all=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_20*flag_chg_20*flag_chh_20*flag_gc_20_all / out=site_count;
proc print data=site_count;
run;


data site_compare;
  merge cg_30 (in=in1) chg_30 (in=in2) chh_30 (in=in3) gc_30 (in=in4);
  by chr pos;
  if in1 then flag_cg_30=1; else flag_cg_30=0;
  if in2 then flag_chg_30=1; else flag_chg_30=0;
  if in3 then flag_chh_30=1; else flag_chh_30=0;
  if in4 then flag_gc_30=1; else flag_gc_30=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_30*flag_chg_30*flag_chh_30*flag_gc_30 / out=site_count;
proc print data=site_count;
run;


data site_compare;
  merge cg_30 (in=in1) chg_30 (in=in2) chh_30 (in=in3) gc_30_all (in=in4);
  by chr pos;
  if in1 then flag_cg_30=1; else flag_cg_30=0;
  if in2 then flag_chg_30=1; else flag_chg_30=0;
  if in3 then flag_chh_30=1; else flag_chh_30=0;
  if in4 then flag_gc_30_all=1; else flag_gc_30_all=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_30*flag_chg_30*flag_chh_30*flag_gc_30_all / out=site_count;
proc print data=site_count;
run;



data site_compare;
  merge cg_50 (in=in1) chg_50 (in=in2) chh_50 (in=in3) gc_50 (in=in4);
  by chr pos;
  if in1 then flag_cg_50=1; else flag_cg_50=0;
  if in2 then flag_chg_50=1; else flag_chg_50=0;
  if in3 then flag_chh_50=1; else flag_chh_50=0;
  if in4 then flag_gc_50=1; else flag_gc_50=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_50*flag_chg_50*flag_chh_50*flag_gc_50 / out=site_count;
proc print data=site_count;
run;


data site_compare;
  merge cg_50 (in=in1) chg_50 (in=in2) chh_50 (in=in3) gc_50_all (in=in4);
  by chr pos;
  if in1 then flag_cg_50=1; else flag_cg_50=0;
  if in2 then flag_chg_50=1; else flag_chg_50=0;
  if in3 then flag_chh_50=1; else flag_chh_50=0;
  if in4 then flag_gc_50_all=1; else flag_gc_50_all=0;
run;

proc freq data=site_compare noprint;
  tables flag_cg_50*flag_chg_50*flag_chh_50*flag_gc_50_all / out=site_count;
proc print data=site_count;
run;

