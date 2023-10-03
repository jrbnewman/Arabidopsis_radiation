/* For 22C and 4C, show and plot the average level of GC accessibility
   (calculated as 100U - 0U ) +/- 1kb around TSS for all genes */

libname cold '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname coldloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

* GC methylation data;
data gc_data;
  set coldloc.methylation_data;
  where site_type="GC" and temperature ne "FANS";
run;

proc sort data=gc_data;
  by chr pos_end condition;
proc transpose data=gc_data out=gc_sbys;
  by chr pos_end;
  id condition;
  var perc_methyl;
run;


proc transpose data=gc_data out=gc_sbys10;
  where total_C >= 10;
  by chr pos_end;
  id condition;
  var perc_methyl;
run;

data gc_sbys_2;
  set gc_sbys;
  if _22C_0U ne . and _22C_100U ne . then diff_22C_100U_0U=_22C_100U - _22C_0U;
  if _4C_0U ne . and _4C_100U ne . then diff_4C_100U_0U=_4C_100U - _4C_0U;
  if _22C_0U ne . and _22C_100U ne . and _4C_0U ne . and _4C_100U ne . then do;
        diff_22C_100U_0U_all=_22C_100U - _22C_0U;
        diff_4C_100U_0U_all=_4C_100U - _4C_0U;
  end;
  rename pos_end=pos;
run;


data gc_sbys10_2;
  set gc_sbys10;
  if _22C_0U ne . and _22C_100U ne . then diff_22C_100U_0U=_22C_100U - _22C_0U;
  if _4C_0U ne . and _4C_100U ne . then diff_4C_100U_0U=_4C_100U - _4C_0U;
  if _22C_0U ne . and _22C_100U ne . and _4C_0U ne . and _4C_100U ne . then do;
        diff_22C_100U_0U_all=_22C_100U - _22C_0U;
        diff_4C_100U_0U_all=_4C_100U - _4C_0U;
  end;
  rename pos_end=pos;
run;

/* Get TSSs and regions */

data exon2gene;
  set tair.tair20_exons_w_info;
  keep chrom start stop strand exon_id gene_id;
run;

*ID TSS;

proc sort data=exon2gene;
  by gene_id chrom start stop ;
run;

data tss_gene_plus;
  set exon2gene;
  by gene_id;
  if first.gene_id and strand="+" then output;
  keep gene_id chrom start strand;
  rename chrom=chr start=tss;
run;

proc sort data=exon2gene;
  by gene_id chrom descending stop descending start;
run;

data tss_gene_minus;
  set exon2gene;
  by gene_id;
  if first.gene_id and strand="-" then output;
  keep gene_id chrom stop strand;
  rename chrom=chr stop=tss;
run;

data tss_gene;
  set tss_gene_plus tss_gene_minus;
    tss_1kb_start=tss-1000;
    tss_1kb_stop=tss+1000;
run;


proc sort data=tss_gene;
  by gene_id;
run;

data tss_pos_index;
  set tss_gene;
  by gene_id ;
  do pos = tss_1kb_start to tss_1kb_stop ;
  pos_abs=pos-tss;
  output;
  end;
 run;

data tss_pos_index2;
  set tss_pos_index;
  if strand="-" then do;
       pos_abs2=pos_abs * -1;
    output;
    end;
  else do;
    pos_abs2=pos_abs;
    output;
    end;
  drop pos_abs;
  rename pos_abs2=pos_abs;
run;


proc sort data=tss_pos_index2;
  by chr pos;
proc sort data=gc_sbys_2;
  by chr pos;
proc sort data=gc_sbys10_2;
  by chr pos;
run;

data tss_pos_gc_all;
  merge gc_sbys_2 (in=in1) tss_pos_index2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data tss_pos_gc_10;
  merge gc_sbys10_2 (in=in1) tss_pos_index2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

/* Also going to group every 10 positions */
data tss_pos_gc_all_2;
  set tss_pos_gc_all;
  grouped_pos=int(pos_abs/10) * 10;
run;

data tss_pos_gc_10_2;
  set tss_pos_gc_10;
  grouped_pos=int(pos_abs/10) * 10;
run;

/* For each, calculate the mean and SD */

proc sort data=tss_pos_gc_all_2;
  by pos_abs grouped_pos;
proc sort data=tss_pos_gc_10_2;
  by pos_abs grouped_pos;
run;

proc means data=tss_pos_gc_all_2 noprint;
  by pos_abs grouped_pos;
  var diff_22C_100U_0U diff_4C_100U_0U diff_22C_100U_0U_all diff_4C_100U_0U_all;
  output out=mean_diff_tss_all_1000
  mean(diff_22C_100U_0U)=mean_22C_all   stddev(diff_22C_100U_0U)=sd_22C_all
  mean(diff_4C_100U_0U)=mean_4C_all   stddev(diff_4C_100U_0U)=sd_4C_all
  mean(diff_22C_100U_0U_all)=mean_22C_common   stddev(diff_22C_100U_0U_all)=sd_22C_common
  mean(diff_4C_100U_0U_all)=mean_4C_common   stddev(diff_4C_100U_0U_all)=sd_4C_common;
run;

proc means data=tss_pos_gc_10_2 noprint;
  by pos_abs grouped_pos;
  var diff_22C_100U_0U diff_4C_100U_0U diff_22C_100U_0U_all diff_4C_100U_0U_all;
  output out=mean_diff_tss_10X_1000
  mean(diff_22C_100U_0U)=mean_22C_all   stddev(diff_22C_100U_0U)=sd_22C_all
  mean(diff_4C_100U_0U)=mean_4C_all   stddev(diff_4C_100U_0U)=sd_4C_all
  mean(diff_22C_100U_0U_all)=mean_22C_common   stddev(diff_22C_100U_0U_all)=sd_22C_common
  mean(diff_4C_100U_0U_all)=mean_4C_common   stddev(diff_4C_100U_0U_all)=sd_4C_common;
run;



proc sort data=mean_diff_tss_all_1000;
  by grouped_pos;
proc sort data=mean_diff_tss_10X_1000;
  by grouped_pos;
run;

proc means data=mean_diff_tss_all_1000 noprint;
  by grouped_pos;
  var mean_22C_all mean_4C_all mean_22C_common mean_4C_common;
  output out=mean_diff_tss_all_100
  mean(mean_22C_all)=mean_22C_all   stddev(mean_22C_all)=sd_22C_all
  mean(mean_4C_all)=mean_4C_all   stddev(mean_4C_all)=sd_4C_all
  mean(mean_22C_common)=mean_22C_common   stddev(mean_22C_common)=sd_22C_common
  mean(mean_4C_common)=mean_4C_common   stddev(mean_4C_common)=sd_4C_common;
run;

proc means data=mean_diff_tss_10X_1000 noprint;
  by grouped_pos;
  var mean_22C_all mean_4C_all mean_22C_common mean_4C_common;
  output out=mean_diff_tss_10X_100
  mean(mean_22C_all)=mean_22C_all   stddev(mean_22C_all)=sd_22C_all
  mean(mean_4C_all)=mean_4C_all   stddev(mean_4C_all)=sd_4C_all
  mean(mean_22C_common)=mean_22C_common   stddev(mean_22C_common)=sd_22C_common
  mean(mean_4C_common)=mean_4C_common   stddev(mean_4C_common)=sd_4C_common;
run;

data mean_diff_tss_all_1000_2;
  set mean_diff_tss_all_1000;
  drop grouped_pos;
  rename pos_abs=pos;
run;

data mean_diff_tss_10X_1000_2;
  set mean_diff_tss_10X_1000;
  drop grouped_pos;
  rename pos_abs=pos;
run;



data mean_diff_tss_all_100_2;
  set mean_diff_tss_all_100;
  rename grouped_pos=pos;
run;

data mean_diff_tss_10X_100_2;
  set mean_diff_tss_10X_100;
  rename grouped_pos=pos;
run;


/* Export and plot in python */

proc export data=mean_diff_tss_all_1000_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/GC_accessibility_1kb_TSS_all_sites_nobin.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_10X_1000_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/GC_accessibility_1kb_TSS_10X_sites_nobin.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_all_100_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/GC_accessibility_1kb_TSS_all_sites_binned.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_10X_100_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/GC_accessibility_1kb_TSS_10X_sites_binned.csv"
   dbms=csv
   replace;
run;

