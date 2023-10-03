/* Cold WGBS QC:

(1) Binned rep-to-rep concordance:
        bin every 100bp, calc mean methylation for rep, rep concordance on windows
        (a) reps for same condition
        (b) 22C 100U vs FANS

(2) Estimate bisulfite conversion ratio:
        average % methylation on CHH sites on Chrom Pt by replicate */

libname cold '!PATCON/DTRA/arabidopsis_wgbs_cold/sas_data';
libname coldloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";
ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


/* Bin sites into 100bp windows */

data meth_data_all;
  set coldloc.methylation_data_by_rep;
  chr_bin=int(pos_end/100) + 1;
  if total_C < 100 then flag_lt100_reads=1;
  else flag_lt100_reads=0;
run;

proc sort data=meth_data_all;
  by sample_id temperature units rep site_type chr chr_bin;
proc means data=meth_data_all noprint;
  by sample_id temperature units rep site_type chr chr_bin;
  var perc_methyl;
  output out=mean_methyl_by_bin_all mean=;
run;

proc means data=meth_data_all noprint;
  by sample_id temperature units rep site_type chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl;
  output out=mean_methyl_by_bin_gt100 mean=;
run;

/* Calc correlations */
ods graphics / ANTIALIASMAX=50000000;

%macro corrPlot(siteType,sample1,sample2);

data sample1_all;
  set mean_methyl_by_bin_all;
  where sample_id = "&sample1." and site_type="&siteType";
  keep chr chr_bin perc_methyl;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample1_10;
  set mean_methyl_by_bin_gt100;
  where sample_id = "&sample1." and site_type="&siteType.";
  keep chr chr_bin perc_methyl;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample2_all;
  set mean_methyl_by_bin_all;
  where sample_id = "&sample2." and site_type="&siteType.";
  keep chr chr_bin perc_methyl;
  rename perc_methyl=perc_meth_&sample2.;
run;

data sample2_10;
  set mean_methyl_by_bin_gt100;
  where sample_id = "&sample2." and site_type="&siteType.";
  keep chr chr_bin perc_methyl;
  rename perc_methyl=perc_meth_&sample2.;
run;

proc sort data=sample1_all;
  by chr chr_bin;
  proc sort data=sample1_10;
  by chr chr_bin;
proc sort data=sample2_all;
  by chr chr_bin;
proc sort data=sample2_10;
  by chr chr_bin;
run;

data sample_1v2_all;
  merge sample1_all (in=in1) sample2_all (in=in2);
  by chr chr_bin;
run;

data sample_1v2_10;
  merge sample1_10 (in=in1) sample2_10 (in=in2);
  by chr chr_bin;
run;

title "Correlation &siteType. &sample1. &sample2. (all sites) ";
ods text="Correlation &siteType. &sample1. &sample2. (all sites) ";

proc corr data=sample_1v2_all pearson spearman polychoric polyserial kendall;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_all;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;


title "Correlation &siteType. &sample1. &sample2. (100X reads) ";
ods text="Correlation &siteType. &sample1. &sample2. (100X reads) ";

proc corr data=sample_1v2_10 pearson spearman polychoric polyserial kendall;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_10;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;


%mend;

%corrPlot(GC,FANS_0p5U_1,FANS_0p5U_2);
%corrPlot(GC,FANS_0p5U_1,FANS_1p5U_1);
%corrPlot(GC,FANS_0p5U_1,FANS_1p5U_2);
%corrPlot(GC,FANS_0p5U_1,FANS_5U_1);
%corrPlot(GC,FANS_0p5U_1,FANS_5U_2);
%corrPlot(GC,FANS_0p5U_1,FANS_25U_1);
%corrPlot(GC,FANS_0p5U_1,FANS_25U_2);
%corrPlot(GC,FANS_0p5U_2,FANS_1p5U_1);
%corrPlot(GC,FANS_0p5U_2,FANS_1p5U_2);
%corrPlot(GC,FANS_0p5U_2,FANS_5U_1);
%corrPlot(GC,FANS_0p5U_2,FANS_5U_2);
%corrPlot(GC,FANS_0p5U_2,FANS_25U_1);
%corrPlot(GC,FANS_0p5U_2,FANS_25U_2);
%corrPlot(GC,FANS_1p5U_1,FANS_1p5U_2);
%corrPlot(GC,FANS_1p5U_1,FANS_5U_1);
%corrPlot(GC,FANS_1p5U_1,FANS_5U_2);
%corrPlot(GC,FANS_1p5U_1,FANS_25U_1);
%corrPlot(GC,FANS_1p5U_1,FANS_25U_2);
%corrPlot(GC,FANS_1p5U_2,FANS_5U_1);
%corrPlot(GC,FANS_1p5U_2,FANS_5U_2);
%corrPlot(GC,FANS_1p5U_2,FANS_25U_1);
%corrPlot(GC,FANS_1p5U_2,FANS_25U_2);
%corrPlot(GC,FANS_5U_1,FANS_5U_2);
%corrPlot(GC,FANS_5U_1,FANS_25U_1);
%corrPlot(GC,FANS_5U_1,FANS_25U_2);
%corrPlot(GC,FANS_5U_2,FANS_25U_1);
%corrPlot(GC,FANS_5U_2,FANS_25U_2);
%corrPlot(GC,FANS_25U_1,FANS_25U_2);

