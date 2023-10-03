ods listing; ods html close;
libname wgbs '!PATCON/arabidopsis_wgbs/sas_data';
libname wgbsloc "/TB14/TB14/sandbox/wgbs_sandbox/sas_data";

/* Flag sites with at least one methylated read in at least one rep 
   Only look at sites with sufficient coverage */

data sites2keep;
  set wgbs.flag_coverage_ge10;
  if site_type="GC" then do;
     if flag_0Gy_0U_coverage_ge10=1 and flag_0Gy_100U_coverage_ge10=1 
    and flag_01Gy_0U_coverage_ge10=1 and flag_01Gy_100U_coverage_ge10=1 
    and flag_1Gy_0U_coverage_ge10=1 and flag_1Gy_100U_coverage_ge10=1 
    then output;
    end;
  else if site_type="GCH" then do;
     if flag_0Gy_0U_coverage_ge10=1 and flag_0Gy_100U_coverage_ge10=1 
    and flag_01Gy_0U_coverage_ge10=1 and flag_01Gy_100U_coverage_ge10=1 
    and flag_1Gy_0U_coverage_ge10=1 and flag_1Gy_100U_coverage_ge10=1 
    then output;
    end;
  else do;
     if flag_0Gy_0U_coverage_ge10=1 and flag_01Gy_0U_coverage_ge10=1
     and flag_1Gy_0U_coverage_ge10=1
    then output;
    end;
 keep chr pos site_type;
run;


data meth_data;
  set wgbsloc.methylation_data;
  if chr="" then delete;
  keep sample_id dose units rep site_type chr pos_end methyl_C total_C;
  rename pos_end=pos;
run;

proc sort data=meth_data;
  by chr pos site_type;
proc sort data=sites2keep;
  by chr pos site_type;
run;

data meth_data2;
  merge sites2keep (in=in1) meth_data (in=in2);
  by chr pos site_type;
  if in1 and in2;
run;

data flag_methylation;
  set meth_data2;
  if methyl_C > 1 then flag_methyl=1;
  else flag_methyl=0;
  run;

proc sort data=flag_methylation;
   by dose units site_type chr pos;
proc means data=flag_methylation noprint;
   by dose units site_type chr pos;
   var flag_methyl;
   output out=perc_methyl mean=;
run;


data flag_perc_methyl_gt0;
  set perc_methyl;
  if flag_methyl > 0 then flag_methylated=1;
  else flag_methylated=0;
   condition=catx("_",dose,units); 
run;


proc sort data=flag_perc_methyl_gt0;
   by chr pos  site_type condition;
proc transpose data=flag_perc_methyl_gt0 out=flag_methyl_sbys;
   by chr pos site_type;
   id condition;
   var flag_methylated;
run;

/* Make permenant */

data wgbs.flag_methylation_gt0;
   set flag_methyl_sbys;
   drop _NAME_;
   rename _01Gy_0U = flag_01Gy_0U_methylation_gt0
          _01Gy_100U = flag_01Gy_100U_methylation_gt0
          _0Gy_0U = flag_0Gy_0U_methylation_gt0
          _0Gy_100U = flag_0Gy_100U_methylation_gt0
          _1Gy_0U = flag_1Gy_0U_methylation_gt0
          _1Gy_100U = flag_1Gy_100U_methylation_gt0;
run;

