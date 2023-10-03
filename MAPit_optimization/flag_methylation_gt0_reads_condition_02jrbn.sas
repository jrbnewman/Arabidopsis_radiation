ods listing; ods html close;
libname wgbs '!PATCON/arabidopsis_wgbs/sas_data';
libname wgbsloc "/TB14/TB14/sandbox/wgbs_sandbox/sas_data";

/* Flag sites that have at least 1 mapped read in at least one condition
For analysis, we will pick the sites that have at least a total of 10 mapped reads at 22C AND 4C
(0U for native, BOTH for 100U), and have at least one methylated read at 22C or 4C

For 22C 100U vs FANS, flag for FANS analysis
*/


data flag_coverage;
  set wgbs.flag_coverage_ge10;
  if site_type="GC" then do;
     if flag_22C_0U_coverage_ge10=1 and flag_22C_100U_coverage_ge10=1 then flag_22C_ge10=1; else flag_22C_ge10=0;
     if flag_4C_0U_coverage_ge10=1 and flag_4C_100U_coverage_ge10=1 then flag_4C_ge10=1; else flag_4C_ge10=0;
     if (flag_FANS_0p5U_coverage_ge10=1 or flag_FANS_1p5U_coverage_ge10=1 or
         flag_FANS_5U_coverage_ge10=1 or flag_FANS_25U_coverage_ge10=1) then flag_FANS_ge10=1; else flag_FANS_ge10=0;
     end;
   else do;
        if flag_22C_0U_coverage_ge10=1 then flag_22C_ge10=1; else flag_22C_ge10=0;
        if  flag_4C_0U_coverage_ge10=1 then flag_4C_ge10=1; else flag_4C_ge10=0;
   end;
run;

data sites2keep sites2keep_fans;
  set flag_coverage;
  if flag_22C_ge10=1 and flag_4C_ge10=1 then output sites2keep;
  if site_type="GC" and flag_22C_ge10=1 and flag_FANS_ge10=1 then output sites2keep_fans;
  keep chr pos site_type ;
run;

data  meth_data;
  set wgbsloc.methylation_data;
  if chr="" then delete;
  if temperature="FANS" then delete;
  keep condition temperature units site_type chr pos_end methyl_C total_C ;
  rename pos_end=pos;
run;

data  meth_data_fans;
  set wgbsloc.methylation_data;
  if chr="" then delete;
  if temperature="FANS" then output;
  else if temperature="22C" and units="100U" then output;
  keep condition temperature units site_type chr pos_end methyl_C total_C ;
  rename pos_end=pos;
run;


proc sort data=meth_data;
  by chr pos site_type;
proc sort data=meth_data_fans;
  by chr pos site_type;
proc sort data=sites2keep;
  by chr pos site_type;
proc sort data=sites2keep_fans;
  by chr pos site_type;
run;

data meth_data2;
  merge sites2keep (in=in1) meth_data (in=in2);
  by chr pos site_type;
  if in1 and in2;
run;

data meth_data_fans2;
  merge sites2keep_fans (in=in1) meth_data_fans (in=in2);
  by chr pos site_type;
  if in1 and in2;
run;

data flag_methylation;
  set meth_data2;
  if total_C < 10 then delete; *only want to count sites if there are at least 10 reads;
  if methyl_C > 1 then flag_methyl=1;
  else flag_methyl=0;
  run;

data flag_methylation_fans;
  set meth_data_fans2;
  if total_C < 10 then delete; *only want to count sites if there are at least 10 reads;
  if methyl_C > 1 then flag_methyl=1;
  else flag_methyl=0;
  run;

proc sort data=flag_methylation;
   by chr pos  site_type condition;
proc transpose data=flag_methylation out=flag_methyl_sbys;
   by chr pos site_type;
   id condition;
   var flag_methyl;
run;

proc sort data=flag_methylation_fans;
   by chr pos  site_type condition;
proc transpose data=flag_methylation_fans out=flag_methyl_sbys_fans;
   by chr pos site_type;
   id condition;
   var flag_methyl;
run;

/* Make permenant */

data wgbs.flag_methylation_gt0;
   set flag_methyl_sbys;
   drop _NAME_;
   rename _22C_0U = flag_22C_0U_methylation_gt0
          _4C_0U = flag_4C_0U_methylation_gt0
          _22C_100U = flag_22C_100U_methylation_gt0
          _4C_100U = flag_4C_100U_methylation_gt0;
run;

data wgbs.flag_methylation_FANS_gt0;
   set flag_methyl_sbys_fans;
   drop _NAME_;
   rename _22C_100U = flag_22C_100U_methylation_gt0
          FANS_0p5U = flag_FANS_0p5U_methylation_gt0
          FANS_1p5U = flag_FANS_1p5U_methylation_gt0
          FANS_25U = flag_FANS_25U_methylation_gt0
          FANS_5U = flag_FANS_5U_methylation_gt0;
run;


