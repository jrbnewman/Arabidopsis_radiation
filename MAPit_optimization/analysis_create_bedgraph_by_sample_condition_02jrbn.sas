libname wgbs '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname wgbsloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';

/* Output filtered bedGraph-like files 

For CG, CHG, CHH sites:
(1) Flag sites <10 coverage in any sample
(2) Create by-sample BEDs, avg-over-condition BEDs, difference BEDs (from perc methyl)
    for each
        -> Flag-filtered and unfiltered
(3) Convert to bigWig

For GC sites:
(1) Flag sites <10 coverage in any sample
(2) Flag sites where 0U > 100U
(3) Create by-sample BEDs, condition average BEDs
    Diffs between 100U and 0U by dose
    Diffs between doses by units
    Diff of diffs
        -> Flag-filtered and unfiltered
(4) Convert all to bigWig

*/
proc datasets lib=work kill noprint;
run; 
quit;

data meth_data;
  set wgbsloc.methylation_data;
  if total_C < 10 then flag_coverage_ge10=0; else flag_coverage_ge10=1;
  keep condition chr pos_end perc_methyl site_type flag_coverage_ge10;
  rename pos_end=pos ;
run;

proc sort data=meth_data;
  by chr pos site_type condition;
run;

proc transpose data=meth_data out=meth_data_sbys;
  by chr pos site_type;
  id condition;
  var perc_methyl;
run;


data sites_native;
  set wgbsloc.dmc_results_combined_annot;
  keep site_type chr pos;
run;

data sites_gc;
  set wgbsloc.dac_results_combined_annot_noGT;
  keep site_type chr pos;
run;

data all_tested_sites;
  set sites_native sites_gc;
run;


proc sort data=all_tested_sites;
  by site_type chr pos ;
proc sort data=meth_data_sbys;
  by site_type chr pos;
run;

data meth_data_sbys2;
  merge meth_data_sbys (in=in1) all_tested_sites (in=in2);
  by site_type chr pos;
  pos2=pos-2;
  pos3=pos-1;
  if in2 then flag_analyzed=1; else flag_analyzed=0;
  if in1 then output;
run;

/* For each create a "all sites" and "analyzed sites" for each condition, plus differences between doses and units */

data cg_22C_0U_all cg_22C_0U_tested
     chg_22C_0U_all chg_22C_0U_tested
     chh_22C_0U_all chh_22C_0U_tested
     gc_22C_0U_all gc_22C_0U_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_0U=. then delete;
   if site_type="CG" then output cg_22C_0U_all;
   if site_type="CG" and flag_analyzed=1 then output cg_22C_0U_tested;
   if site_type="CHG" then output chg_22C_0U_all;
   if site_type="CHG" and flag_analyzed=1 then output chg_22C_0U_tested;
   if site_type="CHH" then output chh_22C_0U_all;
   if site_type="CHH" and flag_analyzed=1 then output chh_22C_0U_tested;
   if site_type="GC" then output gc_22C_0U_all;
   if site_type="GC" and flag_analyzed=1 then output gc_22C_0U_tested;
   keep chr pos2 pos3 _22C_0U;
run;

data cg_4C_0U_all cg_4C_0U_tested
     chg_4C_0U_all chg_4C_0U_tested
     chh_4C_0U_all chh_4C_0U_tested
     gc_4C_0U_all gc_4C_0U_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _4C_0U=. then delete;
   if site_type="CG" then output cg_4C_0U_all;
   if site_type="CG" and flag_analyzed=1 then output cg_4C_0U_tested;
   if site_type="CHG" then output chg_4C_0U_all;
   if site_type="CHG" and flag_analyzed=1 then output chg_4C_0U_tested;
   if site_type="CHH" then output chh_4C_0U_all;
   if site_type="CHH" and flag_analyzed=1 then output chh_4C_0U_tested;
   if site_type="GC" then output gc_4C_0U_all;
   if site_type="GC" and flag_analyzed=1 then output gc_4C_0U_tested;
   keep chr pos2 pos3 _4C_0U;
run;

data gc_22C_100U_all gc_22C_100U_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. then delete;
   if site_type="GC" then output gc_22C_100U_all;
   if site_type="GC" and flag_analyzed=1 then output gc_22C_100U_tested;
   keep chr pos2 pos3 _22C_100U;
run;

data gc_4C_100U_all gc_4C_100U_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _4C_100U=. then delete;
   if site_type="GC" then output gc_4C_100U_all;
   if site_type="GC" and flag_analyzed=1 then output gc_4C_100U_tested;
   keep chr pos3 pos2 _4C_100U;
run;

data cg_4C_22C_0U_all cg_4C_22C_0U_tested
     chg_4C_22C_0U_all chg_4C_22C_0U_tested
     chh_4C_22C_0U_all chh_4C_22C_0U_tested;
     retain chr pos2 pos3;
     set meth_data_sbys2;
     _4C_22C_0U=_4C_0U-_22C_0U;
     if _22C_0U=. or _4C_0U=. then delete;
     if site_type="CG" then output cg_4C_22C_0U_all;
     if site_type="CG" and flag_analyzed=1 then output cg_4C_22C_0U_tested;
     if site_type="CHG" then output chg_4C_22C_0U_all;
     if site_type="CHG" and flag_analyzed=1 then output chg_4C_22C_0U_tested;
     if site_type="CHH" then output chh_4C_22C_0U_all;
     if site_type="CHH" and flag_analyzed=1 then output chh_4C_22C_0U_tested;
     keep chr pos3 pos2 _4C_22C_0U;
run;

data gc_fans_0p5_all gc_fans_0p5_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2; 
   if FANS_0p5U=. then delete;
   if site_type="GC" then output gc_fans_0p5_all;
   if site_type="GC" and flag_analyzed=1 then output gc_fans_0p5_tested;
   keep chr pos3 pos2 FANS_0p5U;
run;

data gc_fans_1p5_all gc_fans_1p5_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2; 
   if FANS_1p5U=. then delete;
   if site_type="GC" then output gc_fans_1p5_all;
   if site_type="GC" and flag_analyzed=1 then output gc_fans_1p5_tested;
   keep chr pos3 pos2 FANS_1p5U;
run;

data gc_fans_5_all gc_fans_5_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2; 
   if FANS_5U=. then delete;
   if site_type="GC" then output gc_fans_5_all;
   if site_type="GC" and flag_analyzed=1 then output gc_fans_5_tested;
   keep chr pos3 pos2 FANS_5U;
run;

data gc_fans_25_all gc_fans_25_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2; 
   if FANS_25U=. then delete;
   if site_type="GC" then output gc_fans_25_all;
   if site_type="GC" and flag_analyzed=1 then output gc_fans_25_tested;
   keep chr pos3 pos2 FANS_25U;
run;

data gc_22C_100U_0U_all gc_22C_100U_0U_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. or _22C_0U=. then delete;
   _22C_100U_0U=_22C_100U-_22C_0U;
   if site_type="GC" then output gc_22C_100U_0U_all;
   if site_type="GC" and flag_analyzed=1 then output gc_22C_100U_0U_tested;
   keep chr pos3 pos2 _22C_100U_0U;
run;

data gc_4C_100U_0U_all gc_4C_100U_0U_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _4C_100U=. or _4C_0U=. then delete;
   _4C_100U_0U=_4C_100U-_4C_0U;
   if site_type="GC" then output gc_4C_100U_0U_all;
   if site_type="GC" and flag_analyzed=1 then output gc_4C_100U_0U_tested;
   keep chr pos3 pos2 _4C_100U_0U;
run;

data gc_4C_22C_all gc_4C_22C_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. or _22C_0U=. or _4C_100U=. or _4C_0U=. then delete;
   _22C_4C_100U_0U=(_4C_100U-_4C_0U)-(_22C_100U-_22C_0U);
   if site_type="GC" then output gc_4C_22C_all;
   if site_type="GC" and flag_analyzed=1 then output gc_4C_22C_tested;
   keep chr pos3 pos2 _22C_4C_100U_0U;
run;

data gc_FANS_0p5_22C_all gc_FANS_0p5_22C_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. or FANS_0p5U=. then delete;
   FANS_0p5U_22C_100U=FANS_0p5U-_22C_100U;
   if site_type="GC" then output gc_FANS_0p5_22C_all;
   if site_type="GC" and flag_analyzed=1 then output gc_FANS_0p5_22C_tested;
   keep chr pos3 pos2 FANS_0p5U_22C_100U;
run;

data gc_FANS_1p5_22C_all gc_FANS_1p5_22C_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. or FANS_1p5U=. then delete;
   FANS_1p5U_22C_100U=FANS_1p5U-_22C_100U;
   if site_type="GC" then output gc_FANS_1p5_22C_all;
   if site_type="GC" and flag_analyzed=1 then output gc_FANS_1p5_22C_tested;
   keep chr pos3 pos2 FANS_1p5U_22C_100U;
run;

data gc_FANS_5_22C_all gc_FANS_5_22C_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. or FANS_5U=. then delete;
   FANS_5U_22C_100U=FANS_5U-_22C_100U;
   if site_type="GC" then output gc_FANS_5_22C_all;
   if site_type="GC" and flag_analyzed=1 then output gc_FANS_5_22C_tested;
   keep chr pos3 pos2 FANS_5U_22C_100U;
run;


data gc_FANS_25_22C_all gc_FANS_25_22C_tested;
   retain chr pos2 pos3;
   set meth_data_sbys2;
   if _22C_100U=. or FANS_25U=. then delete;
   FANS_25U_22C_100U=FANS_25U-_22C_100U;
   if site_type="GC" then output gc_FANS_25_22C_all;
   if site_type="GC" and flag_analyzed=1 then output gc_FANS_25_22C_tested;
   keep chr pos3 pos2 FANS_25U_22C_100U;
run;

%macro exportBED(site,test);

proc export data=&site._&test._all
    outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/bedGraphs/new/&site._&test._all_new.bedGraph"
    dbms=tab replace;
    putnames=no;
run;

proc export data=&site._&test._tested
    outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/bedGraphs/new/&site._&test._analyzed_new.bedGraph"
    dbms=tab replace;
    putnames=no;
run;


%mend;

%exportBED(CG,22C_0U); %exportBED(CHG,22C_0U); %exportBED(CHH,22C_0U); %exportBED(GC,22C_0U);
%exportBED(CG,4C_0U); %exportBED(CHG,4C_0U); %exportBED(CHH,4C_0U); %exportBED(GC,4C_0U);
%exportBED(GC,22C_100U); %exportBED(GC,4C_100U);
%exportBED(CG,4C_22C_0U); %exportBED(CHG,4C_22C_0U); %exportBED(CHH,4C_22C_0U); 
%exportBED(GC,FANS_0p5); %exportBED(GC,FANS_1p5); %exportBED(GC,FANS_5); %exportBED(GC,FANS_25);
%exportBED(GC,22C_100U_0U); %exportBED(GC,4C_100U_0U); %exportBED(GC,4C_22C);
%exportBED(GC,FANS_0p5_22C); %exportBED(GC,FANS_1p5_22C); %exportBED(GC,FANS_5_22C); %exportBED(GC,FANS_25_22C); 

