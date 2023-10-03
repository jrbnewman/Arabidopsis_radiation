libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Output filtered bedGraph-like files for sites with 10X coverage: any that meet criteria, and those that where analyzed 
    by condition, between conditions */

data meth_data;
  set wgbslocA.methylation_data_cg_chg_chh wgbslocA.methylation_data_gc;
  start_pos=start_pos-1;
  stop_pos=stop_pos-1;
  if site_type="GC" then perc_methyl=perc_methyl_norm;
run;

proc sort data=meth_data;
  by site_type chr start_pos stop_pos treatment units ;
proc means data=meth_data noprint;
  by site_Type chr start_pos stop_pos treatment units;
  var total_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_c mean(perc_methyl)=perc_methyl;
run;

data meth_data3;
  set meth_data2;
  if total_C < 10 and perc_methyl ne . then delete;
run;



proc sort data=meth_data3;
   by site_type chr start_pos stop_pos treatment units;
proc transpose data=meth_data3 out=meth_sbys;
   by site_type chr start_pos stop_pos ;
  id treatment units;
  var perc_methyl;
run;


data native_sites_01 native_sites_1;
    set wgbslocA.results_by_dmc_annot;
  start_pos=start_pos-1;
  stop_pos=stop_pos-1;
   if comparison="0Gy_vs_01G" then output native_sites_01;
   if comparison="0Gy_vs_1Gy" then output native_sites_1;
    keep site_type chr start_pos stop_pos;
run;


data gc_sites_01 gc_sites_1;
    set wgbslocA.results_by_dac_annot;
  start_pos=start_pos-1;
  stop_pos=stop_pos-1;
   if comparison="01Gy_0Gy" then output gc_sites_01;
   if comparison="1Gy_0Gy" then output gc_sites_1;
    keep site_type chr start_pos stop_pos;
run;



data all_tested_sites_01;
  set native_sites_01 gc_sites_01;
run;


data all_tested_sites_1;
  set native_sites_1 gc_sites_1;
run;


proc sort data=meth_sbys;
  by site_type chr start_pos stop_pos;
proc sort data=all_tested_sites_01 nodup;
  by site_type chr start_pos stop_pos;
proc sort data=all_tested_sites_1 nodup;
  by site_type chr start_pos stop_pos;
run;

data meth_data_sbys2;
  merge meth_sbys (in=in1) all_tested_sites_01 (in=in2) all_tested_sites_1 (in=in3);
  by site_type chr start_pos stop_pos;
  if in2 then flag_analyzed_01=1; else flag_analyzed_01=0;
  if in3 then flag_analyzed_1=1; else flag_analyzed_1=0;
  if in1 then output;
run;

/* For each create a "all sites" and "analyzed sites" for each condition, plus differences between doses and units */


proc means data=meth_data_sbys2 noprint;
   var _numeric_;
   output out=check_meth max=;
run;

proc means data=meth_data_sbys2 noprint;
   var _numeric_;
   output out=check_meth min=;
run;


%macro makeBedgraph(siteType,condit,filtVar);

data &siteType.&condit._all &siteType.&condit._tested ;
     retain chr start_pos stop_pos;
     set meth_data_sbys2;
     if &condit. = . then delete;
     if site_type="&siteType." then output &siteType&condit._all;
     %if &filtVar. = none %then %do;
     if site_type="&siteType." and (flag_analyzed_01=1 or flag_analyzed_1=1) then output &siteType.&condit._tested;
     %end; %else %do;   
     if site_type="&siteType." and &filtVar.=1 then output &siteType.&condit._tested;
     %end; 
     keep chr start_pos stop_pos &condit.;
run;

proc sort data=&siteType.&condit._all ;
   by chr start_pos stop_pos;
proc sort data=&siteType.&condit._tested ;
   by chr start_pos stop_pos;
run;


%mend;

%makeBedgraph(CG,_0Gy0U,none);  %makeBedgraph(CG,_01Gy0U,flag_analyzed_01); %makeBedgraph(CG,_1Gy0U,flag_analyzed_1); 
%makeBedgraph(CHG,_0Gy0U,none);  %makeBedgraph(CHG,_01Gy0U,flag_analyzed_01); %makeBedgraph(CHG,_1Gy0U,flag_analyzed_1); 
%makeBedgraph(CHH,_0Gy0U,none);  %makeBedgraph(CHH,_01Gy0U,flag_analyzed_01); %makeBedgraph(CHH,_1Gy0U,flag_analyzed_1); 
%makeBedgraph(GC,_0Gy0U,none);  %makeBedgraph(GC,_01Gy0U,flag_analyzed_01); %makeBedgraph(GC,_1Gy0U,flag_analyzed_1); 
%makeBedgraph(GC,_0Gy100U,none);  %makeBedgraph(GC,_01Gy100U,flag_analyzed_01); %makeBedgraph(GC,_1Gy100U,flag_analyzed_1); 




%macro makeCompare(siteType,condit1,condit2,outName,filtVar);

data &outName._all &outName._tested;
   retain chr start_pos stop_pos;
    set meth_data_sbys2;
    where site_type="&siteType." ;
    if &condit1. = . or &condit2. = . then delete;
    meth_diff = &condit1. - &condit2.;
    output &outName._all;
     %if &filtVar.=none %then %do;
     if site_type="&siteType." and (flag_analyzed_01=1 or flag_analyzed_1=1) then output &outName._tested;
     %end;
     %else %do;   
     if site_type="&siteType." and &filtVar.=1 then output &outName._tested;
     %end; 

    keep chr start_pos stop_pos meth_diff;
run;

%mend;



%makeCompare(CG,_01Gy0U,_0Gy0U,CG_01Gy_0Gy,flag_analyzed_01);
%makeCompare(CG,_1Gy0U,_0Gy0U,CG_1Gy_0Gy,flag_analyzed_1);
%makeCompare(CHG,_01Gy0U,_0Gy0U,CHG_01Gy_0Gy,flag_analyzed_01);
%makeCompare(CHG,_1Gy0U,_0Gy0U,CHG_1Gy_0Gy,flag_analyzed_1);
%makeCompare(CHH,_01Gy0U,_0Gy0U,CHH_01Gy_0Gy,flag_analyzed_01);
%makeCompare(CHH,_1Gy0U,_0Gy0U,CHH_1Gy_0Gy,flag_analyzed_1);
%makeCompare(GC,_0Gy100U,_0Gy0U,GC_0Gy_100U_0U,none);
%makeCompare(GC,_01Gy100U,_01Gy0U,GC_01Gy_100U_0U,flag_analyzed_01);
%makeCompare(GC,_1Gy100U,_1Gy0U,GC_1Gy_100U_0U,flag_analyzed_1);




%macro makeGCCompare(siteType,condit1,condit2,condit3,condit4,outName,filtVar);

data &outName._all &outName._tested;
   retain chr start_pos stop_pos;
    set meth_data_sbys2;
    where site_type="&siteType." ;
    if &condit1. = . or &condit2. = . or &condit3. = . or &condit4. = . then delete;
    meth_diff = (&condit1. - &condit2.) -  (&condit3. - &condit4.);
    output &outName._all;
     if site_type="&siteType." and &filtVar.=1 then output &outName._tested;
    keep chr start_pos stop_pos meth_diff;
run;

%mend;


%makeGCCompare(GC,_01Gy100U,_01Gy0U,_0Gy100U,_0Gy0U,GC_01Gy_0Gy,flag_analyzed_01);
%makeGCCompare(GC,_1Gy100U,_1Gy0U,_0Gy100U,_0Gy0U,GC_1Gy_0Gy,flag_analyzed_1);



%macro exportBED(dataSet,outFile);

proc export data=&dataSet.
    outfile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_&outFile..bedGraph"
    dbms=tab replace;
    putnames=no;
run;

%mend;


%exportBED(cg_0Gy0U_all,CG_0Gy_all);
%exportBED(cg_0Gy0U_tested,CG_0Gy_analyzed);
%exportBED(chg_0Gy0U_all,CHG_0Gy_all);
%exportBED(chg_0Gy0U_tested,CHG_0Gy_analyzed);
%exportBED(chh_0Gy0U_all,CHH_0Gy_all);
%exportBED(chh_0Gy0U_tested,CHH_0Gy_analyzed);
%exportBED(gc_0Gy0U_all,GC_0Gy_0U_all);
%exportBED(gc_0Gy0U_tested,GC_0Gy_0U_analyzed);
%exportBED(gc_0Gy100U_all,GC_0Gy_100U_all);
%exportBED(gc_0Gy100U_tested,GC_0Gy_100U_analyzed);

%exportBED(cg_01Gy0U_all,CG_01Gy_all);
%exportBED(cg_01Gy0U_tested,CG_01Gy_analyzed);
%exportBED(chg_01Gy0U_all,CHG_01Gy_all);
%exportBED(chg_01Gy0U_tested,CHG_01Gy_analyzed);
%exportBED(chh_01Gy0U_all,CHH_01Gy_all);
%exportBED(chh_01Gy0U_tested,CHH_01Gy_analyzed);
%exportBED(gc_01Gy0U_all,GC_01Gy_0U_all);
%exportBED(gc_01Gy0U_tested,GC_01Gy_0U_analyzed);
%exportBED(gc_01Gy100U_all,GC_01Gy_100U_all);
%exportBED(gc_01Gy100U_tested,GC_01Gy_100U_analyzed);

%exportBED(cg_1Gy0U_all,CG_1Gy_all);
%exportBED(cg_1Gy0U_tested,CG_1Gy_analyzed);
%exportBED(chg_1Gy0U_all,CHG_1Gy_all);
%exportBED(chg_1Gy0U_tested,CHG_1Gy_analyzed);
%exportBED(chh_1Gy0U_all,CHH_1Gy_all);
%exportBED(chh_1Gy0U_tested,CHH_1Gy_analyzed);
%exportBED(gc_1Gy0U_all,GC_1Gy_0U_all);
%exportBED(gc_1Gy0U_tested,GC_1Gy_0U_analyzed);
%exportBED(gc_1Gy100U_all,GC_1Gy_100U_all);
%exportBED(gc_1Gy100U_tested,GC_1Gy_100U_analyzed);

%exportBED(gc_0Gy_100U_0U_all,GC_0Gy_100U_0U_all);
%exportBED(gc_0Gy_100U_0U_tested,GC_0Gy_100U_0U_analyzed);
%exportBED(gc_01Gy_100U_0U_all,GC_01Gy_100U_0U_all);
%exportBED(gc_01Gy_100U_0U_tested,GC_01Gy_100U_0U_analyzed);
%exportBED(gc_1Gy_100U_0U_all,GC_1Gy_100U_0U_all);
%exportBED(gc_1Gy_100U_0U_tested,GC_1Gy_100U_0U_analyzed);


%exportBED(gc_01Gy_0gy_all,GC_01Gy_0Gy_all);
%exportBED(gc_01Gy_0gy_tested,GC_01Gy_0Gy_analyzed);
%exportBED(gc_1Gy_0gy_all,GC_1Gy_0Gy_all);
%exportBED(gc_1Gy_0gy_tested,GC_1Gy_0Gy_analyzed);

%exportBED(cg_01Gy_0gy_all,CG_01Gy_0Gy_all);
%exportBED(cg_01Gy_0gy_tested,CG_01Gy_0Gy_analyzed);
%exportBED(cg_1Gy_0gy_all,CG_1Gy_0Gy_all);
%exportBED(cg_1Gy_0gy_tested,CG_1Gy_0Gy_analyzed);


%exportBED(chg_01Gy_0gy_all,CHG_01Gy_0Gy_all);
%exportBED(chg_01Gy_0gy_tested,CHG_01Gy_0Gy_analyzed);
%exportBED(chg_1Gy_0gy_all,CHG_1Gy_0Gy_all);
%exportBED(chg_1Gy_0gy_tested,CHG_1Gy_0Gy_analyzed);


%exportBED(chh_01Gy_0gy_all,CHH_01Gy_0Gy_all);
%exportBED(chh_01Gy_0gy_tested,CHH_01Gy_0Gy_analyzed);
%exportBED(chh_1Gy_0gy_all,CHH_1Gy_0Gy_all);
%exportBED(chh_1Gy_0gy_tested,CHH_1Gy_0Gy_analyzed);






%macro findMax(dataSet,outFile);

proc means data=&dataSet.;
 var _numeric_;
  output max=;
run;


%mend;


data check;
   set gc_01gy_0gy_all;
   where abs(meth_diff) > 1;
run;



%findMax(cg_0Gy0U_all,CG_0Gy_all);
%findMax(cg_0Gy0U_tested,CG_0Gy_analyzed);
%findMax(chg_0Gy0U_all,CHG_0Gy_all);
%findMax(chg_0Gy0U_tested,CHG_0Gy_analyzed);
%findMax(chh_0Gy0U_all,CHH_0Gy_all);
%findMax(chh_0Gy0U_tested,CHH_0Gy_analyzed);
%findMax(gc_0Gy0U_all,GC_0Gy_0U_all);
%findMax(gc_0Gy0U_tested,GC_0Gy_0U_analyzed);
%findMax(gc_0Gy100U_all,GC_0Gy_100U_all);
%findMax(gc_0Gy100U_tested,GC_0Gy_100U_analyzed);
%findMax(cg_01Gy0U_all,CG_01Gy_all);
%findMax(cg_01Gy0U_tested,CG_01Gy_analyzed);
%findMax(chg_01Gy0U_all,CHG_01Gy_all);
%findMax(chg_01Gy0U_tested,CHG_01Gy_analyzed);
%findMax(chh_01Gy0U_all,CHH_01Gy_all);
%findMax(chh_01Gy0U_tested,CHH_01Gy_analyzed);
%findMax(gc_01Gy0U_all,GC_01Gy_0U_all);
%findMax(gc_01Gy0U_tested,GC_01Gy_0U_analyzed);
%findMax(gc_01Gy100U_all,GC_01Gy_100U_all);
%findMax(gc_01Gy100U_tested,GC_01Gy_100U_analyzed);
%findMax(cg_1Gy0U_all,CG_1Gy_all);
%findMax(cg_1Gy0U_tested,CG_1Gy_analyzed);
%findMax(chg_1Gy0U_all,CHG_1Gy_all);
%findMax(chg_1Gy0U_tested,CHG_1Gy_analyzed);
%findMax(chh_1Gy0U_all,CHH_1Gy_all);
%findMax(chh_1Gy0U_tested,CHH_1Gy_analyzed);
%findMax(gc_1Gy0U_all,GC_1Gy_0U_all);
%findMax(gc_1Gy0U_tested,GC_1Gy_0U_analyzed);
%findMax(gc_1Gy100U_all,GC_1Gy_100U_all);
%findMax(gc_1Gy100U_tested,GC_1Gy_100U_analyzed);
%findMax(gc_0Gy_100U_0U_all,GC_0Gy_100U_0U_all);
%findMax(gc_0Gy_100U_0U_tested,GC_0Gy_100U_0U_analyzed);
%findMax(gc_01Gy_100U_0U_all,GC_01Gy_100U_0U_all);
%findMax(gc_01Gy_100U_0U_tested,GC_01Gy_100U_0U_analyzed);
%findMax(gc_1Gy_100U_0U_all,GC_1Gy_100U_0U_all);
%findMax(gc_1Gy_100U_0U_tested,GC_1Gy_100U_0U_analyzed);
%findMax(gc_01Gy_0gy_all,GC_01Gy_0Gy_all);
%findMax(gc_01Gy_0gy_tested,GC_01Gy_0Gy_analyzed);
%findMax(gc_1Gy_0gy_all,GC_1Gy_0Gy_all);
%findMax(gc_1Gy_0gy_tested,GC_1Gy_0Gy_analyzed);
%findMax(cg_01Gy_0gy_all,CG_01Gy_0Gy_all);
%findMax(cg_01Gy_0gy_tested,CG_01Gy_0Gy_analyzed);
%findMax(cg_1Gy_0gy_all,CG_1Gy_0Gy_all);
%findMax(cg_1Gy_0gy_tested,CG_1Gy_0Gy_analyzed);
%findMax(chg_01Gy_0gy_all,CHG_01Gy_0Gy_all);
%findMax(chg_01Gy_0gy_tested,CHG_01Gy_0Gy_analyzed);
%findMax(chg_1Gy_0gy_all,CHG_1Gy_0Gy_all);
%findMax(chg_1Gy_0gy_tested,CHG_1Gy_0Gy_analyzed);
%findMax(chh_01Gy_0gy_all,CHH_01Gy_0Gy_all);
%findMax(chh_01Gy_0gy_tested,CHH_01Gy_0Gy_analyzed);
%findMax(chh_1Gy_0gy_all,CHH_1Gy_0Gy_all);
%findMax(chh_1Gy_0gy_tested,CHH_1Gy_0Gy_analyzed);

