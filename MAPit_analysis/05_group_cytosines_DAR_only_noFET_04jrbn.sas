libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;



/* Group DMCs into windows
 Sites must be no more than 100bp away from the next site
 Sites must have the same direction to be included
 Region must have at least 3 sites
 A DMR has at least 3 DMCs (flag 1 and 3)
 */

data dac;
	set wgbslocA.results_by_dac_w_meth_v4;
if flag_meth_diff_10perc = 1 ;
keep comparison site_type chr start_pos stop_pos flag_meth_diff;
run;

proc sort data=dac;
by comparison site_type chr start_pos stop_pos;
run;

data group_dac;
retain region_num;
set dac;
by comparison site_type chr;
   prev_pos=lag1(start_pos);
   prev_direction=lag1(flag_meth_diff);
   if first.chr then region_num=1;
   else do;
   if prev_pos >= start_pos-100 and prev_direction=flag_meth_diff then region_num=region_num;
  else region_num=region_num + 1;
 end;
run;


proc freq data=group_dac noprint;
tables comparison*site_type*chr*region_num / out=sites_per_region;
run;

data flag_region_ge2;
set sites_per_region;
if count >= 2 then flag_num_sites_ge2=1;
else flag_num_sites_ge2=0;
if count >= 3 then flag_num_sites_ge3=1;
else flag_num_sites_ge3=0;
if count >= 5 then flag_num_sites_ge5=1;
else flag_num_sites_ge5=0;
if count >= 10 then flag_num_sites_ge10=1;
else flag_num_sites_ge10=0;

run;

proc freq data=flag_region_ge2;
tables comparison*site_type*flaG_num_sites_ge2 / out=dmr_check1;
tables comparison*site_type*flaG_num_sites_ge3 / out=dmr_check2;
tables comparison*site_type*flaG_num_sites_ge5 / out=dmr_check3;
tables comparison*site_type*flaG_num_sites_ge10 / out=dmr_check4;
proc print data=dmr_check1;
proc print data=dmr_check2;
proc print data=dmr_check3;
proc print data=dmr_check4;
run;

/*
10% diff
            2       3       5       10
01Gy_vs_0Gy 59635   23080   4570    135
1Gy_vs_0Gy  48610   14479   1660    18

20% diff
            2       3       5       10
01Gy_vs_0Gy 40059   13776   2038    29
1Gy_vs_0Gy  19977   4514    308     1

*/


data dac_meth_all2;
	set wgbslocA.results_by_dac_w_meth_v4;
	run;



proc sort data=group_dac;
   by comparison site_type chr start_pos stop_pos;
proc sort data=dac_meth_all2;
   by comparison site_type chr start_pos stop_pos;
run;

data group_dac2;
   merge group_dac (in=in1) dac_meth_all2 (in=in2);
   by comparison site_type chr start_pos stop_pos;
   if in1 and in2;
run;


data region2keep;
set flag_region_ge2;
where flaG_num_sites_ge2=1;
keep comparison site_type chr region_num;
run;

proc sort data=region2keep;
by comparison site_type chr region_num;
proc sort data=group_dac2;
by comparison site_type chr region_num;
run;

data group_dac3;
merge group_dac2 (in=in1) region2keep (in=in2);
by comparison site_type chr region_num;
if in1 and in2;
run;


/* Make permanent */

data wgbslocA.cytosine_to_acc_region_index_v5;
set group_dac2;
run;

data wgbslocA.acc_region_flag_site_ge2_v5;
set flag_region_ge2;
run;


data wgbslocA.cytosine_to_acc_region_ge2_v5;
set group_dac3;
run;


/* Export BED files for various tools */


data dars;
   set wgbslocA.cytosine_to_acc_region_ge2_v5;
   length regionID $100.;
   regionID = catx("_",site_type,comparison,chr,region_num);
   keep site_type comparison chr start_pos stop_pos regionID;
run;

proc sort data=dars;
  by site_type comparison chr regionID;
proc means data=dars noprint;
   by site_type comparison chr regionID;
   var start_pos stop_pos;
   output out=dar_regions min(start_pos)=start_pos max(stop_pos)=stop_pos;
run;

proc sort data=dar_regions;
  by chr start_pos stop_pos comparison regionID;
run;

data dar_gc_01 dar_gc_1;
     retain chr start_pos stop_pos regionID score strand;
   set dar_regions;
     score=0;
     strand="+";
     if site_type = "GC" and comparison = "01Gy_0Gy" then output dar_gc_01;
     if site_type = "GC" and comparison = "1Gy_0Gy" then output dar_gc_1;
     keep chr start_pos stop_pos regionID score strand;
run;



proc export data=dar_gc_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_10cGy_v2.bed" dbms=tab replace; putnames=no; run;
proc export data=dar_gc_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_100cGy_v2.bed" dbms=tab replace; putnames=no; run;




