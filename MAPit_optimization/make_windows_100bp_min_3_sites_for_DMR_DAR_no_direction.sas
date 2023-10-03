/* Group cytosines by sequence and chromosome DMCs into regions/island of the same direction */

ods html close;
ods listing;
libname cold '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname coldloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';

proc datasets lib=work kill noprint;
run;
quit;

/* Group DMCs into windows
   Sites must be no more than 100bp away from the next site
   Sites must have the same direction to be included
   Region must have at least 3 sites
   A DMR has at least 3 DMCs (flag 1 and 3)
 */

data methyl_count;
  set coldloc.methylation_data;
  if total_C < 10 then delete;
  keep site_type chr pos_end condition perc_methyl;
  rename pos_end=pos;
run;

proc sort data=methyl_count;
  by site_type chr pos condition;
proc transpose data=methyl_count out=methyl_sbys;
 by  site_type chr pos;
  id condition;
  var perc_methyl;
run;

data methyl_meth;
  set methyl_sbys;
  where site_type ne "GC";
  diff_4C_22C = _4C_0U - _22C_0U;
  if _4C_0U=0 and _22C_0U=0 then delete;
  if _4C_0U=. or _22C_0U=. then delete;
  keep site_type chr pos diff_4C_22C _4C_0U _22C_0U;
  rename _4C_0U=meth_4C_0U _22C_0U=meth_22C_0U  ;
run;


data methyl_access;
  set methyl_sbys;
  where site_type = "GC";
  diff_22C_100U_0U=_22C_100U - _22C_0U;
  diff_4C_100U_0U=_4C_100U - _4C_0U;
  diff_4C_22C = diff_4C_100U_0U - diff_22C_100U_0U;
  if _4C_0U=0 and _22C_0U=0 and _4C_100U=0 and _22C_100U=0 then delete;
  if _4C_0U=. or _22C_0U=. or _4C_100U=. or _22C_100U=. then delete;
  keep site_type chr pos diff_4C_22C diff_4C_100U_0U diff_22C_100U_0U 
       _4C_0U _22C_0U _4C_100U _22C_100U;
  rename _4C_0U=meth_4C_0U _22C_0U=meth_22C_0U
         _4C_100U=meth_4C_100U _22C_100U=meth_22C_100U  ;
run;



data dmc;
  set methyl_meth;
  if diff_4C_22C=. then delete;
  keep chr pos site_type diff_4C_22C ;
run;

proc sort data=dmc;
  by site_type chr pos;
run;

data group_dmc;
  retain region_num;
  set dmc;
  by site_type chr;
  prev_pos=lag1(pos);
  if first.chr then region_num=1;
  else do;
     if prev_pos >= pos-100 then region_num=region_num;
     else region_num=region_num + 1;
     end;
run;


proc freq data=group_dmc noprint;
  tables site_type*chr*region_num / out=sites_per_region;
run;

data flag_region_ge3;
  set sites_per_region;
  if count >= 3 then flag_num_sites_ge3=1;
  else flag_num_sites_ge3=0;
run;

proc freq data=flag_region_ge3;
  tables site_type*flaG_num_sites_ge3;
run;

/*
 site_type     flag_num_sites_ge3

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
 CG       | 167948 | 120755 | 288703
          |  19.72 |  14.18 |  33.90
          |  58.17 |  41.83 |
          |  35.88 |  31.47 |
 ---------+--------+--------+
 CHG      | 199590 |  81628 | 281218
          |  23.43 |   9.58 |  33.02
          |  70.97 |  29.03 |
          |  42.64 |  21.27 |
 ---------+--------+--------+
 CHH      | 100515 | 181306 | 281821
          |  11.80 |  21.29 |  33.09
          |  35.67 |  64.33 |
          |  21.48 |  47.25 |
 ---------+--------+--------+
 Total      468053   383689   851742
             54.95    45.05   100.00



*/

/* Make permanent */

data cold.cytosine_2_meth_region_index_nd2;
  set group_dmc;
run;



/* DARs */

data dac;
  set methyl_access;
  if diff_22C_100U_0U=. or  diff_4C_100U_0U=. then delete;
  if diff_22C_100U_0U < 0 then flag_22C_0U_gt_100U=1; else flag_22C_0U_gt_100U=0;
  if diff_4C_100U_0U < 0 then flag_4C_0U_gt_100U=1; else flag_4C_0U_gt_100U=0;
  keep site_type chr pos 
  diff_22C_100U_0U diff_4C_100U_0U diff_4C_22C flag_22C_0U_gt_100U flag_4C_0U_gt_100U;
run;

proc sort data=dac;
  by site_type chr pos;
run;

data group_dac;
  retain region_num;
  set dac;
  by site_type  chr;
  prev_pos=lag1(pos);
  if first.chr then region_num=1;
  else do;
     if prev_pos >= pos-100 then region_num=region_num;
     else region_num=region_num + 1;
     end;
run;

proc freq data=group_dac noprint;
  tables site_type*chr*region_num / out=sites_per_region;
run;

data flag_region_ge3;
  set sites_per_region;
  if count >= 3 then flag_num_sites_ge3=1;
  else flag_num_sites_ge3=0;
run;

proc freq data=flag_region_ge3;
  tables site_type*flaG_num_sites_ge3;
run;


/*
Table of site_type by flag_num_sites_ge3


site_type     flag_num_sites_ge3

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
GC       |  53954 | 192610 | 246564
         |  21.88 |  78.12 | 100.00
         |  21.88 |  78.12 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total       53954   192610   246564
            21.88    78.12   100.00

*/



/* Make permanent */

data cold.cytosine_2_acc_region_index_nd2;
  set group_dac;
run;

