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
  if diff_4C_22C < 0 then diff_direction=-1;
  else diff_direction=1;
  keep chr pos site_type diff_4C_22C diff_direction;
run;

proc sort data=dmc;
  by site_type chr pos;
run;

data group_dmc;
  retain region_num;
  set dmc;
  by site_type chr;
  prev_pos=lag1(pos);
  prev_direction=lag1(diff_direction);
  if first.chr then region_num=1;
  else do;
     if prev_pos >= pos-100 and prev_direction=diff_direction then region_num=region_num;
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
CG       | 750754 | 165624 | 916378
         |  19.06 |   4.21 |  23.27
         |  81.93 |  18.07 |
         |  24.49 |  18.99 |
---------+--------+--------+
CHG      | 625192 | 121781 | 746973
         |  15.88 |   3.09 |  18.97
         |  83.70 |  16.30 |
         |  20.39 |  13.96 |
---------+--------+--------+
CHH      |1689626 | 584987 |2274613
         |  42.91 |  14.86 |  57.76
         |  74.28 |  25.72 |
         |  55.12 |  67.06 |
---------+--------+--------+
Total     3065572   872392  3937964
            77.85    22.15   100.00
*/

/* Make permanent */

data cold.cytosine_to_meth_region_index2;
  set group_dmc;
run;

/* DARs */



data dac;
  set methyl_access;
  if diff_22C_100U_0U=. or  diff_4C_100U_0U=. then delete;
  if diff_4C_22C < 0 then diff_direction=-1;
  else diff_direction=1;
  if diff_22C_100U_0U < 0 then flag_22C_0U_gt_100U=1; else flag_22C_0U_gt_100U=0;
  if diff_4C_100U_0U < 0 then flag_4C_0U_gt_100U=1; else flag_4C_0U_gt_100U=0;
  keep site_type chr pos diff_direction
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
  prev_direction=lag1(diff_direction);
  if first.chr then region_num=1;
  else do;
     if prev_pos >= pos-100 and prev_direction=diff_direction then region_num=region_num;
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
site_type     flag_num_sites_ge3

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
GC       |1126065 | 444367 |1570432
         |  71.70 |  28.30 | 100.00
         |  71.70 |  28.30 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total     1126065   444367  1570432
            71.70    28.30   100.00

*/


/* Make permanent */

data cold.cytosine_to_access_region_index2;
  set group_dac;
run;

