libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


/* DARs

define a DAR as: at least one site with diff > 0.2 and FDR < 0.05, with same differential direction, < 100bp in between

 */

proc datasets lib=work kill noprint;
run;
quit;


data all_dac;
   set wgbsloca.results_dac_0gy_100U_0gy_0U_: (in=in1)
       wgbsloca.results_dac_01gy_100U_01gy_0U_: (in=in2)
       wgbsloca.results_dac_1gy_100U_1gy_0U_: (in=in3);
   length comparison $15.;
   if in1 then comparison="0Gy_100U_0U";
   if in2 then comparison="01Gy_100U_0U";
   if in3 then comparison="1Gy_100U_0U";
run;

proc sort data=all_dac nodup;
   by comparison site_type  chr start_pos stop_pos ;
run;


/* Need to get methylation rates and calculate differences */

data methyl_data;
    set wgbslocA.methylation_data_gc;
    if flag_normalized=1;
    if chr="Mt" or chr="Pt" then delete;
run;

proc sort data=methyl_data;
   by chr start_pos stop_pos site_type treatment units;
proc means data=methyl_data noprint;
   by chr start_pos stop_pos site_type treatment units;
   var total_C perc_methyl_norm;
   output out=methyl_perc sum(total_C)=total_C mean(perc_methyl_norm)=perc_methyl;
run;

/* remove sites without at least 10 reasd in wither 100U or 0U
    then BEFORE FDR, I want to remove sites that I can't use for pairwise comparisons
    (e.g. for 0Gy vs 0.1Gy, only if in both, 1Gy status irrelevant)

    So, the ones I'd KEEP for 0.1Gy have to appear in 0.1Gy and 0Gy
                          for 1Gy have to appear in 1Gy and 0Gy
                          for 0Gy has to appear in EITHER [0.1Gy and 0Gy] OR [1Gy and 0Gy]
*/

data methyl_perc2;
  set methyl_perc;
  if total_C < 10 then delete;
run;

proc sort data=methyl_perc2;
   by site_type chr start_pos stop_pos treatment units;
proc transpose data=methyl_perc2 out=methyl_sbys;
  by site_type chr start_pos stop_pos;
  id treatment units;
  var perc_methyl;
run;

data methyl_sbys2;
  set methyl_sbys;
  if _0Gy0U = . or _0Gy100U = . then flag_missing_0Gy=1; else flag_missing_0Gy=0;
  if _01Gy0U = . or _01Gy100U = . then flag_missing_01Gy=1; else flag_missing_01Gy=0;
  if _1Gy0U = . or _1Gy100U = . then flag_missing_1Gy=1; else flag_missing_1Gy=0;

  if _0Gy0U = 0 and _0Gy100U = 0 then flag_nometh_0Gy=1; else flag_nometh_0Gy=0;
  if _01Gy0U = 0 and _01Gy100U = 0 then flag_nometh_01Gy=1; else flag_nometh_01Gy=0;
  if _1Gy0U = 0 and _1Gy100U = 0 then flag_nometh_1Gy=1; else flag_nometh_1Gy=0;
  

  drop _NAME_;
  rename _01Gy0U=methyl_01Gy_0U
         _1Gy0U=methyl_1Gy_0U
         _0Gy0U=methyl_0Gy_0U
         _01Gy100U=methyl_01Gy_100U
         _1Gy100U=methyl_1Gy_100U
         _0Gy100U=methyl_0Gy_100U;
run;

data methyl_sbys2_0Gy;
   set methyl_sbys2;
   length comparison $15.;
   comparison="0Gy_100U_0U";
   if flag_missing_0Gy=0 and (flag_missing_01Gy=0 or flag_missing_1gy=0) and flag_nometh_0Gy=0;
   keep site_type chr start_pos stop_pos comparison;
run;

data methyl_sbys2_01Gy;
   set methyl_sbys2;
   length comparison $15.;
   comparison="01Gy_100U_0U";
   if flag_missing_0Gy=0 and flag_missing_01Gy=0 and flag_nometh_01Gy=0;
   keep site_type chr start_pos stop_pos comparison;
run;

data methyl_sbys2_1Gy;
   set methyl_sbys2;
   length comparison $15.;
   comparison="1Gy_100U_0U";
   if flag_missing_0Gy=0 and flag_missing_1Gy=0 and flag_nometh_1Gy=0;
   keep site_type chr start_pos stop_pos comparison;
run;

data methyl2keep;
   set methyl_sbys2_0Gy methyl_sbys2_01Gy methyl_sbys2_1Gy;
run;

proc sort data=methyl2keep;
   by comparison site_type chr start_pos stop_pos;
proc sort data=all_dac;
   by comparison site_type chr start_pos stop_pos;
run;

data all_dac2;
    merge all_dac (in=in1) methyl2keep (in=in2);
    by comparison site_type chr start_pos stop_pos;
    if in1 and in2;
run;


proc multtest inpvalues(XP2_FISH)=all_dac2 fdr noprint out=all_dac_fdr;
by  comparison site_type;
run;

data all_dac_fdr2;
  set all_dac_fdr;
  rename XP2_FISH=FET_P fdr_p=FET_FDR_P;
run;

proc sort data=all_dac_fdr2;
   by site_type chr start_pos stop_pos comparison;
run;

proc transpose data=all_dac_fdr2 out=all_dac_p_sbys;
   by site_type chr start_pos stop_pos;
   id comparison;
   var FET_P;
run;

proc transpose data=all_dac_fdr2 out=all_dac_fdr_sbys;
   by site_type chr start_pos stop_pos;
   id comparison;
   var FET_FDR_P;
run;

data all_dac_p_sbys2;
   set all_dac_p_sbys;
   rename _01Gy_100U_0U=FET_P_01Gy_100U_0U
          _0Gy_100U_0U=FET_P_0Gy_100U_0U
          _1Gy_100U_0U=FET_P_1Gy_100U_0U;
   drop _NAME_ _LABEL_;
run;

data all_dac_fdr_sbys2;
   set all_dac_fdr_sbys;
   rename _01Gy_100U_0U=FET_FDR_P_01Gy_100U_0U
          _0Gy_100U_0U=FET_FDR_P_0Gy_100U_0U
          _1Gy_100U_0U=FET_FDR_P_1Gy_100U_0U;
   drop _NAME_ _LABEL_;
run;

proc sort data=all_dac_p_sbys2;
  by site_type chr start_pos stop_pos;
proc sort data=all_dac_fdr_sbys2;
  by site_type chr start_pos stop_pos;
run;


data all_dac_sbys3;
   merge all_dac_p_sbys2 all_dac_fdr_sbys2;
  by site_type chr start_pos stop_pos;
run;


/* Need to make a big table of the following:
   site type, chr, start, stop, methyl_TRT_0U methyl_TRT_100U
   methyl_CTL_0U methyl_CTL_100U comparison P/FDR for CTL TRT
   methyl_diff_TRT methyl_diff_CTL methyl_diff_TRT_CTL
   flag_meth_diff flag_meth_diff_20perc */

data methyl_01Gy_0Gy;
   set methyl_sbys2;
   keep site_type chr start_pos stop_pos
        methyl_01gy_0u methyl_01gy_100u
        methyl_0gy_0u methyl_0gy_100u;
   rename methyl_01gy_0u=methyl_TRT_0U
          methyl_01gy_100u=methyl_TRT_100U
          methyl_0gy_0u=methyl_CTL_0U
          methyl_0gy_100u=methyl_CTL_100U;
run;

data methyl_1Gy_0Gy;
   set methyl_sbys2;
   keep site_type chr start_pos stop_pos
        methyl_1gy_0u methyl_1gy_100u
        methyl_0gy_0u methyl_0gy_100u;
   rename methyl_1gy_0u=methyl_TRT_0U
          methyl_1gy_100u=methyl_TRT_100U
          methyl_0gy_0u=methyl_CTL_0U
          methyl_0gy_100u=methyl_CTL_100U;

run;


data dac_01Gy_0gy;
   set all_dac_sbys3;
   if FET_P_01Gy_100U_0U=. or FET_P_0Gy_100U_0U=. then delete;
   keep site_type chr start_pos stop_pos FET_P_01Gy_100U_0U FET_P_0Gy_100U_0U
   FET_FDR_P_01Gy_100U_0U FET_FDR_P_0Gy_100U_0U;
   rename FET_P_01Gy_100U_0U=FET_P_TRT_100U_0U
          FET_P_0Gy_100U_0U =FET_P_CTL_100U_0U
          FET_FDR_P_01Gy_100U_0U=FET_FDR_P_TRT_100U_0U
          FET_FDR_P_0Gy_100U_0U=FET_FDR_P_CTL_100U_0U;
run;


data dac_1Gy_0gy;
   set all_dac_sbys3;
   if FET_P_1Gy_100U_0U=. or FET_P_0Gy_100U_0U=. then delete;
   keep site_type chr start_pos stop_pos FET_P_1Gy_100U_0U FET_P_0Gy_100U_0U
   FET_FDR_P_1Gy_100U_0U FET_FDR_P_0Gy_100U_0U;
   rename FET_P_1Gy_100U_0U=FET_P_TRT_100U_0U
          FET_P_0Gy_100U_0U =FET_P_CTL_100U_0U
          FET_FDR_P_1Gy_100U_0U=FET_FDR_P_TRT_100U_0U
          FET_FDR_P_0Gy_100U_0U=FET_FDR_P_CTL_100U_0U;
run;

proc sort data=dac_01gy_0gy;
   by site_type chr start_pos stop_pos;
proc sort data=dac_1gy_0gy;
   by site_type chr start_pos stop_pos;
proc sort data=methyl_01Gy_0Gy;
   by site_type chr start_pos stop_pos;
proc sort data=methyl_1Gy_0Gy;
   by site_type chr start_pos stop_pos;
run;

data dac_meth_01gy_0gy;
  merge dac_01gy_0gy (in=in1) methyl_01Gy_0Gy (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;


data dac_meth_1gy_0gy;
  merge dac_1gy_0gy (in=in1) methyl_1Gy_0Gy (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

data dac_meth_all;
   set dac_meth_01Gy_0gy (in=in1)
       dac_meth_1Gy_0gy (in=in2);
   length comparison $15.;
   if in1 then comparison="01Gy_0Gy";
   if in2 then comparison="1Gy_0Gy";
run;

data dac_meth_all2;
   set dac_meth_all;
   methyl_diff_TRT=methyl_TRT_100U-methyl_TRT_0U;
   methyl_diff_CTL=methyl_CTL_100U-methyl_CTL_0U;
   methyl_diff_TRT_CTL=methyl_diff_TRT-methyl_diff_CTL;
   if methyl_diff_TRT_CTL = .  then flag_meth_diff=.;
   else if methyl_diff_TRT_CTL < 0  then flag_meth_diff=-1;
   else if methyl_diff_TRT_CTL > 0  then flag_meth_diff=1;
   else  flag_meth_diff=0;

   if abs(methyl_diff_TRT_CTL) >= 0.2 then flag_meth_diff_10perc=1;
   else flag_meth_diff_10perc=0;
   if FET_FDR_P_TRT_100U_0U=1 or FET_FDR_P_CTL_100U_0U=1 then flag_FDR05_CTL_or_TRT=1;
   else flag_FDR05_CTL_or_TRT=0;
run;




/* Group DMCs into windows
 Sites must be no more than 100bp away from the next site
 Sites must have the same direction to be included
 Region must have at least 3 sites
 A DMR has at least 3 DMCs (flag 1 and 3)
 */

data dac;
set dac_meth_all2;
if flag_meth_diff_10perc ne 1 then delete;
if flag_FDR05_CTL_or_TRT ne 1 then delete;
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
   if prev_pos >= start_pos-200 and prev_direction=flag_meth_diff then region_num=region_num;
  else region_num=region_num + 1;
 end;
run;


proc freq data=group_dac noprint;
tables comparison*site_type*chr*region_num / out=sites_per_region;
run;

data flag_region_ge2;
set sites_per_region;
if count >= 3 then flag_num_sites_ge3=1;
else flag_num_sites_ge3=0;
run;

proc freq data=flag_region_ge2;
tables comparison*site_type*flaG_num_sites_ge3 / out=dmr_check;
proc print data=dmr_check;
run;

/*
                        site_    flag_num_
   Obs    comparison    type     sites_ge2    COUNT    PERCENT

    1      01Gy_0Gy      GC          0        89128    53.2241
    2      01Gy_0Gy      GC          1        23301    13.9145
    3      1Gy_0Gy       GC          0        49508    29.5644
    4      1Gy_0Gy       GC          1         5521     3.2969


*/


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
where flaG_num_sites_ge3=1;
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

*data wgbslocA.cytosine_to_acc_region_index3;
*set group_dac2;
*run;

*data wgbslocA.acc_region_flag_site_ge3;
*set flag_region_ge2;
*run;


*data wgbslocA.cytosine_to_acc_region_ge3;
*set group_dac3;
*run;

