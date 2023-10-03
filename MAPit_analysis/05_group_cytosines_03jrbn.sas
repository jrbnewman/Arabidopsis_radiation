libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* FDR - DMCs and DACs */


data all_dmc;
   set wgbsloca.results_by_dmc_0gy_01gy_: (in=in1)
       wgbsloca.results_by_dmc_0gy_1gy_: (in=in2);
   length comparison $10.;
   if in1 then comparison="0Gy_vs_01Gy";
   if in2 then comparison="0Gy_vs_1Gy";
run;

proc sort data=all_dmc nodup;
   by comparison site_type XP2_FISH chr start_pos stop_pos ;
run;

proc multtest inpvalues(XP2_FISH)=all_dmc fdr noprint out=all_dmc_fdr;
by  comparison site_type;
run;

data all_dmc_fdr2;
  set all_dmc_fdr;
  rename XP2_FISH=FET_P fdr_p=FET_FDR_P;
run;

/* Need to get methylation rates and calculate differences */

data methyl_data;
    set wgbslocA.methylation_data_cg_chg_chh;
    if chr="Mt" or chr="Pt" then delete;
run;

proc sort data=methyl_data;
   by chr start_pos stop_pos site_type treatment units;
proc means data=methyl_data noprint;
   by chr start_pos stop_pos site_type treatment units;
   var total_C perc_methyl;
   output out=methyl_perc sum(total_C)=total_C mean(perc_methyl)=perc_methyl;
run;

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
  if _0Gy0U=. then delete;
  if _01Gy0U=. and _1Gy0U=. then delete;
  if sum(_0Gy0U, _01Gy0U, _1Gy0U) = 0 then delete;
  drop _NAME_;
  rename _01Gy0U=methyl_01Gy
         _1Gy0U=methyl_1Gy
         _0Gy0U=methyl_0Gy;
run;

data methyl_sbys2_01gy;
   set methyl_sbys2;
   length comparison $10.;
   if methyl_0Gy=0 and methyl_01Gy=0 then delete;
   if methyl_0Gy=. or methyl_01Gy=. then delete;
   comparison="0Gy_vs_01Gy";
   methyl_diff=methyl_01Gy-methyl_0Gy;
   drop methyl_1Gy;
   rename methyl_0Gy=methyl_CTL
          methyl_01Gy=methyl_TRT;
 run;

data methyl_sbys2_1gy;
   set methyl_sbys2;
   length comparison $10.;
   if methyl_0Gy=0 and methyl_1Gy=0 then delete;
   if methyl_0Gy=. or methyl_1Gy=. then delete;
   comparison="0Gy_vs_1Gy";
   methyl_diff=methyl_1Gy-methyl_0Gy;
   drop methyl_01Gy;
   rename methyl_0Gy=methyl_CTL
          methyl_1Gy=methyl_TRT;
 run;


data methyl_sbys2_stack;
  set methyl_Sbys2_01Gy methyl_Sbys2_1Gy;
run;



proc sort data=methyl_sbys2_stack;
    by comparison site_type chr start_pos stop_pos;
proc sort data=all_dmc_fdr2;
    by comparison site_type chr start_pos stop_pos;
run;

data methyl_data_all ;
  merge methyl_sbys2_stack (in=in1) all_dmc_fdr2 (in=in2);
  by comparison site_type chr start_pos stop_pos;
  if in1 ;
run;

data methyl_data_all2;
  set methyl_data_all;
  if methyl_diff = . then flag_meth_diff=.;
  else if methyl_diff > 0 then flag_meth_diff=1;
  else if methyl_diff < 0 then flag_meth_diff=-1;
  else flag_meth_diff=0;
  if abs(methyl_diff) >= 0.1 then flag_meth_diff_10perc=1;
  else  flag_meth_diff_10perc=0;
  if abs(methyl_diff) >= 0.2 then flag_meth_diff_20perc=1;
  else  flag_meth_diff_20perc=0;
run;

data wgbslocA.results_by_DMC_w_meth_v4;
    set methyl_Data_all2;
run;




/* Group DMCs into windows
 Sites must be no more than 100bp away from the next site
 Sites must have the same direction to be included
 Region must have at least 3 sites
 A DMR has at least 3 DMCs (flag 1 and 3)
 */

data dmc;
set  wgbslocA.results_by_DMC_w_meth_v4;
where flag_meth_diff_10perc=1;
run;

proc sort data=dmc;
by comparison site_type chr start_pos stop_pos;
run;

data group_dmc;
retain region_num;
set dmc;
by comparison site_type chr;
   prev_pos=lag1(start_pos);
   prev_direction=lag1(flag_meth_diff);
   if first.chr then region_num=1;
   else do;
   if prev_pos >= start_pos-100 and prev_direction=flag_meth_diff then region_num=region_num;
  else region_num=region_num + 1;
 end;
run;


proc freq data=group_dmc noprint;
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

proc freq data=flag_region_ge2 noprint;
tables comparison*site_type*flaG_num_sites_ge2 / out=dmr_check1;
tables comparison*site_type*flaG_num_sites_ge3 / out=dmr_check2;
tables comparison*site_type*flaG_num_sites_ge5 / out=dmr_check3;
tables comparison*site_type*flaG_num_sites_ge10 / out=dmr_check4;
run;

proc print data=dmr_check1;
proc print data=dmr_check2;
proc print data=dmr_check3;
proc print data=dmr_check4;
run;


/* 10%:


Comparison  siteType    2       3       5       10
01Gy_vs_0Gy CG          99645   36710   7009    351
01Gy_vs_0Gy CHG         107861  45783   12952   1152
01Gy_vs_0Gy CHH         439582  224007  82058   19538
1Gy_vs_0Gy  CG          86854   30649   5421    214
1Gy_vs_0Gy  CHG         97412   39254   9365    504
1Gy_vs_0Gy  CHH         421940  215662  78056   14585

20%:

Comparison  siteType    2       3       5       10
01Gy_vs_0Gy CG          13353   3275    416     33
01Gy_vs_0Gy CHG         20812   7027    1108    49
01Gy_vs_0Gy CHH         65920   36219   16486   4390
1Gy_vs_0Gy  CG          11695   2638    233     8
1Gy_vs_0Gy  CHG         20093   6680    957     19
1Gy_vs_0Gy  CHH         71421   37332   15046   3224


 
   
*/

data region2keep;
set flag_region_ge2;
where flaG_num_sites_ge2=1;
keep comparison site_type chr region_num;
run;

proc sort data=region2keep;
by comparison site_type chr region_num;
proc sort data=group_dmc;
by comparison site_type chr region_num;
run;

data group_dmc2;
merge group_dmc (in=in1) region2keep (in=in2);
by comparison site_type chr region_num;
if in1 and in2;
run;


/* Make permanent */

data wgbslocA.cytosine_to_meth_region_index_v4;
set group_dmc;
run;

data wgbslocA.meth_region_flag_site_ge2_v4;
set flag_region_ge2;
run;


data wgbslocA.cytosine_to_meth_region_ge2_v4;
set group_dmc2;
run;



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
   if abs(methyl_diff_TRT_CTL) >= 0.1 then flag_meth_diff_10perc=1;
   else flag_meth_diff_10perc=0;
   if abs(methyl_diff_TRT_CTL) >= 0.2 then flag_meth_diff_20perc=1;
   else flag_meth_diff_20perc=0;
   if FET_FDR_P_TRT_100U_0U=1 or FET_FDR_P_CTL_100U_0U=1 then flag_FDR05_CTL_or_TRT=1;
   else flag_FDR05_CTL_or_TRT=0;
run;


data wgbslocA.results_by_dac_w_meth_v4;
  set dac_meth_all2;
run;


/* Group DMCs into windows
 Sites must be no more than 100bp away from the next site
 Sites must have the same direction to be included
 Region must have at least 3 sites
 A DMR has at least 3 DMCs (flag 1 and 3)
 */

data dac;
set dac_meth_all2;
if flag_meth_diff_10perc = 1 and flag_FDR05_CTL_or_TRT = 1;
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

data wgbslocA.cytosine_to_acc_region_index_v4;
set group_dac2;
run;

data wgbslocA.acc_region_flag_site_ge2_v4;
set flag_region_ge2;
run;


data wgbslocA.cytosine_to_acc_region_ge2_v4;
set group_dac3;
run;


/* Export BED files for various tools */


data dmrs;
   set wgbslocA.cytosine_to_meth_region_ge2_v4;
   length regionID $100.;
   regionID = catx("_",site_type,comparison,chr,region_num);
   keep site_type comparison chr start_pos stop_pos regionID;
run;

proc sort data=dmrs;
  by site_type comparison chr regionID;
proc means data=dmrs noprint;
   by site_type comparison chr regionID;
   var start_pos stop_pos;
   output out=dmr_regions min(start_pos)=start_pos max(stop_pos)=stop_pos;
run;

proc sort data=dmr_regions;
  by chr start_pos stop_pos comparison regionID;
run;



data dmr_cg_01 dmr_cg_1 
     dmr_chg_01 dmr_chg_1
     dmr_chh_01 dmr_chh_1;
     retain chr start_pos stop_pos regionID score strand;
     set dmr_regions;
     score=0;
     strand="+";
     if site_type = "CG" and comparison = "0Gy_vs_01G" then output dmr_cg_01;
     if site_type = "CG" and comparison = "0Gy_vs_1Gy" then output dmr_cg_1;
     if site_type = "CHG" and comparison = "0Gy_vs_01G" then output dmr_chg_01;
     if site_type = "CHG" and comparison = "0Gy_vs_1Gy" then output dmr_chg_1;
     if site_type = "CHH" and comparison = "0Gy_vs_01G" then output dmr_chh_01;
     if site_type = "CHH" and comparison = "0Gy_vs_1Gy" then output dmr_chh_1;
     keep chr start_pos stop_pos regionID score strand;
run;






data dars;
   set wgbslocA.cytosine_to_acc_region_ge2_v4;
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


proc export data=dmr_cg_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_10cGy.bed" dbms=tab replace; putnames=no; run;
proc export data=dmr_cg_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_100cGy.bed" dbms=tab replace; putnames=no; run;

proc export data=dmr_chg_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_10cGy.bed" dbms=tab replace; putnames=no; run;
proc export data=dmr_chg_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_100cGy.bed" dbms=tab replace; putnames=no; run;

proc export data=dmr_chh_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHH_10cGy.bed" dbms=tab replace; putnames=no; run;
proc export data=dmr_chh_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHH_100cGy.bed" dbms=tab replace; putnames=no; run;

proc export data=dar_gc_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_10cGy.bed" dbms=tab replace; putnames=no; run;
proc export data=dar_gc_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_100cGy.bed" dbms=tab replace; putnames=no; run;




