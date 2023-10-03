ods html close;
ods listing;
libname wgbs "!PATCON/arabidopsis_wgbs_cold/sas_data";
libname wgbsloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";

/* Merge in all possible testable methylation sites and flag is analyzable or not */

proc datasets lib=work kill noprint;
run;
quit;



%macro flagDMCs(site);

/* Get list of all sites of type */
data all_sites;
  set wgbs.flag_methylation_gt0;
  where site_type="&site.";
      if flag_4C_0U_methylation_gt0=1 
	  or flag_22C_0U_methylation_gt0=1 
	  then output;  
	 keep chr pos;
	run;

/* Get average methylayion per group */


data contrast;
   set wgbsloc.dmap2_&site._4C_0U_22C_0U_dmc;
   rename _4C_0U=meth_4C_0U _22C_0U=meth_22C_0U diff=diff_4C_22C;
run;

proc sort data=all_sites;
  by chr pos;
proc sort data=contrast;
  by chr pos;
  run;

data meth_sites2;
  merge all_sites (in=in1) contrast (in=in2);
  by chr pos;
  if in1 and in2;
  run;

data pos2region;
   set wgbsloc.meth_pos2region;
   length gene_id $16.;   
   length region_type $16.;
   do i=1 by 1 while(scan(gene_id_cat,i,"|") ^= "");
       gene_id=scan(gene_id_cat,i,"|");
		region_type=scan(region_type_cat,i,"|");
		output;
		end;
keep region_type gene_id chr pos;
run;

proc sort data=meth_sites2;
  by chr pos;
proc sort data=pos2region nodup;
  by chr pos gene_id region_type;
  run;

 data meth_sites2region;
   merge meth_sites2 (in=in1) pos2region (in=in2);
   by chr pos;
   if in2 then flag_orphan_site=0;
   else flag_orphan_site=1;
   if in1;
   run;

  data meth_sites2region2;
    set meth_sites2region;
	if abs(diff_4C_22C) >= 0.1 then flag_diff10=1; else flag_diff10=0;
	if abs(diff_4C_22C) >= 0.2 then flag_diff20=1; else flag_diff20=0;
	if abs(diff_4C_22C) >= 0.3 then flag_diff30=1; else flag_diff30=0;
	if abs(diff_4C_22C) >= 0.5 then flag_diff50=1; else flag_diff50=0;
	run;

	data wgbsloc.results_by_dmc_&site._4C_22C;
	   set meth_sites2region2;
	   run;

	%mend;


%flagDMCs(CG);
%flagDMCs(CHG);
%flagDMCs(CHH);


/* Flag DACs -- cold tolerance */

/* Get list of all sites of type */
data all_sites;
  set wgbs.flag_methylation_gt0;
  where site_type="GC";
      if flag_4C_0U_methylation_gt0=1 or flag_4C_100U_methylation_gt0=1 
	  or flag_22C_0U_methylation_gt0=1 or flag_22C_100U_methylation_gt0=1 
	  then output;
	 keep chr pos;
	run;


data contrast1;
  set wgbsloc.dmap2_gc_4c_100u_4c_0u_dmc;
  length temperature $3.;
  temperature="4C";
  rename _4C_100U=meth_100U _4C_0U=meth_0U fdr_p=fdr_p_old
   flag_fdr05=flag_fdr05_old;
run;

data contrast2;
  set wgbsloc.dmap2_gc_22c_100u_22c_0u_dmc;
  length temperature $3.;
  temperature="22C";
  rename _22C_100U=meth_100U _22C_0U=meth_0U fdr_p=fdr_p_old
   flag_fdr05=flag_fdr05_old;
run;

data contrasts;
  set contrast1 contrast2;
run;


proc multtest inpvalues(p_value)=contrasts fdr noprint out=contrasts_fdr;
run;
quit;


data contrast_22C;
  set contrasts_fdr;
  where temperature = "22C";
  methdiff_22C_100U_0U=diff;
  if fdr_P = . then flag_FDR05_22C_100U_0U=.;
  else if fdr_P < 0.05 then flag_FDR05_22C_100U_0U=1;
  else flag_FDR05_22C_100U_0U=0;
  keep chr pos meth_0U meth_100U p_value fdr_p flag_p05 flag_FDR05_22C_100U_0U methdiff_22C_100U_0U;
  rename p_value=P_22C_100U_0U
  fdr_p=fdr_22C_100U_0U
  flag_p05=flag_P05_22C_100U_0U
  meth_0U=meth_22C_0U meth_100U=meth_22C_100U;
run;

data contrast_4C;
  set contrasts_fdr;
  where temperature = "4C";
  methdiff_4C_100U_0U=diff;
  if fdr_P = . then flag_FDR05_4C_100U_0U=.;
  else if fdr_P < 0.05 then flag_FDR05_4C_100U_0U=1;
  else flag_FDR05_4C_100U_0U=0;
  keep chr pos meth_0U meth_100U p_value fdr_p flag_p05 flag_FDR05_4C_100U_0U methdiff_4C_100U_0U;
  rename p_value=P_4C_100U_0U
  fdr_p=fdr_4C_100U_0U
  flag_p05=flag_P05_4C_100U_0U
  meth_0U=meth_4C_0U meth_100U=meth_4C_100U;
run;




proc sort data=all_sites;
  by chr pos;
proc sort data=contrast_22C;
  by chr pos;
proc sort data=contrast_4C;
  by chr pos;
  run;

data meth_sites2;
  merge all_sites (in=in1) contrast_22C (in=in2) contrast_4C (in=in3);
  by chr pos;
  if in1;
  run;

data meth_sites3;
  set meth_sites2;
  if meth_4C_0U > meth_4C_100U then flag_4C_0U_gt_100U=1; else flag_4C_0U_gt_100U=0;
  if meth_22C_0U > meth_22C_100U then flag_22C_0U_gt_100U=1; else flag_22C_0U_gt_100U=0;
run;


data pos2region;
   set wgbsloc.meth_pos2region;
   length gene_id $16.;   
   length region_type $16.;
   do i=1 by 1 while(scan(gene_id_cat,i,"|") ^= "");
       gene_id=scan(gene_id_cat,i,"|");
		region_type=scan(region_type_cat,i,"|");
		output;
		end;
keep region_type gene_id chr pos;
run;

proc sort data=meth_sites3;
  by chr pos;
proc sort data=pos2region nodup;
  by chr pos gene_id region_type;
  run;

 data meth_sites2region;
   merge meth_sites3 (in=in1) pos2region (in=in2);
   by chr pos;
   if in2 then flag_orphan_site=0;
   else flag_orphan_site=1;
   if in1;
   run;

  data meth_sites2region2;
    set meth_sites2region;
	    diff_4C_100U_0U=meth_4C_100U-meth_4C_0U;
	    diff_22C_100U_0U=meth_22C_100U-meth_22C_0U;
		diff_4C_22C=diff_4C_100U_0U-diff_22C_100U_0U;
        if meth_4C_100U < meth_4C_0U then flag_4C_100U_lt_0U=1;
        else flag_4C_100U_lt_0U=0;
        if meth_22C_100U < meth_22C_0U then flag_22C_100U_lt_0U=1;
        else flag_22C_100U_lt_0U=0;

    if abs(diff_22C_100U_0U) >= 0.1 then flag_22C_diff10=1; else flag_22C_diff10=0;
    if abs(diff_22C_100U_0U) >= 0.2 then flag_22C_diff20=1; else flag_22C_diff20=0;
    if abs(diff_22C_100U_0U) >= 0.3 then flag_22C_diff30=1; else flag_22C_diff30=0;
    if abs(diff_22C_100U_0U) >= 0.5 then flag_22C_diff50=1; else flag_22C_diff50=0;

    if abs(diff_4C_100U_0U) >= 0.1 then flag_4C_diff10=1; else flag_4C_diff10=0;
    if abs(diff_4C_100U_0U) >= 0.2 then flag_4C_diff20=1; else flag_4C_diff20=0;
    if abs(diff_4C_100U_0U) >= 0.3 then flag_4C_diff30=1; else flag_4C_diff30=0;
    if abs(diff_4C_100U_0U) >= 0.5 then flag_4C_diff50=1; else flag_4C_diff50=0;

	if abs(diff_4C_22C) >= 0.1 then flag_diff10=1; else flag_diff10=0;
	if abs(diff_4C_22C) >= 0.2 then flag_diff20=1; else flag_diff20=0;
	if abs(diff_4C_22C) >= 0.3 then flag_diff30=1; else flag_diff30=0;
	if abs(diff_4C_22C) >= 0.5 then flag_diff50=1; else flag_diff50=0;
	run;

	data wgbsloc.results_by_dmc_GC_4c_22c;
	   set meth_sites2region2;
	   run;



/* Flag DACs -- 22C 100U vs FANS 

   Going to stack everything together for ease */

/* Get list of all sites of type */
data all_sites;
  set wgbs.flag_methylation_fans_gt0;
      if flag_22C_100U_methylation_gt0=1 or flag_FANS_0p5U_methylation_gt0=1
      or flag_FANS_1p5U_methylation_gt0=1 or flag_FANS_5U_methylation_gt0=1
      or flag_FANS_25U_methylation_gt0=1
      then output;
	  keep chr pos;
   run;


data contrasts;
  set wgbsloc.dmap2_fans_gc_dmcs;
  length FANS_units $4.;
  meth_22C_100U=_22C_100U;
  if _FANS_0p5U ne . then do; 
     FANS_units="0.5U";
     meth_FANS=_FANS_0p5U;
     end;
  else if _FANS_1p5U ne . then do; 
     FANS_units="1.5U";
     meth_FANS=_FANS_1p5U;
     end;
  else if _FANS_5U ne . then do; 
     FANS_units="5U";
     meth_FANS=_FANS_5U;
     end;
  else if _FANS_25U ne . then do; 
     FANS_units="25U";
     meth_FANS=_FANS_25U;
     end;
   else delete;
  diff_FANS_22C_100U=diff;
   drop _22C_100U _FANS_0p5U _FANS_1p5U _FANS_5U _FANS_25U pval call diff fdr_p_old flag_fdr05_old ;
run;




proc sort data=all_sites;
  by chr pos;
proc sort data=contrasts;
  by chr pos;
  run;

data meth_sites2;
  merge all_sites (in=in1) contrasts (in=in2) ;
  by chr pos;
  if in1;
  run;

data pos2region;
   set wgbsloc.meth_pos2region;
keep gene_id_cat region_type_cat chr pos;
run;

proc sort data=meth_sites2;
  by chr pos;
proc sort data=pos2region nodup;
  by chr pos ;
  run;

 data meth_sites2region;
   merge meth_sites2 (in=in1) pos2region (in=in2);
   by chr pos;
   if in2 then flag_orphan_site=0;
   else flag_orphan_site=1;
   if in1;
   run;


data meth_sites2region2;
   set meth_sites2region;
   length gene_id $16.;   
   length region_type $16.;
   do i=1 by 1 while(scan(gene_id_cat,i,"|") ^= "");
       gene_id=scan(gene_id_cat,i,"|");
		region_type=scan(region_type_cat,i,"|");
		output;
		end;
drop gene_id_cat region_type_cat i;
run;

  data meth_sites2region3;
    set meth_sites2region2;
	if abs(diff_FANS_22C_100U) >= 0.1 then flag_diff10=1; else flag_diff10=0;
	if abs(diff_FANS_22C_100U) >= 0.2 then flag_diff20=1; else flag_diff20=0;
	if abs(diff_FANS_22C_100U) >= 0.3 then flag_diff30=1; else flag_diff30=0;
	if abs(diff_FANS_22C_100U) >= 0.5 then flag_diff50=1; else flag_diff50=0;
	run;

	data wgbsloc.results_by_dmc_GC_FANS_22c;
	   set meth_sites2region2;
	   run;


/* Correlation between FANS and 22C */

ods graphics / ANTIALIASMAX=4161600;
proc sort data=meth_sites2;
  by FANS_units chr pos;
proc sgplot data=meth_sites2;
   by FANS_units;
  scatter x=meth_22C_100U y=meth_FANS;
run;

proc corr data=meth_sites2 pearson;
  by FANS_units;
  var meth_22C_100U meth_FANS;
run;

/* Only GC sites that are methylated in at least one condition (FANS unit, or 22C 100U)
   with at least 10 mapped reads:

-------------------------------------- FANS_units=0.5U -----------------------------------------

                                      The CORR Procedure

                          2  Variables:    meth_22C_100U meth_FANS


                                      Simple Statistics

ariable                N          Mean       Std Dev           Sum       Minimum       Maximum

eth_22C_100U     4154944       0.47392       0.17847       1969095             0       1.00000
eth_FANS         4154944       0.21582       0.19231        896734             0       1.00000


                        Pearson Correlation Coefficients, N = 4154944
                                  Prob > |r| under H0: Rho=0

                                                meth_
                                             22C_100U      meth_FANS

                          meth_22C_100U       1.00000        0.62091
                                                              <.0001

                          meth_FANS           0.62091        1.00000
                                               <.0001


--------------------------------------- FANS_units=1.5U ------------------------------------------

                                       The CORR Procedure

                           2  Variables:    meth_22C_100U meth_FANS


                                       Simple Statistics

Variable                N          Mean       Std Dev           Sum       Minimum       Maximum

meth_22C_100U     4161572       0.47342       0.17897       1970161             0       1.00000
meth_FANS         4161572       0.39375       0.17590       1638632             0       1.00000


                         Pearson Correlation Coefficients, N = 4161572
                                   Prob > |r| under H0: Rho=0

                                                 meth_
                                              22C_100U      meth_FANS

                           meth_22C_100U       1.00000        0.59664
                                                               <.0001

                           meth_FANS           0.59664        1.00000
                                                <.0001

---------------------------------------- FANS_units=5U -------------------------------------------

                                       The CORR Procedure

                           2  Variables:    meth_22C_100U meth_FANS


                                       Simple Statistics

Variable                N          Mean       Std Dev           Sum       Minimum       Maximum

meth_22C_100U     4154826       0.47335       0.17893       1966677             0       1.00000
meth_FANS         4154826       0.50623       0.17430       2103290             0       1.00000


                         Pearson Correlation Coefficients, N = 4154826
                                   Prob > |r| under H0: Rho=0

                                                 meth_
                                              22C_100U      meth_FANS

                           meth_22C_100U       1.00000        0.56563
                                                               <.0001

                           meth_FANS           0.56563        1.00000
                                                <.0001

----------------------------------------- FANS_units=25U -----------------------------------------

                                        The CORR Procedure

                            2  Variables:    meth_22C_100U meth_FANS


                                        Simple Statistics

 Variable                N          Mean       Std Dev           Sum       Minimum       Maximum

 meth_22C_100U     4158287       0.47328       0.17892       1968035             0       1.00000
 meth_FANS         4158287       0.57658       0.14300       2397586             0       1.00000


                          Pearson Correlation Coefficients, N = 4158287
                                    Prob > |r| under H0: Rho=0

                                                  meth_
                                               22C_100U      meth_FANS

                            meth_22C_100U       1.00000        0.48047
                                                                <.0001

                            meth_FANS           0.48047        1.00000
                                                 <.0001

*/


