libname wgbs '/home/jrbnewman/concannon/arabidopsis_wgbs_cold/sas_data';
libname wgbsloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';

ods listing; ods html close;

/* DMC processing:
  (1) Import and fitler DMC results
  (2) FDR corrections across sites
  (3) Compare (by region, site type, dose) with DMR/DAR results
  P<0.05
  FDR P<0.05
*/

proc datasets lib=work kill noprint;
run;
quit;

data good_cov;
  set wgbs.flag_coverage_ge10;
  if site_type ="GC" then do;
     if flag_22C_0U_coverage_ge10=1 and flag_22C_100U_coverage_ge10=1
     and flag_4C_0U_coverage_ge10=1 and flag_4C_100U_coverage_ge10=1
     then output;
     end;
  else do;
  if flag_22C_0U_coverage_ge10=1 and flag_4C_0U_coverage_ge10=1 then output;
     end;
  keep chr pos site_type;
run;


data good_cov_fans_05 good_cov_Fans_15 good_cov_fans_25 good_cov_fans_5;
  set wgbs.flag_coverage_ge10;
  where site_type="GC";
  if flag_22C_100U_coverage_ge10=1 and flag_FANS_0p5U_coverage_ge10=1 then output good_cov_fans_05;
  if flag_22C_100U_coverage_ge10=1 and flag_FANS_1p5U_coverage_ge10=1 then output good_cov_fans_15;
  if flag_22C_100U_coverage_ge10=1 and flag_FANS_5U_coverage_ge10=1 then output good_cov_fans_5;
  if flag_22C_100U_coverage_ge10=1 and flag_FANS_25U_coverage_ge10=1 then output good_cov_fans_25;
  keep chr pos site_type;
run;

data methyl;
set wgbs.flag_methylation_gt0;
    if site_type = "GC" then do;
     if flag_22C_0U_methylation_gt0=1 or flag_22C_100U_methylation_gt0=1
     or flag_4C_0U_methylation_gt0=1 or flag_4C_100U_methylation_gt0=1
     then output;
     end;
  else do;
  if flag_22C_0U_methylation_gt0=1 
     or flag_4C_0U_methylation_gt0=1 then output;
     end;
  keep chr pos site_type;
run;
  

data methyl_fans_05 methyl_fans_15 methyl_fans_5 methyl_fans_25;
set wgbs.flag_methylation_FANS_gt0;
   if flag_22C_100U_methylation_gt0=1 or flag_FANS_0p5U_methylation_gt0=1 then output methyl_fans_05;
   if flag_22C_100U_methylation_gt0=1 or flag_FANS_1p5U_methylation_gt0=1 then output methyl_fans_15;
   if flag_22C_100U_methylation_gt0=1 or flag_FANS_5U_methylation_gt0=1 then output methyl_fans_5;
   if flag_22C_100U_methylation_gt0=1 or flag_FANS_25U_methylation_gt0=1 then output methyl_fans_25;
  keep chr pos site_type;
run;


proc sort data=good_cov nodup;  by chr pos site_type;
proc sort data=good_cov_fans_05; by chr pos site_type;
proc sort data=good_cov_fans_15; by chr pos site_type;
proc sort data=good_cov_fans_5; by chr pos site_type;
proc sort data=good_cov_fans_25; by chr pos site_type;

proc sort data=methyl; by chr pos site_type;
proc sort data=methyl_fans_05; by chr pos site_type;
proc sort data=methyl_fans_15; by chr pos site_type;
proc sort data=methyl_fans_5; by chr pos site_type;
proc sort data=methyl_fans_25; by chr pos site_type;
 run;



 data meth_data_all;
   merge good_cov (in=in1) methyl (in=in2);
   by chr pos site_type;
   if in1 and in2;
   run;

   data meth_data_all_fans_05;
   merge good_cov_fans_05 (in=in1) methyl_fans_05 (in=in2);                                                                                                                       by chr pos site_type;
    if in1 and in2;
       run;

   data meth_data_all_fans_15;
   merge good_cov_fans_15 (in=in1) methyl_fans_15 (in=in2);                                                                                                                       by chr pos site_type;                                                                                                                           if in1 and in2;
          run;
   data meth_data_all_fans_5;
   merge good_cov_fans_5 (in=in1) methyl_fans_5 (in=in2);                                                                                                                       by chr pos site_type;                                                                                                                           if in1 and in2;
          run;
   data meth_data_all_fans_25;
   merge good_cov_fans_25 (in=in1) methyl_fans_25 (in=in2);                                                                                                                       by chr pos site_type;                                                                                                                           if in1 and in2;
          run;


%macro importDMR(inFile,site,group1,group2);

    data WORK.input_data    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!PATCON/arabidopsis_wgbs_cold/analysis_output/dmap2_output/&inFile."
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chr $6. ;
       informat pos best32. ;
       informat _&group1. best32. ;
       informat _&group2.  best32. ;
       informat diff best32. ;
       informat pval $32. ;
       informat call $3. ;
       format chr $6. ;
       format pos best12. ;
       format _&group1. best12. ;
       format _&group2.  best12. ;
       format diff best12. ;
       format pval $32. ;
       format call $3. ;
       input
                   chr $
                   pos
                   _&group1.
                   _&group2.
                   diff
                   pval $
                   call $
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;


/* Subset testable sites */

data input_data1;
  set input_data;
  pos2=pos+1;
  diff2=diff*-1;
  p_value=pval * 1;
  if p_value=. then do;
     if find(pval,'e-3','i') ge 1 then p_value=0;
        else p_value=.;
           end;	drop pos diff;
  rename pos2=pos diff2=diff;
  run;


data meth_data;
  set meth_data_all;
  where site_type="&site.";
  keep chr pos;
run;

proc sort data=input_data1;
  by chr pos;
proc sort data=meth_data nodup;
  by chr pos;
run;

data input_data2;
   merge meth_data (in=in1) input_data1 (in=in2);
   by chr pos;
   if in1 and in2;
   run;


proc multtest inpvalues(p_value)=input_data2 fdr
              out=sig_fdr noprint;
run; quit;

data wgbsloc.dmap2_&site._&group1._&group2._dmc;
   set sig_fdr;
   if p_value=. then flag_p05=.;
   else if p_value < 0.05 then flag_p05=1;
   else flag_p05=0;

   if fdr_p=. then flag_fdr05=.;
   else if fdr_p < 0.05 then flag_fdr05=1;
   else flag_fdr05=0;
run;

proc freq data=wgbsloc.dmap2_&site._&group1._&group2._dmc;
  tables flag_p05 flag_fdr05;
  run;
  
%mend;



%importDMR(CG_22C_0U.vs.4C_0U.full.csv,CG,4C_0U,22C_0U);
%importDMR(CHG_22C_0U.vs.4C_0U.full.csv,CHG,4C_0U,22C_0U);
%importDMR(CHH_22C_0U.vs.4C_0U.full.csv,CHH,4C_0U,22C_0U);

%importDMR(GC_22C_0U.vs.22C_100U.full.csv,GC,22C_100U,22C_0U);
%importDMR(GC_4C_0U.vs.4C_100U.full.csv,GC,4C_100U,4C_0U);


%macro importFANS(inFile,site,FANS,group1,group2);

    data WORK.input_data    ;
        %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!PATCON/arabidopsis_wgbs_cold/analysis_output/dmap2_output/&inFile."
    delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
           informat chr $6. ;
          informat pos best32. ;
         informat _&group1. best32. ;
        informat _&group2.  best32. ;
       informat diff best32. ;
              informat pval $32. ;
             informat call $3. ;
            format chr $6. ;
           format pos best12. ;
          format _&group1. best12. ;
         format _&group2.  best12. ;
        format diff best12. ;
       format pval $32. ;
              format call $3. ;
             input
                        chr $
                   pos
                      _&group1.
                         _&group2.
                    diff
                       pval $
                          call $
         ;
        if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;


       /* Subset testable sites */
       data input_data1;
         set input_data;
   pos2=pos+1;
     diff2=diff*-1;
     p_value=pval*1;
       if p_value=. then do;
           if find(pval,'e-3','i') ge 1 then p_value=0;
	       else p_value=.;
	           end;
       drop pos diff;
         rename pos2=pos diff2=diff;
   run;


   data meth_data;
     set meth_data_all_&FANS.;
       where site_type="&site.";
         keep chr pos;
 run;

 proc sort data=input_data1;
   by chr pos;
   proc sort data=meth_data nodup;
     by chr pos;
     run;

     data input_data2;
        merge meth_data (in=in1) input_data1 (in=in2);
   by chr pos;
      if in1 and in2;
         run;


 proc multtest inpvalues(p_value)=input_data2 fdr
               out=sig_fdr noprint;
       run; quit;

       data wgbsloc.dmap2_&site._&group1._&group2._dmc;
          set sig_fdr;
     if p_value=. then flag_p05=.;
        else if p_value < 0.05 then flag_p05=1;
   else flag_p05=0;

      if fdr_p=. then flag_fdr05=.;
         else if fdr_p < 0.05 then flag_fdr05=1;
    else flag_fdr05=0;
    run;

    proc freq data=wgbsloc.dmap2_&site._&group1._&group2._dmc;
      tables flag_p05 flag_fdr05;
        run;

%mend;



%importFANS(GC_22C_100U.vs.FANS_0p5U.full.csv,GC,fans_05,FANS_0p5U,22C_100U);
%importFANS(GC_22C_100U.vs.FANS_1p5U.full.csv,GC,fans_15,FANS_1p5U,22C_100U);
%importFANS(GC_22C_100U.vs.FANS_5U.full.csv,GC,fans_5,FANS_5U,22C_100U);
%importFANS(GC_22C_100U.vs.FANS_25U.full.csv,GC,fans_25,FANS_25U,22C_100U);




/* For DARs, count intersection between dose and units */


data dmap2_FANS_DMCs;                                    
   set wgbsloc.dmap2_GC_FANS_0p5U_22C_100U_dmc
   wgbsloc.dmap2_GC_FANS_1p5U_22C_100U_dmc
   wgbsloc.dmap2_GC_FANS_5U_22C_100U_dmc
   wgbsloc.dmap2_GC_FANS_25U_22C_100U_dmc;
   rename fdr_p = fdr_p_old flag_fdr05=flag_fdr05_old;
run;


proc multtest inpvalues(p_value)=dmap2_FANS_DMCs fdr
              out=sig_fdr noprint;
run; quit;

data wgbsloc.dmap2_FANS_GC_DMCs;
   set sig_fdr;
   if fdr_p=. then flag_fdr05=.;
   else if fdr_p < 0.05 then flag_fdr05=1;
   else flag_fdr05=0;
run;




