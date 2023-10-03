libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* iterdataset macro */

%macro iterdataset(dataset=,function=);
    %local dsid now total rows cols rc;
    %let dsid = %sysfunc(open(&dataset));
    %let now = 0;
    %let rows = %sysfunc(attrn(&dsid, nobs));
    %let cols = %sysfunc(attrn(&dsid, nvars));

    %do %while(%sysfunc(fetch(&dsid)) = 0); %* outer loop across rows;
        %let now = %eval(&now + 1);

        %do i = 1 %to &cols; %* inner loop across coloumns;
            %local v t;
            %let v=%sysfunc(varname(&dsid,&i));
            %local &v;
            %let t = %sysfunc(vartype(&dsid,&i));
            %let &v = %sysfunc(getvar&t(&dsid,&i));
        %end;

        %unquote(&function);

    %end;
    %let rc = %sysfunc(close(&dsid));
%mend;

/* Import raw BED files */

data design;
  set wgbsA.design_file;
run;

data design1;
  set design;
  where units="0U";
run;



%macro importBED(trt,unit,rep);

%macro siteType(siteType);


     data WORK.&siteType._&trt._&unit._&rep.    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 "!HOME/concannon/DTRA/radiation_bed_files/arabidopsis/by_rep/at_rad_&siteType._all_&trt._&unit._&rep..bed"
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat chr $2. ;
        informat start_pos best32. ;
        informat stop_pos best32. ;
        informat perc_methyl best32. ;
        informat total_C best32. ;
        informat methyl_C best32. ;
        informat VAR7 $1. ;
        informat VAR8 best32. ;
        informat VAR9 $1. ;
        informat VAR10 best32. ;
        informat VAR11 best32. ;
        informat VAR12 $1. ;
        informat VAR13 best32. ;
        informat VAR14 best32. ;
        format chr $2. ;
        format start_pos best12. ;
        format stop_pos best12. ;
        format perc_methyl best12. ;
        format total_C best12. ;
        format methyl_C best12. ;
        format VAR7 $1. ;
        format VAR8 best12. ;
        format VAR9 $1. ;
        format VAR10 best12. ;
        format VAR11 best12. ;
        format VAR12 $1. ;
        format VAR13 best12. ;
        format VAR14 best12. ;
     input
               chr $
               start_pos
               stop_pos
               perc_methyl
               total_C
               methyl_C
               VAR7 $
               VAR8
               VAR9 $
               VAR10
               VAR11
               VAR12 $
               VAR13
               VAR14
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;



data &siteType._&trt._&unit._&rep._2;
  set &siteType._&trt._&unit._&rep.;
  length treatment $10.;
  length units $10.;
  format rep best12.;
  length site_type $3.;
  treatment="&trt.";
  units="&unit.";
  rep=&rep.;
  site_type=upcase("&siteType.");
  drop VAR7-VAR14;
run;


%mend;

%siteType(cg);
%siteType(chg);
%siteType(chh);
%siteType(gc);


%mend;

%iterdataset(dataset=design, function=%nrstr(%importBED(&treatment,&units,&rep)));



/* Import raw and normalized GC 100U methylation:
   I am going to only use the normalized methylation rate here (but keep the total and methyl counts)
*/



%macro importBED(trt,unit,rep);

     data WORK.GC_&trt._&unit._&rep.    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 "!HOME/concannon/DTRA/radiation_bed_files/arabidopsis/by_rep/at_rad_gc_all_&trt._&unit._&rep..bed"
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat chr $2. ;
        informat start_pos best32. ;
        informat stop_pos best32. ;
        informat perc_methyl best32. ;
        informat total_C best32. ;
        informat methyl_C best32. ;
        informat VAR7 $1. ;
        informat VAR8 best32. ;
        informat VAR9 $1. ;
        informat VAR10 best32. ;
        informat VAR11 best32. ;
        informat VAR12 $1. ;
        informat VAR13 best32. ;
        informat VAR14 best32. ;
        format chr $2. ;
        format start_pos best12. ;
        format stop_pos best12. ;
        format perc_methyl best12. ;
        format total_C best12. ;
        format methyl_C best12. ;
        format VAR7 $1. ;
        format VAR8 best12. ;
        format VAR9 $1. ;
        format VAR10 best12. ;
        format VAR11 best12. ;
        format VAR12 $1. ;
        format VAR13 best12. ;
        format VAR14 best12. ;
     input
               chr $
               start_pos
               stop_pos
               perc_methyl
               total_C
               methyl_C
               VAR7 $
               VAR8
               VAR9 $
               VAR10
               VAR11
               VAR12 $
               VAR13
               VAR14
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;



data gc_&trt._&unit._&rep._2;
  set gc_&trt._&unit._&rep.;
  length treatment $10.;
  length units $10.;
  format rep best12.;
  length site_type $3.;
  treatment="&trt.";
  units="&unit.";
  rep=&rep.;
  site_type=upcase("gc");
  drop VAR7-VAR14;
run;


     data WORK.GCN_&trt._&unit._&rep.    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile

 "!HOME/concannon/DTRA/radiation_bed_files/plantNormalization/at_gch_&trt._&unit._rerun_&rep._NormalizedByGroup.BED"
 delimiter=' ' MISSOVER DSD lrecl=32767 ;
        informat chr $2. ;
        informat start_pos best32. ;
        informat stop_pos best32. ;
        informat total_C_norm best32. ;
        informat methyl_C_norm best32. ;
        informat perc_methyl_norm best32. ;
        format chr $2. ;
        format start_pos best12. ;
        format stop_pos best12. ;
        format total_C_norm best12. ;
        format methyl_C_norm best12. ;
        format perc_methyl_norm best12. ;
     input
               chr $
               start_pos
               stop_pos
               total_C_norm
               methyl_C_norm
               perc_methyl_norm

   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

data gcN_&trt._&unit._&rep._2;
  set gcN_&trt._&unit._&rep.;
  length treatment $10.;
  length units $10.;
  format rep best12.;
  length site_type $3.;
  treatment="&trt.";
  units="&unit.";
  rep=&rep.;
  site_type=upcase("gc");
run;

%mend;

%iterdataset(dataset=design, function=%nrstr(%importBED(&treatment,&units,&rep)));

/* Stack data */

data native_meth;
   set cg_01gy_0u_1_2 cg_01gy_0u_2_2 cg_0gy_0u_1_2 cg_0gy_0u_2_2  cg_1gy_0u_1_2 cg_1gy_0u_2_2
       chg_01gy_0u_1_2 chg_01gy_0u_2_2 chg_0gy_0u_1_2 chg_0gy_0u_2_2  chg_1gy_0u_1_2 chg_1gy_0u_2_2
       chh_01gy_0u_1_2 chh_01gy_0u_2_2 chh_0gy_0u_1_2 chh_0gy_0u_2_2  chh_1gy_0u_1_2 chh_1gy_0u_2_2 

cg_01gy_100u_1_2 cg_01gy_100u_2_2 cg_0gy_100u_1_2 cg_0gy_100u_2_2  cg_1gy_100u_1_2 cg_1gy_100u_2_2
       chg_01gy_100u_1_2 chg_01gy_100u_2_2 chg_0gy_100u_1_2 chg_0gy_100u_2_2  chg_1gy_100u_1_2 chg_1gy_100u_2_2
       chh_01gy_100u_1_2 chh_01gy_100u_2_2 chh_0gy_100u_1_2 chh_0gy_100u_2_2  chh_1gy_100u_1_2 chh_1gy_100u_2_2 ;
run;

data gc_raw;
   set gc_01gy_0u_1_2 gc_01gy_0u_2_2 gc_0gy_0u_1_2 gc_0gy_0u_2_2 gc_1gy_0u_1_2 gc_1gy_0u_2_2 
       gc_01gy_100u_1_2 gc_01gy_100u_2_2 gc_0gy_100u_1_2 gc_0gy_100u_2_2 gc_1gy_100u_1_2 gc_1gy_100u_2_2 ;
run;

data gc_norm;
   set gcn_01gy_0u_1_2 gcn_01gy_0u_2_2 gcn_0gy_0u_1_2 gcn_0gy_0u_2_2 gcn_1gy_0u_1_2 gcn_1gy_0u_2_2 
       gcn_01gy_100u_1_2 gcn_01gy_100u_2_2 gcn_0gy_100u_1_2 gcn_0gy_100u_2_2 gcn_1gy_100u_1_2 gcn_1gy_100u_2_2 ;
run;

proc sort data=gc_raw;
   by treatment units rep chr start_pos stop_pos;
proc sort data=gc_norm;
   by treatment units rep chr start_pos stop_pos;
run;

data gc_all no_raw_oops;
  merge gc_raw (in=in1) gc_norm (in=in2);
  by treatment units rep chr start_pos stop_pos;
  if in1 and in2 then flag_normalized=1; else flag_normalized=0;
  if in1 then output gc_all;
  else output no_raw_oops;
run;

proc sort data=gc_all;
  by chr start_pos stop_pos treatment units rep;
run;

data gc_all2;
 set gc_all;
 norm_factor = total_c / total_c_norm;
 if total_C > 10 then flag_rep_gt10=1;
 else flag_rep_gt10=0;
run;


proc sort data=gc_all2;
  by chr start_pos stop_pos treatment units;
proc means data=gc_all2 noprint;
  by chr start_pos stop_pos treatment units;
  var total_C flag_rep_gt10 flag_normalized;
  output out=gc_summary sum(total_C)=summed_total_C 
          sum(flag_rep_gt10)=num_rep_gt10
          sum(flag_normalized) = num_reps_norm
          max(flag_normalized) = flag_site_norm;
run;


data gc_summary2;
  set gc_summary;
  if summed_total_C >= 10 then flag_sum_ge10=1;
  else flag_sum_ge10=0;
  if summed_total_C > 10 then flag_sum_gt10=1;
  else flag_sum_gt10=0;
run;

proc freq data=gc_summary2 noprint;
  tables _FREQ_*flag_site_norm*num_reps_norm*num_rep_gt10*flag_sum_ge10*flag_sum_gt10 / out=norm_check;
run;

data check;
  set gc_summary2;
  where flag_sum_ge10 = 1;
run;

data endo_site;
  set native_meth;
  keep chr start_pos;
  rename start_pos=stop_pos;
run;

data gc_site;
  set gc_all2;
  keep chr stop_pos;
  rename stop_pos=start_pos;
run;

proc sort data=native_meth;
  by chr start_pos;
proc sort data=gc_site nodup;
  by chr start_pos;
run;

data native_meth2;
  merge native_meth (in=in1) gc_site (in=in2);
  by chr start_pos;
  if in2 then flag_gc=1;
  else flag_gc=0;
  if in1;
run;


proc sort data=gc_all2;
  by chr stop_pos;
proc sort data=endo_site  nodup;
  by chr stop_pos;
run;

data gc_all3;
  merge gc_all2 (in=in1) endo_site (in=in2);
  by chr stop_pos;
  if in2 then flag_endo=1;
  else flag_endo=0;
  if in1;
run;


proc freq data=native_meth2 noprint;
  where total_C >= 10 and chr in ("1","2","3","4","5");
  tables treatment*units*rep*site_type*flag_gc / out=count_endo;
run;

proc freq data=gc_All3 noprint;
  where total_C >= 10 and chr in ("1","2","3","4","5");
  tables treatment*units*rep*site_type*flag_endo / out=count_gc;
run;

proc sort data=native_meth2;
   by treatment units site_type flag_gc chr start_pos stop_pos;
proc means data=native_meth2 noprint;
  by treatment units site_type flag_gc chr start_pos stop_pos;
  var total_C;
  output out=native_meth3 sum(total_C)=total_C;
run;

proc sort data=gc_All3;
   by treatment units site_type flag_endo chr start_pos stop_pos;
proc means data=gc_All3 noprint;
  by treatment units site_type flag_endo chr start_pos stop_pos;
  var total_C;
  output out=gc_All4 sum(total_C)=total_C;
run;


proc freq data=native_meth3 noprint;
  where total_C >= 10 and chr in ("1","2","3","4","5");
  tables treatment*units*site_type*flag_gc / out=count_endo_sum;
run;

proc freq data=gc_All4 noprint;
  where total_C >= 10 and chr in ("1","2","3","4","5");
  tables treatment*units*site_type*flag_endo / out=count_gc_sum;
run;


proc print data=count_endo; run;

/*

                                       site_
 treatment    units             rep    type     flag_gc      COUNT

   01Gy       0U                  1     CG         0        2661572
   01Gy       0U                  1     CG         1         536525
   01Gy       0U                  1     CHG        0        3075190
   01Gy       0U                  1     CHG        1         421349
   01Gy       0U                  1     CHH        0       13735248
   01Gy       0U                  1     CHH        1        1652306
   01Gy       0U                  2     CG         0        2513566
   01Gy       0U                  2     CG         1         492486
   01Gy       0U                  2     CHG        0        2889088
   01Gy       0U                  2     CHG        1         386037
   01Gy       0U                  2     CHH        0       13682592
   01Gy       0U                  2     CHH        1        1569287
   01Gy       100U                1     CG         0        2591129
   01Gy       100U                1     CG         1         530314
   01Gy       100U                1     CHG        0        2996261
   01Gy       100U                1     CHG        1         416677
   01Gy       100U                1     CHH        0       13609233
   01Gy       100U                1     CHH        1        1657645
   01Gy       100U                2     CG         0        2911959
   01Gy       100U                2     CG         1         601205
   01Gy       100U                2     CHG        0        3400193
   01Gy       100U                2     CHG        1         470920
   01Gy       100U                2     CHH        0       15521567
   01Gy       100U                2     CHH        1        1881592
   0Gy        0U                  1     CG         0        2230265
   0Gy        0U                  1     CG         1         486203
   0Gy        0U                  1     CHG        0        2506454
   0Gy        0U                  1     CHG        1         375958
   0Gy        0U                  1     CHH        0        9513492
   0Gy        0U                  1     CHH        1        1300705
   0Gy        0U                  2     CG         0        1889637
   0Gy        0U                  2     CG         1         381207
   0Gy        0U                  2     CHG        0        2133194
   0Gy        0U                  2     CHG        1         297664
   0Gy        0U                  2     CHH        0        9160970
   0Gy        0U                  2     CHH        1        1130664
   0Gy        100U                1     CG         0        2076712
   0Gy        100U                1     CG         1         460860
   0Gy       100U                1     CHG        0       2330347
   0Gy       100U                1     CHG        1        359443
   0Gy       100U                1     CHH        0       8657974
   0Gy       100U                1     CHH        1       1227378
   0Gy       100U                2     CG         0       1537302
   0Gy       100U                2     CG         1        320713
   0Gy       100U                2     CHG        0       1716798
   0Gy       100U                2     CHG        1        253871
   0Gy       100U                2     CHH        0       6961074
   0Gy       100U                2     CHH        1        918669
   1Gy       0U                  1     CG         0       1381982
   1Gy       0U                  1     CG         1        272578
   1Gy       0U                  1     CHG        0       1525944
   1Gy       0U                  1     CHG        1        211796
   1Gy       0U                  1     CHH        0       6715196
   1Gy       0U                  1     CHH        1        811484
   1Gy       0U                  2     CG         0       2069797
   1Gy       0U                  2     CG         1        426639
   1Gy       0U                  2     CHG        0       2340059
   1Gy       0U                  2     CHG        1        330789
   1Gy       0U                  2     CHH        0       9935874
   1Gy       0U                  2     CHH        1       1236074
   1Gy       100U                1     CG         0       1793272
   1Gy       100U                1     CG         1        385661
   1Gy       100U                1     CHG        0       2013615
   1Gy       100U                1     CHG        1        301155
   1Gy       100U                1     CHH        0       8098138
   1Gy       100U                1     CHH        1       1086048
   1Gy       100U                2     CG         0       1815221
   1Gy       100U                2     CG         1        380096
   1Gy       100U                2     CHG        0       2058043
   1Gy       100U                2     CHG        1        299857
   1Gy       100U                2     CHH        0       8448979
   1Gy       100U                2     CHH        1       1106427

*/

proc print data=count_gc; run;

/*

                                      site_    flag_
treatment    units             rep    type      endo     COUNT

  01Gy       0U                  1     GC        0      1452305
  01Gy       0U                  1     GC        1      2655312
  01Gy       0U                  2     GC        0      1376115
  01Gy       0U                  2     GC        1      2488954
  01Gy       100U                1     GC        0      1453485
  01Gy       100U                1     GC        1      2649997
  01Gy       100U                2     GC        0      1639680
  01Gy       100U                2     GC        1      2999682
  0Gy        0U                  1     GC        0      1164361
  0Gy        0U                  1     GC        1      2215483
  0Gy        0U                  2     GC        0      1007260
  0Gy        0U                  2     GC        1      1848152
  0Gy        100U                1     GC        0      1104200
  0Gy        100U                1     GC        1      2100425
  0Gy        100U                2     GC        0       831339
  0Gy        100U                2     GC        1      1532433
  1Gy        0U                  1     GC        0       729953
  1Gy        0U                  1     GC        1      1324850
  1Gy        0U                  2     GC        0      1098067
  1Gy        0U                  2     GC        1      2033067
  1Gy        100U                1     GC        0       974193
  1Gy        100U                1     GC        1      1814191
  1Gy        100U                2     GC        0       992160
  1Gy        100U                2     GC        1      1827574


*/

proc print data=count_endo_sum; run;

/*

                        site_
  treatment    units    type     flag_gc      COUNT

    01Gy       0U        CG         0        3719187
    01Gy       0U        CG         1         756045
    01Gy       0U        CHG        0        4356355
    01Gy       0U        CHG        1         576396
    01Gy       0U        CHH        0       21257029
    01Gy       0U        CHH        1        2381546
    01Gy       100U      CG         0        3785777
    01Gy       100U      CG         1         783705
    01Gy       100U      CHG        0        4447866
    01Gy       100U      CHG        1         596377
    01Gy       100U      CHH        0       21512017
    01Gy       100U      CHH        1        2458634
    0Gy        0U        CG         0        3233143
    0Gy        0U        CG         1         679180
    0Gy        0U        CHG        0        3740869
    0Gy        0U        CHG        1         521658
    0Gy        0U        CHH        0       16252485
    0Gy        0U        CHH        1        1997298
    0Gy        100U      CG         0        3006376
    0Gy        100U      CG         1         643528
    0Gy        100U      CHG        0        3474205
    0Gy        100U      CHG        1         498619
    0Gy        100U      CHH        0       14568674
    0Gy        100U      CHH        1        1869857
    1Gy        0U        CG         0        2949063
    1Gy        0U        CG         1         603746
    1Gy        0U        CHG        0        3403440
    1Gy        0U        CHG        1         464670
    1Gy        0U        CHH        0       15486427
    1Gy        0U        CHH        1        1826942
    1Gy        100U      CG         0        2982333
    1Gy        100U      CG         1         629852
    1Gy        100U      CHG        0        3462176
    1Gy        100U      CHG        1         487358
    1Gy        100U      CHH        0       15135421
    1Gy        100U      CHH        1        1881338


*/

proc print data=count_gc_sum; run;

/*

                       site_    flag_
 treatment    units    type      endo     COUNT

   01Gy       0U        GC        0      2037417
   01Gy       0U        GC        1      3748850
   01Gy       100U      GC        0      2098859
   01Gy       100U      GC        1      3871634
   0Gy        0U        GC        0      1738941
   0Gy        0U        GC        1      3244066
   0Gy        100U      GC        0      1636704
   0Gy        100U      GC        1      3061496
   1Gy        0U        GC        0      1591834
   1Gy        0U        GC        1      2938102
   1Gy        100U      GC        0      1641166
   1Gy        100U      GC        1      3044637


*/

data gc_all3a;
  set gc_all3;
  norm_factor_methylC=methyl_C / methyl_C_norm;
  norm_factor_methylC_2= methyl_C_norm / methyl_C;
  norm_factor_perc=perc_methyl / perc_methyl_norm;

  norm_factor_totalC_2=total_C_norm / total_C;
  norm_factor_perc_2=perc_methyl_norm/ perc_methyl;
run;



proc sort data= gc_all3;
  by  chr start_pos stop_pos treatment units rep;
proc transpose data=gc_all3 out=gc_all3_sbys_total;
  by  chr start_pos stop_pos treatment units;
  id rep;
  var total_C;
run;


proc transpose data=gc_all3 out=gc_all3_sbys_methyl;
  by  chr start_pos stop_pos treatment units;
  id rep;
  var methyl_C;
run;

proc transpose data=gc_all3 out=gc_all3_sbys_perc;
  by  chr start_pos stop_pos treatment units;
  id rep;
  var perc_methyl;
run;

data gc_all3_sbys_total2;
  set gc_all3_sbys_total;
  _1_over_2 = _1 / _2;
  _2_over_1 = _2 / _1;
  _1_total = _1 / (_1 + _2);
  _2_total = _2 / (_1 + _2);
  sum_coverage=sum(_1, _2);
  if sum_coverage >= 10  then flag_10X=1; else flag_10X=0;
run;

data gc2keep;
  set gc_all3_sbys_total2;
  keep chr start_pos stop_pos treatment units sum_coverage flag_10X;
run;


data gc_all3_sbys_methyl2;
  set gc_all3_sbys_methyl;
  _1_over_2 = _1 / _2;
  _2_over_1 = _2 / _1;
  _1_total = _1 / (_1 + _2);
  _2_total = _2 / (_1 + _2);
run;

data gc_all3_sbys_perc2;
  set gc_all3_sbys_perc;
  _1_over_2 = _1 / _2;
  _2_over_1 = _2 / _1;
  _1_total = _1 / (_1 + _2);
  _2_total = _2 / (_1 + _2);
run;

proc sort data=gc_all3_sbys_methyl2;
  by chr start_pos stop_pos treatment units;
proc sort data=gc_all3_sbys_perc2;
  by chr start_pos stop_pos treatment units;
proc sort data=gc2keep;
  by chr start_pos stop_pos treatment units;
run;

data gc_all3_sbys_methyl3;
  merge gc_all3_sbys_methyl2 (in=in1) gc2keep (in=in2);
  by chr start_pos stop_pos treatment units;
  if in1 and in2;
run;

data gc_all3_sbys_perc3;
  merge gc_all3_sbys_perc2 (in=in1) gc2keep (in=in2);
  by chr start_pos stop_pos treatment units;
  if in1 and in2;
run;


proc sort data=gc_all3_sbys_total2;
  by treatment units ;
proc sort data=gc_all3_sbys_methyl3;
  by treatment units ;
proc sort data=gc_all3_sbys_perc3;
  by treatment units ;
run;

proc means data=gc_all3_sbys_total2 noprint;
  by treatment units ;
  where flag_10X = 1 and _1 ne . and _2 ne .;
  var _1_over_2 _2_over_1 _1_total _2_total;
  output out=gc_total_ratio mean=;
run;

proc means data=gc_all3_sbys_methyl3 noprint;
  by treatment units ;
  where flag_10X = 1 and _1 ne . and _2 ne .;
  var _1_over_2 _2_over_1 _1_total _2_total;
  output out=gc_methyl_ratio mean=;
run;

proc means data=gc_all3_sbys_perc3 noprint;
  by treatment units ;
  where flag_10X = 1 and _1 ne . and _2 ne .;
    var _1_over_2 _2_over_1 _1_total _2_total;
  output out=gc_perc_ratio mean=;
run;


proc sort data=gc_All3a;
  by treatment units rep;
proc means data=gc_all3a noprint;
  by treatment units rep;
  where flag_normalized=1;
  var norm_Factor norm_Factor_methylC norm_Factor_perc norm_factor_methylC_2
norm_factor_totalC_2 norm_factor_perc_2;
  output out=norm_factor_distrib
         mean(norm_factor)=mean_norm_factor_totalC
         stddev(norm_factor)=sd_norm_factor_totalC
         mean(norm_factor_methylC)=mean_norm_factor_methylC
         stddev(norm_factor_methylC)=sd_norm_factor_methylC
         mean(norm_factor_methylC_2)=mean_norm_factor_methylC_2
         stddev(norm_factor_methylC_2)=sd_norm_factor_methylC_2
         mean(norm_factor_perc)=mean_norm_factor_perc
         stddev(norm_factor_perc)=sd_norm_factor_perc
         mean(norm_factor_totalC_2)=mean_norm_factor_totalC_2
         stddev(norm_factor_totalC_2)=sd_norm_factor_totalC_2
         mean(norm_factor_perc_2)=mean_norm_factor_perc_2
         stddev(norm_factor_perc_2)=sd_norm_factor_perc_2;
run;

proc print data=norm_factor_distrib;
run;


