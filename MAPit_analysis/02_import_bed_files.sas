libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
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
 "!PATCON/DTRA/radiation_bed_files/arabidopsis/by_rep/at_rad_&siteType._all_&trt._&unit._&rep..bed"
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

%iterdataset(dataset=design1, function=%nrstr(%importBED(&treatment,&units,&rep)));



/* Import raw and normalized GC 100U methylation:
   I am going to only use the normalized methylation rate here (but keep the total and methyl counts)
*/



%macro importBED(trt,unit,rep);

     data WORK.GC_&trt._&unit._&rep.    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 "!PATCON/DTRA/radiation_bed_files/arabidopsis/by_rep/at_rad_gc_all_&trt._&unit._&rep..bed"
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

 "!PATCON/DTRA/radiation_bed_files/plantNormalization/at_gch_&trt._&unit._rerun_&rep._NormalizedByGroup.BED"
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
       chh_01gy_0u_1_2 chh_01gy_0u_2_2 chh_0gy_0u_1_2 chh_0gy_0u_2_2  chh_1gy_0u_1_2 chh_1gy_0u_2_2 ;
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

/* make permament */

data wgbslocA.methylation_data_CG_CHG_CHH;
   set native_meth;
run;

data wgbslocA.methylation_data_GC;
   set gc_all;
run;






