/* Cold WGBS QC:

(1) Binned rep-to-rep concordance:
        bin every 100bp, calc mean methylation for rep, rep concordance on windows
        (a) reps for same condition
        (b) 22C 100U vs FANS

(2) Estimate bisulfite conversion ratio:
        average % methylation on CHH sites on Chrom Pt by replicate */

libname arab '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


/* import CHH data for 100U samples */

%macro importBED(trt,unit,rep,siteType);

     data WORK.input_data    ;
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



data &siteType._&trt._&unit._&rep.;
  length condition $15.;
  length dose $4.;
  length units $4.;
  length site_type $3.;
  length rep $4.;
  set input_data;
  dose="&trt.";
  rep="&rep.";
  units="&unit.";
  condition=catx("_",dose,units);
  site_type=upcase("&siteType.");
  keep chr  start_pos stop_pos perc_methyl rep total_C methyl_C  condition dose units site_type;
run;

%mend;

%importBED(0Gy,100U,1,chh);
%importBED(0Gy,100U,2,chh);
%importBED(01Gy,100U,1,chh);
%importBED(01Gy,100U,2,chh);
%importBED(1Gy,100U,1,chh);
%importBED(1Gy,100U,2,chh);


data chh_100u;
   set chh_: ;
    where chr="Pt";
run;

data gc_sites;
   set arab.methylation_data_gc;
   where chr="Pt";
   keep chr start_pos stop_pos;
run;

proc sort data=gc_sites nodup;
  by chr start_pos stop_pos;
proc sort data=chh_100u;
  by chr start_pos stop_pos;
run;

data hchh_100u;
  merge chh_100u (in=in1) gc_sites (in=in2);
  by chr start_pos stop_pos;
  if in2 then delete;
run;

proc sort data=hchh_100u;
   by dose units rep;
proc means data=hchh_100u noprint;
   by dose units rep;
   var perc_methyl;
   output out=mean_meth_HCHH_Pt_by_rep mean=mean_methyl_all stddev=sd_methyl_all;
run;

proc print data=mean_meth_HCHH_Pt_by_rep;
run;

/*
                                              mean_methyl_      sd_methyl_
  dose    units    rep    _TYPE_    _FREQ_        all                  all

  01Gy    100U      1        0       22782    0.0452758381    0.0566418001
  01Gy    100U      2        0       22762    0.0493652411    0.0689724972
  0Gy     100U      1        0       23887    0.0284024473    0.0484747028
  0Gy     100U      2        0       23067    0.0287800745    0.0451179667
  1Gy     100U      1        0       23523    0.0350333003    0.0513646279
  1Gy     100U      2        0       23436    0.0276555258    0.0432934304


*/



/* bisulfite conversion ratio */
data meth_data;
  set arab.methylation_data_cg_chg_chh;
  where chr="Pt" and site_type="CHH";
run;

proc sort data=meth_data;
   by treatment units rep;
proc means data=meth_data noprint;
   by treatment units rep;
   var perc_methyl;
   output out=mean_meth_CHH_Pt_by_rep mean=mean_methyl_all stddev=sd_methyl_all;
run;

proc print data=mean_meth_CHH_Pt_by_rep;
run;


/*

                                                           mean_methyl_      sd_methyl_
 treatment    units             rep    _TYPE_    _FREQ_        all                  all

   01Gy        0U                 1       0       26617    0.0340701666    0.0514181153
   01Gy        0U                 2       0       26173    0.0319882793    0.0443313611
   0Gy         0U                 1       0       27822    0.0271857812    0.0480980497
   0Gy         0U                 2       0       27406    0.0292746323    0.0430725983
   1Gy         0U                 1       0       27075    0.0244003449    0.0397991234
   1Gy         0U                 2       0       26838    0.0467796456    0.0948773482




*/

