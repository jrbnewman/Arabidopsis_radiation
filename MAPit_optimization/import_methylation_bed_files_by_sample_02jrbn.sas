/* Import BED files for each site and sample */

ods listing; ods html close;
libname cold '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname coldloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";

proc datasets lib=work kill noprint;
run;
quit;

%macro importBED(site,temp,units,rep);

   data WORK.input_data    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile "!PATCON/arabidopsis_wgbs_cold/analysis_output/dmap2_output/bed_files_by_sample/&site./&temp._&units._&rep..bed" delimiter='09'x MISSOVER DSD lrecl=32767 ;
      informat VAR1 $3. ;
      informat VAR2 best32. ;
      informat VAR3 best32. ;
      informat VAR4 best32. ;
      informat VAR5 best32. ;
      informat VAR6 best32. ;
      informat VAR7 $1. ;
      informat VAR8 best32. ;
      informat VAR9 $1. ;
      informat VAR10 best32. ;
      informat VAR11 best32. ;
      informat VAR12 $1. ;
      informat VAR13 best32. ;
      informat VAR14 best32. ;
      format VAR1 $3. ;
      format VAR2 best12. ;
      format VAR3 best12. ;
      format VAR4 best12. ;
      format VAR5 best12. ;
      format VAR6 best12. ;
      format VAR7 $1. ;
      format VAR8 best12. ;
      format VAR9 $1. ;
      format VAR10 best12. ;
      format VAR11 best12. ;
      format VAR12 $1. ;
      format VAR13 best12. ;
      format VAR14 best12. ;
     input
                 VAR1
                 VAR2
                 VAR3
                 VAR4
                 VAR5
                 VAR6
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

data meth_&site._&temp._&units._&rep.;
  length sample_id $15.;
  length temperature $4.;
  length units $4.;
  length rep $1.;
  length site_type $3.;
  set input_data;
  temperature="&temp.";
  units="&units.";
  rep="&rep.";
  sample_id=catx("_",temperature,units,rep);
  site_type=upcase("&site.");
  keep VAR1-VAR6 sample_id temperature units rep site_type;
  rename VAR1=chr
         VAR2=pos
         VAR3=pos_end
         VAR4=perc_methyl
         VAR5=total_C
         VAR6=methyl_C;
run;

%mend;

%importBED(CG,22C,0U,1A);
%importBED(CG,22C,0U,2A);
%importBED(CG,4C,0U,1A);
%importBED(CG,4C,0U,2A);
%importBED(CHG,22C,0U,1A);
%importBED(CHG,22C,0U,2A);
%importBED(CHG,4C,0U,1A);
%importBED(CHG,4C,0U,2A);
%importBED(CHH,22C,0U,1A);
%importBED(CHH,22C,0U,2A);
%importBED(CHH,4C,0U,1A);
%importBED(CHH,4C,0U,2A);
%importBED(GC,22C,0U,1A);
%importBED(GC,22C,0U,2A);
%importBED(GC,4C,0U,1A);
%importBED(GC,4C,0U,2A);
%importBED(GC,22C,100U,1A);
%importBED(GC,22C,100U,2A);
%importBED(GC,4C,100U,1A);
%importBED(GC,4C,100U,2A);
%importBED(GC,FANS,0p5U,1A);
%importBED(GC,FANS,0p5U,2A);
%importBED(GC,FANS,1p5U,1A);
%importBED(GC,FANS,1p5U,2A);
%importBED(GC,FANS,5U,1A);
%importBED(GC,FANS,5U,2A);
%importBED(GC,FANS,25U,1A);
%importBED(GC,FANS,25U,2A);



/* Combine */

data methyl_data;
  set meth_: ;
run;


/* Make permenant */

data coldloc.methylation_data_by_rep;
  set methyl_data;
run;

