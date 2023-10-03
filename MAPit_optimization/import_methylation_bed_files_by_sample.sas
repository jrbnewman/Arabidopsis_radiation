/* Import BED files for each site and sample */

ods listing; ods html close;
libname wgbs '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname wgbsloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";

%macro importBED(site,temp,methylase);

    data WORK.input_data    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/TB14/TB14/sandbox/wgbs_sandbox/bed_files/at_cold_&site._all_&temp._&methylase..bed"
delimiter='09'x MISSOVER DSD lrecl=32767 ;
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
                   VAR1 $ 
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

data meth_&site._&temp._&methylase.;
  length condition $15.;
  length temperature $4.;
  length units $4.;
  length site_type $3.;
  set input_data;
  temperature="&temp.";
  units="&methylase.";
  condition=catx("_",temperature,units);
  site_type=upcase("&site.");
  keep VAR1-VAR6 condition temperature units site_type;
  rename VAR1=chr
         VAR2=pos
         VAR3=pos_end
         VAR4=perc_methyl
         VAR5=total_C
         VAR6=methyl_C;
run;
%mend;

%importBED(CG,22C,0U);
%importBED(CG,4C,0U);
%importBED(CHG,22C,0U);
%importBED(CHG,4C,0U);
%importBED(CHH,22C,0U);
%importBED(CHH,4C,0U);

%importBED(GC,22C,0U);
%importBED(GC,4C,0U);
%importBED(GC,22C,100U);
%importBED(GC,4C,100U);

%importBED(GC,FANS,0p5U);
%importBED(GC,FANS,1p5U);
%importBED(GC,FANS,5U);
%importBED(GC,FANS,25U);


/* Combine */

data methyl_data;
  set meth_: ;
run;


/* Make permenant */

data wgbsloc.methylation_data;
  set methyl_data;
run;


/* Export unique positions as a BED file: I want to identify the exact cytosine for each site */

data uniq_sites;
  set methyl_data;
  keep chr pos pos_end;
run;

proc sort data=uniq_sites nodup;
  by chr pos pos_end;
run;

data uniq_sites2;
  set uniq_sites;
  pos2=pos-1;
  pos_end2=pos_end-1;
  length id $100.;
  if chr="" then delete;
  id = catx("_",chr,pos,pos_end);
  keep chr pos2 pos_end2 id;
run;

proc export data=uniq_sites2
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/unique_sites.bed"
     dbms=tab replace;
     putnames=no;
run;

/* pos+1 is the cytosine */



