/* Cold WGBS QC:

(1) Binned rep-to-rep concordance:
        bin every 100bp, calc mean methylation for rep, rep concordance on windows
        (a) reps for same condition
        (b) 22C 100U vs FANS

(2) Estimate bisulfite conversion ratio:
        average % methylation on CHH sites on Chrom Pt by replicate */

libname cold '!PATCON/DTRA/arabidopsis_wgbs_cold/sas_data';
libname coldloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";
ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;



%macro importBED(site,temp,methylase,rep);

    data WORK.input_data    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/TB14/TB14/sandbox/dtra_sandbox/cold_bed_&site./&temp._&methylase._&rep..bed"
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

data meth_&site._&temp._&methylase._&rep.;
  length condition $15.;
  length temperature $4.;
  length units $4.;
  length site_type $3.;
  length rep $3.;
  set input_data;
  temperature="&temp.";
  units="&methylase.";
  rep="&rep.";
  condition=catx("_",temperature,units);
  site_type=upcase("&site.");
  keep VAR1-VAR6 condition temperature units site_type rep;
  rename VAR1=chr
         VAR2=pos
         VAR3=pos_end
         VAR4=perc_methyl
         VAR5=total_C
         VAR6=methyl_C;
run;
%mend;

%importBED(chh,22C,100U,1A);
%importBED(chh,22C,100U,2A);
%importBED(chh,4C,100U,1A);
%importBED(chh,4C,100U,2A);

%importBED(chg,22C,100U,1);
%importBED(chg,22C,100U,2);
%importBED(chg,4C,100U,1);
%importBED(chg,4C,100U,2);

%importBED(cg,22C,100U,1);
%importBED(cg,22C,100U,2);
%importBED(cg,4C,100U,1);
%importBED(cg,4C,100U,2);


data meth_data_100U;
   set meth_cg_: meth_chg_: meth_chh_: ;
run;


data meth_data_0u;
  set coldloc.methylation_data_by_rep;
  where site_type ne "GC" and units="0U";
run;

data meth_data_all;
    set meth_data_0U  meth_data_100U;
   keep temperature units site_type chr pos pos_end perc_methyl total_C;
run;

proc sort data=meth_data_all;
   by temperature units site_type chr pos pos_end;
proc means data=meth_data_all noprint;
   by temperature units site_type chr pos pos_end;
   var perc_methyl total_C;
   output out=meth_data_all2 mean(perc_methyl)=perc_methyl sum(total_C)=total_C;
run;

data gc_sites;
  set coldloc.methylation_data_by_rep;
  where site_type="GC";
  keep chr pos;
run;

proc sort data=gc_sites nodup;
  by chr pos;
proc sort data=meth_data_all2;
   by chr pos;
run;

data meth_data_all3;
  merge meth_data_all2 gc_sites (in=in2);
  by chr pos;
  if in2 then delete;
run;



data meth_data_all4;
  set meth_data_all3;
  chr_bin=int(pos_end/100) + 1;
  if total_C < 100 then delete;
run;

proc sort data=meth_data_all4;
  by temperature units chr chr_bin ;
proc means data=meth_data_all4 noprint;
  by  temperature units chr chr_bin ;
  var perc_methyl;
  output out=mean_methyl_by_bin_all mean=;
run;


/* Calc correlations */
ods graphics / ANTIALIASMAX=50000000;

%macro corrPlot(temp1,temp2,units1,units2);

data sample1_all;
  set mean_methyl_by_bin_all;
  where temperature = "&temp1." and units="&units1.";
  keep chr chr_bin perc_methyl;
  rename perc_methyl=perc_meth_&temp1._&units1.;
run;


data sample2_all;
  set mean_methyl_by_bin_all;
  where temperature = "&temp2." and units="&units2.";
  keep chr chr_bin perc_methyl;
  rename perc_methyl=perc_meth_&temp2._&units2.;
run;


proc sort data=sample1_all;
  by chr chr_bin;
proc sort data=sample2_all;
  by chr chr_bin;
run;

data sample_1v2_all;
  merge sample1_all (in=in1) sample2_all (in=in2);
  by chr chr_bin;
run;

title "Correlation  &temp1. &units1. &temp2. &units2.  (100X sites) ";
ods text="Correlation  &temp1. &units1. &temp2. &units2.  (100X sites)";

proc corr data=sample_1v2_all pearson spearman polychoric polyserial kendall;
  var perc_meth_&temp1._&units1. perc_meth_&temp2._&units2.;
run;

proc sgplot data=sample_1v2_all;
   scatter x=perc_meth_&temp1._&units1. y= perc_meth_&temp2._&units2.;
run;

%mend;


%corrPlot(22C,22C,0U,100U);
%corrPlot(4C,4C,0U,100U);


/*

 Variable                     N          Mean       Std Dev        Median       Minimum       Maximum

 perc_meth_22C_0U          4412       0.10609       0.14153       0.04469             0       0.91928
 perc_meth_22C_100U        4197       0.08656       0.11459       0.03708             0       0.92657


                                  Pearson Correlation Coefficients
                                     Prob > |r| under H0: Rho=0
                                        Number of Observations

                                                       perc_         perc_
                                                       meth_         meth_
                                                      22C_0U      22C_100U

                            perc_meth_22C_0U         1.00000       0.93728
                                                                    <.0001
                                                        4412          4064



                                          Simple Statistics

 Variable                    N          Mean       Std Dev        Median       Minimum       Maximum

 perc_meth_4C_0U          4545       0.08506       0.12717       0.02335             0       0.94283
 perc_meth_4C_100U        4215       0.08527       0.11281       0.03693             0       0.90770


                                  Pearson Correlation Coefficients
                                     Prob > |r| under H0: Rho=0
                                       Number of Observations

                                                      perc_         perc_
                                                      meth_      meth_4C_
                                                      4C_0U          100U

                            perc_meth_4C_0U         1.00000       0.92695
                                                                   <.0001
                                                       4545          4154



*/

