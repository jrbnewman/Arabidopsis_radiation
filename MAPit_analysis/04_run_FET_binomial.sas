libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

*libname wgbslocA '/blue/concannon/share/jnewman/mingqi_arab/sas_data';
ods listing;
ods html close;
options nosyntaxcheck;



proc datasets lib=work kill noprint;
run;
quit;


/* DMCs:

(1) Subset only sites with at least 10 mapped reads in both conditions, and methylated in at least one condition
(2) Run a Fisher's exact test on each site
(3) FDR correction on FETs
(4) Calculate methylation difference
(5) Flag FDR 5% and diff >10%

*/

/* Prep data by condition -- may speed things up a little */
%macro byChrom(chrom);
%macro prepData(trt);

data counts;
   set wgbslocA.methylation_data_cg_chg_chh;
   if chr = "Mt" then delete;
   if chr = "Pt" then delete;
   where units="0U" and treatment="&trt." and chr="&chrom.";
  rename methyl_C=methyl_C_&trt. total_C=total_C_&trt. ;
run;  

data cov;
   set wgbslocA.flag_coverage_10x_cg_chg_chh;
   where flag_&trt._coverage_ge_10x=1;
   keep site_type chr start_pos stop_pos;
run;

proc sort data=cov;
  by site_type chr start_pos stop_pos;
proc sort data=counts;
  by site_type chr start_pos stop_pos;
run;

data counts_&trt.;
   merge counts (in=in1) cov (in=in2);
   by site_type chr start_pos stop_pos;
   if in1 and in2;
run;

%mend;

%prepData(0Gy);
%prepData(01Gy);
%prepData(1Gy);


/* Stack data, prep, and run FETs */
%macro runBinomial(condit1, condit2);

data sites2keep1;
  set counts_&condit1.;
  where methyl_C_&condit1.  > 0;
  keep site_type chr start_pos stop_pos ;
run;
 data sites2keep2;
   set counts_&condit2.;
   where methyl_C_&condit2.  > 0;
   keep site_type chr start_pos stop_pos;
 run;




proc sort data=sites2keep1 nodup;
 by site_Type chr start_pos stop_pos;
proc sort data=sites2keep2  nodup;
 by site_Type chr start_pos stop_pos;
run;

data sites2keep;
  merge sites2keep1 (in=in1) sites2keep2 (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;


data counts_for_binomial;
  set counts_&condit1. (in=in1) counts_&condit2. (in=in2);
  if in1 then do;
    methyl_C = methyl_C_&condit1.;
    total_C = total_C_&condit1.;
    end;
  if in2 then do;
    methyl_C = methyl_C_&condit2.;
    total_C = total_C_&condit2.;
    end;
run;


proc sort data=counts_for_binomial;
  by site_Type chr start_pos stop_pos;
proc sort data=sites2keep;
  by site_type chr start_pos stop_pos;
run;

data counts_for_binomial2;
  merge sites2keep (in=in1) counts_for_binomial (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

proc printto log="!HOME/binomial_DMC_arab_&condit1._&condit2._&chrom..log";
run;


ods listing close;
proc logistic data=counts_for_binomial2;
   by site_type chr start_pos stop_pos;
   class treatment;
   model methyl_C / total_C = treatment / link=logit ;
   output out=gmxout pred=pred ;
   ods output modelanova=anova_out FitStatistics=fitstats_out;
run;
quit;

proc printto;
run;


data wgbsA.binom_DMC_&condit1._&condit2._&chrom.;
   set anova_out;
run;

data wgbsA.binft_DMC_&condit1._&condit2._&chrom.;
   set fitstats_out;
run;


%mend;


%runBinomial(0Gy, 01Gy);
%runBinomial(0Gy, 1Gy);



%macro prepData(trt,units);

data counts;
   set wgbslocA.methylation_data_gc;
      if chr = "Mt" then delete;
         if chr = "Pt" then delete;
	    where units="&units." and treatment="&trt." and chr="&chrom.";
	    rename methyl_C=methyl_C_&trt._&units. total_C=total_C_&trt._&units.;
	    run;

			data cov;
			   set wgbslocA.flag_coverage_10x_gc;
			      where flag_&trt._&units._coverage_ge_10x=1;
			         keep site_type chr start_pos stop_pos;
				 run;

				 proc sort data=cov;
				   by site_type chr start_pos stop_pos;
				   proc sort data=counts;
				     by site_type chr start_pos stop_pos;
				     run;

				     data counts_&trt._&units.;
				        merge counts (in=in1) cov (in=in2);
					   by site_type chr start_pos stop_pos;
					      if in1 and in2;
					      run;

%mend;



%prepData(0Gy,0U);
%prepData(01Gy,0U);
%prepData(1Gy,0U);
%prepData(0Gy,100U);
%prepData(01Gy,100U);
%prepData(1Gy,100U);


/* Stack data, prep, and run FETs */



/* Run binomial DAC */

%macro runBinomial(condit1, condit2);

data counts_for_binomial;
  set counts_&condit1._0U (in=in1) counts_&condit1._100U (in=in2)
      counts_&condit2._0U (in=in3) counts_&condit2._100U (in=in4);
  if in1 then do;
    methyl_C = methyl_C_&condit1._0U;
    total_C = total_C_&condit1._0U;
    end;
  if in2 then do;
    methyl_C = methyl_C_&condit1._100U;
    total_C = total_C_&condit1._100U;
    end;
  if in3 then do;
    methyl_C = methyl_C_&condit2._0U;
    total_C = total_C_&condit2._0U;
    end;
  if in4 then do;
    methyl_C = methyl_C_&condit2._100U;
    total_C = total_C_&condit2._100U;
    end;
run;


proc sort data=counts_for_binomial;
  by site_Type chr start_pos stop_pos;
run;


data sites2keep1;
  set counts_&condit1._0U;
  where methyl_C_&condit1._0U > 0;
  keep site_type chr start_pos stop_pos ;
run;
 data sites2keep2;
   set counts_&condit2._0U;
  where methyl_C_&condit2._0U > 0;
   keep site_type chr start_pos stop_pos;
 run;

data sites2keep3;
  set counts_&condit1._100U;
  where methyl_C_&condit1._100U > 0;
  keep site_type chr start_pos stop_pos ;
run;

 data sites2keep4;
   set counts_&condit2._100U;
  where methyl_C_&condit2._100U > 0;
   keep site_type chr start_pos stop_pos ;
 run;

proc sort data=sites2keep1 nodup;
 by site_Type chr start_pos stop_pos;
proc sort data=sites2keep2 nodup;
 by site_Type chr start_pos stop_pos;
proc sort data=sites2keep3 nodup;
 by site_Type chr start_pos stop_pos;
proc sort data=sites2keep4 nodup;
 by site_Type chr start_pos stop_pos;
run;

data sites2keep;
  merge sites2keep1 (in=in1) sites2keep2 (in=in2)
        sites2keep3 (in=in3) sites2keep4 (in=in4);
  by site_type chr start_pos stop_pos;
  if in1 and in2 and in3 and in4;
run;

proc sort data=sites2keep nodup;
  by site_type chr start_pos stop_pos;
run;

data counts_for_binomial2;
  merge sites2keep (in=in1) counts_for_binomial (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

proc printto log="!HOME/binomial_DAC_arab_&condit1._&condit2._&chrom..log";
run;


ods listing close;
proc logistic data=counts_for_binomial2;
   by site_type chr start_pos stop_pos;
   class treatment units;
   model methyl_C / total_C = treatment units treatment*units / link=logit ;
   output out=gmxout pred=pred ;
   ods output modelanova=anova_out FitStatistics=fitstats_out;
run;
quit;

proc printto ;
run;


data wgbsA.binom_DAC_&condit1._&condit2._&chrom.;
   set anova_out;
run;

data wgbsA.binft_DAC_&condit1._&condit2._&chrom.;
   set fitstats_out;
run;



%mend;

%runBinomial(0Gy, 01Gy);
%runBinomial(0Gy, 1Gy);

%mend;


%byChrom(1);
%byChrom(2);
%byChrom(3);
%byChrom(4);
%byChrom(5);



