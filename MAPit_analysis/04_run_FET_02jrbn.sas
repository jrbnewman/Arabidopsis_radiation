*libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
*libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

libname wgbslocA '/blue/concannon/share/jnewman/mingqi_arab/sas_data';
ods listing;
ods html close;

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
%let chrom=4;

/* Prep data by condition -- may speed things up a little */
%macro byChrom(chrom);
%macro prepData(trt);

data counts;
   set wgbslocA.methylation_data_cg_chg_chh;
   if chr = "Mt" then delete;
   if chr = "Pt" then delete;
   where units="0U" and treatment="&trt." and chr="&chrom.";
run;  

proc sort data=counts;
   by chr start_pos stop_pos site_type treatment units;
proc means data=counts noprint;
   by chr start_pos stop_pos site_type treatment units;
   var total_C methyl_C perc_methyl;
   output out=counts2 sum(total_C)=total_C_&trt. sum(methyl_C)=methyl_C_&trt. mean(perc_methyl)=perc_methyl_&trt.;
run;

data cov;
   set wgbslocA.flag_coverage_10x_cg_chg_chh;
   where flag_&trt._coverage_ge_10x=1;
   keep site_type chr start_pos stop_pos;
run;

proc sort data=cov;
  by site_type chr start_pos stop_pos;
proc sort data=counts2;
  by site_type chr start_pos stop_pos;
run;

data counts_&trt.;
   merge counts2 (in=in1) cov (in=in2);
   by site_type chr start_pos stop_pos;
   if in1 and in2;
run;

%mend;

%prepData(0Gy);
%prepData(01Gy);
%prepData(1Gy);


/* Stack data, prep, and run FETs */

%macro runFET(condit1,condit2);


%let condit1=0Gy;
%let condit2=01Gy;

proc sort data=counts_&condit1.;
   by site_type chr start_pos stop_pos;
proc sort data=counts_&condit2.;
   by site_type chr start_pos stop_pos;
run;

data methyl_data;
  merge counts_&condit1. (in=in1) counts_&condit2. (in=in2);
   by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

   data methyl_data2;
    set methyl_data;
    unmethyl_C_&condit1. = total_C_&condit1. - methyl_C_&condit1.;
    unmethyl_C_&condit2. = total_C_&condit2. - methyl_C_&condit2.;
 run;


 data methyl_data_fet;
 set methyl_data2;
  length condition $10.;
  by site_type chr start_pos stop_pos;
  i=0;
   j=0;
   k=0;
    l=0;
    if methyl_C_&condit1. > 0 then do while (i < methyl_C_&condit1.);
     condition="&condit1.";
     flag_methylated=1;
      i + 1;
     output;
    end;

    if unmethyl_C_&condit1. > 0 then do while (j < unmethyl_C_&condit1.);
  condition="&condit1.";
   flag_methylated=0;
     j + 1;
  output;
    end;

    if methyl_C_&condit2. > 0 then do while (k < methyl_C_&condit2.);
  condition="&condit2.";
   flag_methylated=1;
     k + 1;
  output;
    end;

    if unmethyl_C_&condit2. > 0 then do while (l < unmethyl_C_&condit2.);
  condition="&condit2.";
   flag_methylated=0;
     l + 1;
  output;
    end;

    run;


proc printto log="/blue/concannon/share/jnewman/mingqi_arab/FET_arab_&condit1._&condit2..log";
run;


proc sort data=methyl_data_fet;
    by site_type chr start_pos stop_pos;
run;

    proc freq data=methyl_data_fet noprint;
    by site_type chr start_pos stop_pos;
 exact fisher;
 tables condition * flag_methylated / chisq cmh ;
  output out=fet_dmc fisher chisq cmh;
  run;


proc printto ;
run;

/* Make permament */

data wgbslocA.results_by_DMC_&condit1._&condit2._&chrom.;
   set fet_dmc;
run;

%mend;

%runFET(0Gy, 01Gy);
%runFET(0Gy, 1Gy);


%macro prepData(trt,units);

data counts;
   set wgbslocA.methylation_data_gc;
      if chr = "Mt" then delete;
         if chr = "Pt" then delete;
	    where units="&units." and treatment="&trt." and chr="&chrom.";
	    run;

	    proc sort data=counts;
	       by chr start_pos stop_pos site_type treatment units;
	       proc means data=counts noprint;
	          by chr start_pos stop_pos site_type treatment units;
		     var total_C methyl_C perc_methyl;
		        output out=counts2 sum(total_C)=total_C_&trt._&units. sum(methyl_C)=methyl_C_&trt._&units. mean(perc_methyl)=perc_methyl_&trt._&units.;
			run;

			data cov;
			   set wgbslocA.flag_coverage_10x_gc;
			      where flag_&trt._&units._coverage_ge_10x=1;
			         keep site_type chr start_pos stop_pos;
				 run;

				 proc sort data=cov;
				   by site_type chr start_pos stop_pos;
				   proc sort data=counts2;
				     by site_type chr start_pos stop_pos;
				     run;

				     data counts_&trt._&units.;
				        merge counts2 (in=in1) cov (in=in2);
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

%macro runFET(condit1, condit2);
proc sort data=counts_&condit1.;
by site_type chr start_pos stop_pos;
proc sort data=counts_&condit2.;
by site_type chr start_pos stop_pos;
run;

data methyl_data;
merge counts_&condit1. (in=in1) counts_&condit2. (in=in2);
by site_type chr start_pos stop_pos;
if in1 and in2;
run;

data methyl_data2;
set methyl_data;
unmethyl_C_&condit1. = total_C_&condit1. - methyl_C_&condit1.;
unmethyl_C_&condit2. = total_C_&condit2. - methyl_C_&condit2.;
 run;


data methyl_data_fet;
set methyl_data2;
length condition $10.;
by site_type chr start_pos stop_pos;
i=0;
j=0;
k=0;
l=0;
if methyl_C_&condit1. > 0 then do while (i < methyl_C_&condit1.);
condition="&condit1.";
flag_methylated=1;
i + 1;
output;
end;

if unmethyl_C_&condit1. > 0 then do while (j < unmethyl_C_&condit1.);
condition="&condit1.";
flag_methylated=0;
j + 1;
output;
end;

if methyl_C_&condit2. > 0 then do while (k < methyl_C_&condit2.);
condition="&condit2.";
flag_methylated=1;
k + 1;
output;
end;

if unmethyl_C_&condit2. > 0 then do while (l < unmethyl_C_&condit2.);
condition="&condit2.";
flag_methylated=0;
l + 1;
output;
end;
run;


proc printto log="/blue/concannon/share/jnewman/mingqi_arab/FET_DAC_arab_&condit1._&condit2..log";
run;


proc sort data=methyl_data_fet;
    by site_type chr start_pos stop_pos;
run;
proc freq data=methyl_data_fet noprint;
by site_type chr start_pos stop_pos;
exact fisher;
tables condition * flag_methylated / chisq cmh ;
output out=fet_dmc fisher chisq cmh;
run;


proc printto ;
run;

/* Make permament */

data wgbslocA.results_DAC_&condit1._&condit2._&chrom.;
 set fet_dmc;
run;

%mend;

%runFET(0Gy_100U,0Gy_0U);
%runFET(01Gy_100U,01Gy_0U);
%runFET(1Gy_100U,1Gy_0U);


/* 2x2x2 DAC test */




%macro runFET(condit1, condit2);
%let condit1=0Gy;
%let condit2=01Gy;

proc sort data=counts_&condit1._0U;
by site_type chr start_pos stop_pos;
proc sort data=counts_&condit1._100U;
by site_type chr start_pos stop_pos;
proc sort data=counts_&condit2._0U;
by site_type chr start_pos stop_pos;
proc sort data=counts_&condit2._100U;
by site_type chr start_pos stop_pos;
run;

data methyl_data_1;
merge counts_&condit1._100U (in=in1) counts_&condit1._0U (in=in2);
by site_type chr start_pos stop_pos;
if in1 and in2;
run;

data methyl_data_2;
merge counts_&condit2._100U (in=in1) counts_&condit2._0U (in=in2);
by site_type chr start_pos stop_pos;
if in1 and in2;
run;


data methyl_data2_1;
set methyl_data_1;
unmethyl_C_&condit1._100U = total_C_&condit1._100U - methyl_C_&condit1._100U;
unmethyl_C_&condit1._0U = total_C_&condit1._0U - methyl_C_&condit1._0U;
 run;

data methyl_data2_2;
set methyl_data_2;
unmethyl_C_&condit2._100U = total_C_&condit2._100U - methyl_C_&condit2._100U;
unmethyl_C_&condit2._0U = total_C_&condit2._0U - methyl_C_&condit2._0U;
 run;


data sites1;
  set methyl_Data2_1;
   keep site_type chr start_pos stop_pos;
run;

data sites2;
  set methyl_Data2_2;
   keep site_type chr start_pos stop_pos;
run;

proc sort data=sites1 nodup;
  by site_type chr start_pos stop_pos;
proc sort data=sites2 nodup;
  by site_type chr start_pos stop_pos;
run;


data sites2keep;
  merge sites1 (in=in1) sites2 (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;


proc sort data=sites2keep;
  by site_type chr start_pos stop_pos;
proc sort data=methyl_data2_1;
  by site_type chr start_pos stop_pos;
proc sort data=methyl_data2_2;
  by site_type chr start_pos stop_pos;
run;


data methyl_data2_1a;
  merge sites2keep (in=in1) methyl_data2_1 (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;

data methyl_data2_2a;
  merge sites2keep (in=in1) methyl_data2_2 (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;



data methyl_data_fet1;
set methyl_data2_1a;
length condition $10.;
length unit $4.;
by site_type chr start_pos stop_pos;
i=0;
j=0;
k=0;
l=0;
if methyl_C_&condit1._100U > 0 then do while (i < methyl_C_&condit1._100U);
condition="&condit1.";
unit="100U";
flag_methylated=1;
i + 1;
output;
end;

if unmethyl_C_&condit1._100U > 0 then do while (j < unmethyl_C_&condit1._100U);
condition="&condit1.";
unit="100U";
flag_methylated=0;
j + 1;
output;
end;

if methyl_C_&condit1._0U > 0 then do while (k < methyl_C_&condit1._0U);
condition="&condit1.";
unit="0U";
flag_methylated=1;
k + 1;
output;
end;

if unmethyl_C_&condit1._0U > 0 then do while (l < unmethyl_C_&condit1._0U);
condition="&condit1.";
unit="0U";
flag_methylated=0;
l + 1;
output;
end;
run;



data methyl_data_fet2;
set methyl_data2_2a;
length condition $10.;
length unit $4.;
by site_type chr start_pos stop_pos;
i=0;
j=0;
k=0;
l=0;
if methyl_C_&condit2._100U > 0 then do while (i < methyl_C_&condit2._100U);
condition="&condit2.";
unit="100U";
flag_methylated=1;
i + 1;
output;
end;

if unmethyl_C_&condit2._100U > 0 then do while (j < unmethyl_C_&condit2._100U);
condition="&condit2.";
unit="100U";
flag_methylated=0;
j + 1;
output;
end;

if methyl_C_&condit2._0U > 0 then do while (k < methyl_C_&condit2._0U);
condition="&condit2.";
unit="0U";
flag_methylated=1;
k + 1;
output;
end;

if unmethyl_C_&condit2._0U > 0 then do while (l < unmethyl_C_&condit2._0U);
condition="&condit2.";
unit="0U";
flag_methylated=0;
l + 1;
output;
end;
run;


data methyl_data_fet;
  set methyl_data_fet1 methyl_data_fet2;
run;


proc printto log="/blue/concannon/share/jnewman/mingqi_arab/FET_DAC3_arab_&condit1._&condit2..log";
run;


proc freq data=methyl_data_fet noprint;
by site_type chr start_pos stop_pos;
exact fisher;
tables condition * unit * flag_methylated / chisq cmh ;
output out=fet_dmc fisher chisq cmh;
run;


proc printto ;
run;

/* Make permament */

data wgbslocA.results_DAC3_&condit1._&condit2._&chrom.;
 set fet_dmc;
run;

%mend;

%runFET(01Gy,0Gy);
%runFET(1Gy,0Gy);


%mend;


%byChrom(1);
%byChrom(2);
%byChrom(3);
%byChrom(4);
%byChrom(5);



