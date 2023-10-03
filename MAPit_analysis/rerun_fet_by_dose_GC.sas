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

%macro runFET(condit1, condit2);

proc sort data=counts_&condit1.;
by site_type chr start_pos stop_pos;
proc sort data=counts_&condit2.;
by site_type chr start_pos stop_pos;
run;

data methyl_data_1;
merge counts_&condit1. (in=in1) counts_&condit2. (in=in2);
by site_type chr start_pos stop_pos;
if in1 and in2;
run;

data methyl_data2_1;
set methyl_data_1;
unmethyl_C_&condit1. = total_C_&condit1. - methyl_C_&condit1.;
unmethyl_C_&condit2. = total_C_&condit2. - methyl_C_&condit2.;
 run;


proc sort data=methyl_data2_1;
  by site_type chr start_pos stop_pos;

proc transpose data=methyl_data2_1 out=methyl_data2_1_sbys;
  by site_type chr start_pos stop_pos;
  var unmethyl_C_&condit1. methyl_C_&condit1. 
      unmethyl_C_&condit2. methyl_C_&condit2.;
run;


data methyl_data2_1_sbys2;
  set methyl_data2_1_sbys;
  length condition $10.;
  if _NAME_ = "unmethyl_C_&condit1." then do;
        condition = "&condit1.";
        flag_methylated=0;
        end;
  if _NAME_ = "methyl_C_&condit1." then do;
        condition = "&condit1.";
        flag_methylated=1;
        end;
  if _NAME_ = "unmethyl_C_&condit2." then do;
        condition = "&condit2.";
        flag_methylated=0;
        end;
  if _NAME_ = "methyl_C_&condit2." then do;
        condition = "&condit2.";
        flag_methylated=1;
        end;
  rename col1=count;
run;



proc printto log="!HOME/FETbyDose_DAC_arab_&condit1._&condit2._&chrom..log";
run;


proc sort data=methyl_data2_1_sbys2;
  by site_type chr start_pos stop_pos condition flag_methylated;
run;

proc freq data=methyl_data2_1_sbys2 noprint;
by site_type chr start_pos stop_pos;
    weight count;
    exact fisher;
   tables condition * flag_methylated / chisq cmh fisher;
   output out=fet_dmc fisher chisq cmh agree ;
    run;


proc printto ;
run;

/* Make permament */

data wgbsA.FET_DAC_&condit1._&condit2._&chrom.;
 set fet_dmc;
run;

%mend;

%runFET(0Gy_100U,0Gy_0U);
%runFET(01Gy_100U,01Gy_0U);
%runFET(1Gy_100U,1Gy_0U);



%mend;


%byChrom(1);
%byChrom(2);
%byChrom(3);
%byChrom(4);
%byChrom(5);



