/* Running model -- by fusion, time: time*treatment, plate as random */

ods listing; ods html close;

libname rs "!HOME/concannon/DTRA/arabidopsis/sas_data";


/* Get list of "on" fusions */
* as we are now doing contrasts that involve mock-only or treatment-only comparisons, I am changing;
* this from "on in all treatments" to "on in any treatment";
* Can always drop fusions later, but before FDR;

data on_gene;
   set rs.arab_flag_gene_on_cpm_gt0;
   if flag_gene_on_mock_cpm0=1 or flag_gene_on_01gy_cpm0=1 or flag_gene_on_1gy_cpm0=1 ;
   keep gene_id;
run;


/* Merge in with count data */


proc sort data=on_gene;
   by gene_id;
proc sort data=cpm_norm_counts_by_xs;
   by gene_id;
run;

data gene_counts_on;
   merge cpm_norm_counts_by_xs (in=in1) on_gene (in=in2);
   by gene_id;
   if in1 and in2;
   log_cpm = log(cpm+1);
run;

data xs_on;
   set rs.arab_flag_transcript_on_gt0;
   if flag_xscript_on_01gy_apn0=1 or  flag_xscript_on_1gy_apn0=1 or flag_xscript_on_Mock_apn0=1;
   keep transcript_id;
run;

proc sort data=gene_counts_on;
  by transcript_id;
proc sort data=xs_on;
  by transcript_id;
run;


data gene_counts_on1;
  merge gene_counts_on (in=in1) xs_on (in=in2);
  by transcript_id;
  if in1 and in2;
run;


proc sort data=gene_counts_on1;
  by gene_id transcript_id;
run;


/* Merge in design file
   Want variables: treatment, time */

/* Model for contrasts */

*Cat together genotype, time and treatment, then find order for contrasts;

data gene_counts_on2;
  set gene_counts_on1;
  length time_trt $20.;
  time_trt = catx('_', time, treatment);
run;

proc freq data=gene_counts_on2 noprint;
  tables time_trt / out=order;
run;

proc print data=order;
run;


*12 levels;
/* ORDER IS:

1_0.1gy
1_1gy
1_Mock
24_0.1gy
24_1gy
24_Mock
3_0.1gy
3_1gy
3_Mock
72_0.1gy
72_1gy
72_Mock

*/

proc sort data=gene_counts_on2;
   by gene_id time;
   run;

ods listing close;


%macro runModels(measure,outname);
ods listing close;
proc mixed data=gene_counts_on2;
  by gene_id;
  class time treatment time_trt transcript_id;
  model log_cpm = time_trt| transcript_id / htype=1;


ods output tests1=tests1 ;
run;
quit;


data sig_ds;
  set tests1;
  where effect ? "*transcr";
run;

proc multtest inpvalues(probf)=sig_ds fdr out=sig_ds_fdr noprint;
run;

data sig_ds_fdr2;
  set sig_ds_fdr;
  if fdr_p = . then flag_fdr05=.;
  else if fdr_p < 0.05 then flag_fdr05=1; 
  else flag_fdr05=0;
run;

ods listing;
proc freq data=sig_ds_fdr2;
  tables flag_fdr05;
run;


/* Make permenant */


data rs.arab_gene_cntrs_tests1_&outname. ;
  set tests1 ;
  run ;

data rs.arab_gene_cntrs_lsmeans_&outname. ;
  set lsmeans ;
  run ;

data rs.arab_gene_cntrs_constr_&outname. ;
  set contrasts ;
  run ;

data rs.arab_gene_cntrs_estim_&outname. ;
  set estimates ;
  run ;


%mend;

%runModels(log_cpm,lcpm);
%runModels(cpm,cpm);
   


