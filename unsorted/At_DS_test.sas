ods listing; ods html close;
libname arabRNA "!HOME/concannon/DTRA/arabidopsis/sas_data";

/* Normalize gene read counts using "copies per million" CPM: reads aligning to gene / 1,000,000
   [reads aligning to gene] * 10^6 / [total number of aligned reads in one sample] */

data counts_by_isoform;
  set arabRNA.counts_by_isoform;
run;

/* sum expected_count by sample */

proc sort data=counts_by_isoform;
  by sample_id;
proc means data=counts_by_isoform noprint;
  by sample_id;
  var expected_count;
  output out=total_counts_per_sample sum=total_counts;
run;

proc sort data=total_counts_per_sample;
  by sample_id;
proc sort data=counts_by_isoform;
  by sample_id;
run;

data counts_by_isoform2;
  merge counts_by_isoform (in=in1) total_counts_per_sample (in=in2);
  by sample_id;
  if in1 and in2;
run;

/* Calculate CPM from expected_counts */

data calc_cpm;
  set counts_by_isoform2;
  cpm = (expected_count * 1000000)/total_counts;
run;

data cpm_norm_counts_by_xs;
  set calc_cpm;
  drop _TYPE_ _FREQ_;
run;



/* Get list of "on" fusions */
* as we are now doing contrasts that involve mock-only or treatment-only comparisons, I am changing;
* this from "on in all treatments" to "on in any treatment";
* Can always drop fusions later, but before FDR;

data on_gene;
   set arabRNA.arab_flag_gene_on_cpm_gt0;
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
   set arabRNA.arab_flag_transcript_on_gt0;
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


/* Split analysis so each comparison is a separate analysis. This is NOT ideal, but I can't think of
another way to model a dynamically-size interaction terms like treatment*time*transcript */


data counts_mock_1 counts_mock_3 counts_mock_24 counts_mock_72
     counts_10_1 counts_10_3 counts_10_24 counts_10_72
     counts_100_1 counts_100_3 counts_100_24 counts_100_72;
   set gene_counts_on2;
   if treatment="Mock" and time=1 then output counts_mock_1;
   if treatment="Mock" and time=3 then output counts_mock_3;
   if treatment="Mock" and time=24 then output counts_mock_24;
   if treatment="Mock" and time=72 then output counts_mock_72;

   if treatment="0.1gy" and time=1 then output counts_10_1;
   if treatment="0.1gy" and time=3 then output counts_10_3;
   if treatment="0.1gy" and time=24 then output counts_10_24;
   if treatment="0.1gy" and time=72 then output counts_10_72;

   if treatment="1gy" and time=1 then output counts_100_1;
   if treatment="1gy" and time=3 then output counts_100_3;
   if treatment="1gy" and time=24 then output counts_100_24;
   if treatment="1gy" and time=72 then output counts_100_72;
run;

data counts_mock_10_1;
   set counts_mock_1 (in=in1) counts_10_1  (in=in2);
   length comparison $32.;
   comparison = "10_vs_Mock_1h";
run;

data counts_mock_10_3;
   set counts_mock_3 (in=in1) counts_10_3  (in=in2);
   length comparison $32.;
   comparison = "10_vs_Mock_3h";
run;

data counts_mock_10_24;
   set counts_mock_24 (in=in1) counts_10_24 (in=in2);
   length comparison $32.;
   comparison = "10_vs_Mock_24h";
run;

data counts_mock_10_72;
   set counts_mock_72 (in=in1) counts_10_72  (in=in2);
   length comparison $32.;
   comparison = "10_vs_Mock_72h";
run;


data counts_mock_100_1;
   set counts_mock_1 (in=in1) counts_100_1  (in=in2);
   length comparison $32.;
   comparison = "100_vs_Mock_1h";
run;

data counts_mock_100_3;
   set counts_mock_3 (in=in1) counts_100_3  (in=in2);
   length comparison $32.;
   comparison = "100_vs_Mock_3h";
run;

data counts_mock_100_24;
   set counts_mock_24 (in=in1) counts_100_24 (in=in2);
   length comparison $32.;
   comparison = "100_vs_Mock_24h";
run;

data counts_mock_100_72;
   set counts_mock_72 (in=in1) counts_100_72  (in=in2);
   length comparison $32.;
   comparison = "100_vs_Mock_72h";
run;

data counts_all_stack;
  set counts_mock_10_1 counts_mock_10_3 counts_mock_10_24 counts_mock_10_72
       counts_mock_100_1 counts_mock_100_3 counts_mock_100_24 counts_mock_100_72;
run;

proc sort data=counts_all_stack;
  by comparison gene_id treatment transcript_id ;
run;



ods listing close;
proc mixed data=counts_all_stack;
  by comparison gene_id;
  class treatment  transcript_id;
  model log_cpm = treatment | transcript_id / htype=1;
ods output tests1=tests1 ;
run;
quit;


data sig_ds;
  set tests1;
  where effect ? "*transcr";
run;

proc multtest inpvalues(probf)=sig_ds fdr out=sig_ds_fdr_by_comp noprint;
by comparison ;
run;

data sig_ds_fdr_by_comp2;
  set sig_ds_fdr_by_comp;
  if fdr_p = . then flag_fdr05=.;
  else if fdr_p < 0.05 then flag_fdr05=1; 
  else flag_fdr05=0;
run;

ods listing;
proc freq data=sig_ds_fdr_by_comp2;
  tables comparison*flag_fdr05;
run;

data sig_ds_fdr_by_comp3;
  set sig_ds_fdr_by_comp2;
  where flag_fdr05=1;
run;

proc sort data=sig_ds_fdr_by_comp3;
  by gene_id comparison;
run;

proc freq data=sig_ds_fdr_by_comp3 noprint;
  tables gene_id / out=gene_cnt;
run;

data gene_cnt2;
  set gene_cnt;
  where count > 1;
  keep gene_id;
run;
 
proc print data=gene_cnt2;
run;

proc sort data=gene_cnt2;
  by gene_id;
proc sort data=sig_ds_fdr_by_comp3;
  by gene_id;
run;

data check;
  merge gene_cnt2 (in=in1) sig_ds_fdr_by_comp3 (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc print data=check;
run;


                                                                                       Num     Den                                     flag_
    gene_id            COUNT    PERCENT      comparison              Effect             DF      DF     FValue     ProbF      fdr_p     fdr05

    AT1G16480            2      1.06383    100_vs_Mock_1h     treatment*transcript       1       8      90.11    <.0001    0.008711      1
    AT1G16480            2      1.06383    10_vs_Mock_1h      treatment*transcript       1       8      58.05    <.0001    0.033655      1
    AT1G63770            2      1.06383    100_vs_Mock_24h    treatment*transcript       5      24      22.18    <.0001    0.000245      1
    AT1G63770            2      1.06383    10_vs_Mock_24h     treatment*transcript       5      24      24.65    <.0001    0.000029      1
    AT1G69840            2      1.06383    100_vs_Mock_1h     treatment*transcript       5      24      12.09    <.0001    0.006197      1
    AT1G69840            2      1.06383    10_vs_Mock_1h      treatment*transcript       5      24       8.28    0.0001    0.035262      1
    AT2G06025            2      1.06383    100_vs_Mock_1h     treatment*transcript       6      28       7.65    <.0001    0.026365      1
    AT2G06025            2      1.06383    10_vs_Mock_1h      treatment*transcript       6      28       6.52    0.0002    0.048550      1
    AT2G22990            3      1.59574    100_vs_Mock_24h    treatment*transcript       5      24      19.62    <.0001    0.000392      1
    AT2G22990            3      1.59574    10_vs_Mock_24h     treatment*transcript       5      24      30.03    <.0001    0.000006      1
    AT2G22990            3      1.59574    10_vs_Mock_3h      treatment*transcript       5      24      17.57    <.0001    0.000436      1
    AT2G32700            2      1.06383    100_vs_Mock_1h     treatment*transcript       5      24      16.70    <.0001    0.001030      1
    AT2G32700            2      1.06383    10_vs_Mock_1h      treatment*transcript       5      24       8.70    <.0001    0.033655      1
    AT2G43070            2      1.06383    10_vs_Mock_24h     treatment*transcript       2      12     191.95    <.0001    0.000006      1
    AT2G43070            2      1.06383    10_vs_Mock_72h     treatment*transcript       2      12     258.30    <.0001    0.000001      1
    AT2G43700            2      1.06383    100_vs_Mock_24h    treatment*transcript       3      16      15.21    <.0001    0.048229      1
    AT2G43700            2      1.06383    10_vs_Mock_24h     treatment*transcript       3      16      15.28    <.0001    0.028629      1
    AT4G16990            2      1.06383    100_vs_Mock_1h     treatment*transcript      10      44       6.56    <.0001    0.005145      1
    AT4G16990            2      1.06383    10_vs_Mock_24h     treatment*transcript      10      44       7.34    <.0001    0.001569      1
    AT5G14740            2      1.06383    100_vs_Mock_3h     treatment*transcript       8      36       7.38    <.0001    0.015839      1
    AT5G14740            2      1.06383    10_vs_Mock_3h      treatment*transcript       8      36      10.61    <.0001    0.000371      1
    AT5G24160            2      1.06383    100_vs_Mock_1h     treatment*transcript       2      12      19.58    0.0002    0.048242      1
    AT5G24160            2      1.06383    10_vs_Mock_1h      treatment*transcript       2      12      23.38    <.0001    0.033655      1
    AT5G38460            2      1.06383    100_vs_Mock_72h    treatment*transcript       1       8    1990.18    <.0001    0.000001      1
    AT5G38460            2      1.06383    10_vs_Mock_72h     treatment*transcript       1       8      91.50    <.0001    0.014906      1
    AT5G45060            2      1.06383    100_vs_Mock_1h     treatment*transcript       1       8      72.61    <.0001    0.016172      1
    AT5G45060            2      1.06383    10_vs_Mock_1h      treatment*transcript       1       8      50.73    <.0001    0.033655      1
    AT5G48655            3      1.59574    100_vs_Mock_1h     treatment*transcript       4      20      10.08    0.0001    0.037397      1
    AT5G48655            3      1.59574    100_vs_Mock_72h    treatment*transcript       4      20      20.17    <.0001    0.003781      1
    AT5G48655            3      1.59574    10_vs_Mock_3h      treatment*transcript       4      20      24.72    <.0001    0.000371      1
    AT5G57180            2      1.06383    100_vs_Mock_24h    treatment*transcript       3      16      17.33    <.0001    0.031702      1
    AT5G57180            2      1.06383    10_vs_Mock_24h     treatment*transcript       3      16      27.26    <.0001    0.001986      1
    AT5G57940            2      1.06383    100_vs_Mock_1h     treatment*transcript       1       8      43.48    0.0002    0.048242      1
    AT5G57940            2      1.06383    10_vs_Mock_1h      treatment*transcript       1       8     383.23    <.0001    0.000141      1
    AT5G63370            2      1.06383    100_vs_Mock_1h     treatment*transcript       7      32      18.94    <.0001    0.000009      1
    AT5G63370            2      1.06383    10_vs_Mock_1h      treatment*transcript       7      32      21.34    <.0001    0.000002      1


