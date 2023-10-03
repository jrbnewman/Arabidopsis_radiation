ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";

/* Normalize gene read counts using "copies per million" CPM: reads aligning to gene / 1,000,000
   [reads aligning to gene] * 10^6 / [total number of aligned reads in one sample] */

data counts_by_gene;
  set rs.counts_by_gene;
run;

/* sum expected_count by sample */

proc sort data=counts_by_gene;
  by sample_id;
proc means data=counts_by_gene noprint;
  by sample_id;
  var expected_count;
  output out=total_counts_per_sample sum=total_counts;
run;

proc sort data=total_counts_per_sample;
  by sample_id;
proc sort data=counts_by_gene;
  by sample_id;
run;

data counts_by_gene2;
  merge counts_by_gene (in=in1) total_counts_per_sample (in=in2);
  by sample_id;
  if in1 and in2;
run;

/* Calculate CPM from expected_counts */

data calc_cpm;
  set counts_by_gene2;
  cpm = (expected_count * 1000000)/total_counts;
run;

data rs.cpm_norm_counts_by_gene;
  set calc_cpm;
  drop _TYPE_ _FREQ_;
run;


