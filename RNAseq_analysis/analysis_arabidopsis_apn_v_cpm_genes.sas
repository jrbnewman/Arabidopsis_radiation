/* Check concordance of results between gene-level analyses */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

data old_counts;
  set rs.arab_gene_mean_apn_by_trt_time;
run;

data new_counts;
  set rs.arab_gene_mean_cpm_by_trt_time;
run;

proc sort data=old_counts;
  by gene_id;
proc sort data=new_counts;
  by gene_id;
run;

data old_v_new_counts;
  merge new_counts (in=in1) old_counts (in=in2);
  by gene_id;
  if in1 and in2;
run;

%macro calcCorr(var1,var2);
data old_v_new_counts2;
  set old_v_new_counts;
  log_&var1.=log(&var1. + 1);
  log_&var2.=log(&var2. + 1);
run;

proc corr data=old_v_new_counts2 pearson;
  var log_&var1. log_&var2.;
run;

%mend;

%calcCorr(mean_q3apn_mock_1h,mean_cpm_mock_1h);
%calcCorr(mean_q3apn_mock_3h,mean_cpm_mock_3h);
%calcCorr(mean_q3apn_mock_24h,mean_cpm_mock_24h);
%calcCorr(mean_q3apn_mock_72h,mean_cpm_mock_72h);

%calcCorr(mean_q3apn_01gy_1h,mean_cpm_01gy_1h);
%calcCorr(mean_q3apn_01gy_3h,mean_cpm_01gy_3h);
%calcCorr(mean_q3apn_01gy_24h,mean_cpm_01gy_24h);
%calcCorr(mean_q3apn_01gy_72h,mean_cpm_01gy_72h);

%calcCorr(mean_q3apn_1gy_1h,mean_cpm_1gy_1h);
%calcCorr(mean_q3apn_1gy_3h,mean_cpm_1gy_3h);
%calcCorr(mean_q3apn_1gy_24h,mean_cpm_1gy_24h);
%calcCorr(mean_q3apn_1gy_72h,mean_cpm_1gy_72h);

/*
			Rvalue	Pvalue
MOCK 1h		0.95743	<.0001
MOCK 3h		0.95673	<.0001
MOCK 24h	0.95557	<.0001
MOCK 72h	0.95634	<.0001
0.1gy 1h	0.95749	<.0001
0.1gy 3h	0.95702	<.0001
0.1gy 24h	0.95490	<.0001
0.1gy 72h	0.95652	<.0001
1gy 1h		0.95717	<.0001
1gy 3h		0.95652	<.0001
1gy 24h		0.95544	<.0001
1gy 72h		0.95610	<.0001
*/

proc export data=old_v_new_counts outfile="!PATCON/arabidopsis/analysis_output/gene_counts_q3apn_v_cpm.csv"
  dbms=csv replace;
run;


