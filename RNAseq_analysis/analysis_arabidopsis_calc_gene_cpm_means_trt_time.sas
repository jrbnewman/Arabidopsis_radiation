ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";

/* Calculating the means for treatment*time. */

/* Get coverage counts. Going to calculate means on the non-log, UQ-normalized APN */

data coverage;
   set rs.cpm_norm_counts_by_gene;
   keep sample_id gene_id transcript_id_s_ cpm;
run;

/* Get timepoints and treatment for samples */
   
data sample_key;
  set rs.arab_design_file;
   length sample_id2 $50.;
   if grays=0 then sample_id2=catx("-",sample_number,"M",time,replicate);
   else if grays=0.1 then sample_id2=catx("-",sample_number,"0-1",time,replicate);
   else if grays=1 then sample_id2=catx("-",sample_number,"1",time,replicate);
  keep sample_id2 treatment time;
  rename sample_id2=sample_id;
run;

proc sort data=sample_key nodup;
   by sample_id;
proc sort data=coverage;
   by sample_id;
run;

data gene_counts_w_key;
   merge sample_key (in=in1) coverage (in=in2);
   by sample_id;
   if in1 and in2;
   run;

/* Calc means for treatment*time by fusion */

proc sort data=gene_counts_w_key;
   by gene_id transcript_id_s_ treatment time;
proc means data=gene_counts_w_key noprint;
   by gene_id transcript_id_s_ treatment time;
   var cpm;
   output out=mean_cpm_by_trt_time mean=;
run;


/* Transpose */

data mean_cpm_by_trt_time_2;
  set mean_cpm_by_trt_time;
  trt_time=catx('_',treatment,time);
run;

proc sort data=mean_cpm_by_trt_time_2;
  by gene_id transcript_id_s_ trt_time;
proc transpose data=mean_cpm_by_trt_time_2 out=mean_cpm_sbys;
   by gene_id transcript_id_s_;
   var cpm;
   id trt_time;
run;

/* Make permenant */

data rs.arab_gene_mean_cpm_by_trt_time;
   set mean_cpm_sbys;
   drop _NAME_;
   rename _0_1gy_1=mean_cpm_01gy_1h
          _0_1gy_3=mean_cpm_01gy_3h
          _0_1gy_24=mean_cpm_01gy_24h
          _0_1gy_72=mean_cpm_01gy_72h
          _1gy_1=mean_cpm_1gy_1h
          _1gy_3=mean_cpm_1gy_3h
          _1gy_24=mean_cpm_1gy_24h
          _1gy_72=mean_cpm_1gy_72h
          Mock_1=mean_cpm_mock_1h
          Mock_3=mean_cpm_mock_3h
          Mock_24=mean_cpm_mock_24h
          Mock_72=mean_cpm_mock_72h ;
run;


