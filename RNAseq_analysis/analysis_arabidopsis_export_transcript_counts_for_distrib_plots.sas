ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Prep data for exporting so I can plot the distribution of coverage. We can use this to decide if
   there are low-coverage samples that need dropping

   I want to include: sampleID, fusionID, log_apn, plate, condition (trt*time) */

data counts_w_key;
   set rs.counts_by_isoform ;
run;

data calc_log;
  set counts_w_key;
  length condition $50.;
  log_tpm=log(tpm+1);
  condition=catx("_", treatment,time);
  keep sample_id transcript_id treatment time condition tpm log_tpm;
run;

proc sort data=calc_log;
  by sample_id transcript_id;
run;


data drop_zero;
  set calc_log;
  where tpm > 0;
run;

/* Export for plotting coverage distribution */

proc export data=calc_log
     outfile="!PATCON/arabidopsis/analysis_output/distrib_logtpm_transcripts_arabidopsis.csv"
     dbms=csv replace;
run;

proc export data=drop_zero
     outfile="!PATCON/arabidopsis/analysis_output/distrib_logapn_transcripts_arabidopsis_gt0.csv"
     dbms=csv replace;
run;

