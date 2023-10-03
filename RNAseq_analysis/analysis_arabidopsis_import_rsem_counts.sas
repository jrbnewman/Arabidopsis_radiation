ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";

/* Design file to iterate through */

data design;
  set rs.arab_design_file;
  length sample_id2 $50.;
  if grays=0 then sample_id2=catx("-",sample_number,"M",time,replicate);
  else if grays=0.1 then sample_id2=catx("-",sample_number,"0-1",time,replicate);
  else sample_id2=catx("-",sample_number,"1",time,replicate);
  keep sample_id2 sample_number treatment time replicate;
  rename sample_id2=sample_id;
run;

proc sort data=design nodup;
  by sample_id ;
run;


%include '!PATCON/arabidopsis/sas_programs/iterdataset.sas';

/* Import RSEM counts */

%macro importRSEM(sample,sampnum);

proc import datafile="!PATCON/arabidopsis/alignment_output/rsem_output/&sample._expr.isoforms.results"
     out=rsem_counts dbms=tab replace;
     guessingrows=55160;
run;

data rsem_cnts_&sampnum.;
   length sample_id $50.;
   set rsem_counts;
   sample_id="&sample.";
run;

%mend;
%iterdataset(dataset=design, function=%nrstr(%importRSEM(&sample_id, &sample_number)));

data rsem_all_counts;
   set rsem_cnts_: ;
run;

proc freq data=rsem_all_counts noprint;
  tables transcript_id / out=xs_count_check;
run;

data xs_count_check2;
  set xs_count_check;
  where count ne 36;
run; *looks good, expect all transcripts to have 36 observations;

proc sort data=design;
  by sample_id;
proc sort data=rsem_all_counts;
  by sample_id;
run;

data rsem_counts_w_key;
  merge design (in=in1) rsem_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;


/* Make permenant */

data rs.counts_by_isoform;
  set rsem_counts_w_key;
run;


