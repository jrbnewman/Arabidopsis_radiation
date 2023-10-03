ods listing; ods html close;
libname arabRNA "!HOME/concannon/DTRA/arabidopsis/sas_data";
libname tair "!HOME/concannon/useful_arabidopsis_data/tair10/sas_data";

/* Design file to iterate through */

data design;
  set arabRNA.arab_design_file;
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


%include '!HOME/concannon/DTRA/arabidopsis/sas_programs/iterdataset.sas';

/* Import RSEM counts */

%macro importRSEM(sample,sampnum);

proc import datafile="!HOME/concannon/DTRA/arabidopsis/TEcount_output_uniq/&sample..cntTable"
     out=uniq dbms=tab replace; getnames=no; datarow=2;
     guessingrows=max;
run;

proc import datafile="!HOME/concannon/DTRA/arabidopsis/TEcount_output_multi/&sample..cntTable"
     out=multi dbms=tab replace; getnames=no; datarow=2;
     guessingrows=max;
run;


proc import datafile="!HOME/concannon/DTRA/arabidopsis/TEcount_output_uniq2/&sample..cntTable"
     out=uniq2 dbms=tab replace; getnames=no; datarow=2;
     guessingrows=max;
run;

proc import datafile="!HOME/concannon/DTRA/arabidopsis/TEcount_output_multi2/&sample..cntTable"
     out=multi2 dbms=tab replace; getnames=no; datarow=2;
     guessingrows=max;
run;



proc import datafile="!HOME/concannon/DTRA/arabidopsis/TElocal_output_uniq/&sample..cntTable"
     out=uniq_local dbms=tab replace; getnames=no; datarow=2;
     guessingrows=max;
run;

proc import datafile="!HOME/concannon/DTRA/arabidopsis/TElocal_output_multi/&sample..cntTable"
     out=multi_local dbms=tab replace; getnames=no; datarow=2;
     guessingrows=max;
run;




data WORK.BWA_UNIQ_&sampnum.    ;
	%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
	infile "!HOME/concannon/DTRA/arabidopsis/TE_elements_uniq/cvrg_cnts_&sample._uniq.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
	informat sample_id $11. ;
	informat element_id $23. ;
	informat mapped_reads best32. ;
	informat read_length best32. ;
	informat region_length best32. ;
	informat region_depth best32. ;
	informat reads_in_region best32. ;
	informat apn best32. ;
	informat rpkm best32. ;
	informat mean best32. ;
	informat std best32. ;
	informat cv best32. ;
	format sample_id $11. ;
	format element_id $23. ;
	format mapped_reads best12. ;
	format read_length best12. ;
	format region_length best12. ;
	format region_depth best12. ;
	format reads_in_region best12. ;
	format apn best12. ;
	format rpkm best12. ;
	format mean best12. ;
	format std best12. ;
	format cv best12. ;
	input
	sample_id  $
	element_id  $
	mapped_reads
	read_length
	region_length
	region_depth
	reads_in_region
	apn
	rpkm
	mean
	std
	cv
	;
	if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;

data WORK.BWA_ALL_&sampnum.    ;
	%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
	infile "!HOME/concannon/DTRA/arabidopsis/TE_elements_all/cvrg_cnts_&sample._all.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
	informat sample_id $11. ;
	informat element_id $23. ;
	informat mapped_reads best32. ;
	informat read_length best32. ;
	informat region_length best32. ;
	informat region_depth best32. ;
	informat reads_in_region best32. ;
	informat apn best32. ;
	informat rpkm best32. ;
	informat mean best32. ;
	informat std best32. ;
	informat cv best32. ;
	format sample_id $11. ;
	format element_id $23. ;
	format mapped_reads best12. ;
	format read_length best12. ;
	format region_length best12. ;
	format region_depth best12. ;
	format reads_in_region best12. ;
	format apn best12. ;
	format rpkm best12. ;
	format mean best12. ;
	format std best12. ;
	format cv best12. ;
	input
	sample_id  $
	element_id  $
	mapped_reads
	read_length
	region_length
	region_depth
	reads_in_region
	apn
	rpkm
	mean
	std
	cv
	;
	if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;



data te_uniq_&sampnum.;
   length sample_id $50.;
   set uniq;
   sample_id="&sample.";
   rename VAR1=gene_id VAR2=count;
run;

data te_multi_&sampnum.;
   length sample_id $50.;
   set multi;
   sample_id="&sample.";
   rename VAR1=gene_id VAR2=count;
run;

data te_uniq2_&sampnum.;
   length sample_id $50.;
   set uniq2;
   sample_id="&sample.";
   rename VAR1=gene_id VAR2=count;
run;

data te_multi2_&sampnum.;
   length sample_id $50.;
   set multi2;
   sample_id="&sample.";
   rename VAR1=gene_id VAR2=count;
run;


data te_uniqlocal_&sampnum.;
   length sample_id $50.;
   set uniq_local;
   sample_id="&sample.";
   rename VAR1=gene_id VAR2=count;
run;

data te_multilocal_&sampnum.;
   length sample_id $50.;
   set multi_local;
   sample_id="&sample.";
   rename VAR1=gene_id VAR2=count;
run;

%mend;
%iterdataset(dataset=design, function=%nrstr(%importRSEM(&sample_id, &sample_number)));

data uniq_all_counts;
   set te_uniq_: ;
run;

data multi_all_counts;
   set te_multi_: ;
run;

data uniqlocal_all_counts;
   set te_uniqlocal_: ;
run;

data multilocal_all_counts;
   set te_multilocal_: ;
run;

data uniq2_all_counts;
   set te_uniq2_: ;
run;

data multi2_all_counts;
   set te_multi2_: ;
run;

data uniq_bwa_all_counts;
   set bwa_uniq_: ;
run;

data multi_bwa_all_counts;
   set bwa_all_: ;
run;


proc freq data=uniq_all_counts noprint;
  tables gene_id / out=gene_count_check;
run;

data gene_count_check2;
  set gene_count_check;
  where count ne 36;
run; *looks good, expect all genes to have 36 observations;

proc freq data=multi_all_counts noprint;
  tables gene_id / out=gene_count_check;
run;

data gene_count_check2;
  set gene_count_check;
  where count ne 36;
run; *looks good, expect all genes to have 36 observations;

proc sort data=design;
  by sample_id;
proc sort data=uniq_all_counts;
  by sample_id;
proc sort data=multi_all_counts;
  by sample_id;
proc sort data=uniqlocal_all_counts;
  by sample_id;
proc sort data=multilocal_all_counts;
  by sample_id;
proc sort data=uniq2_all_counts;
  by sample_id;
proc sort data=multi2_all_counts;
  by sample_id;
proc sort data=uniq_bwa_all_counts;
  by sample_id;
proc sort data=multi_bwa_all_counts;
  by sample_id;
run;

data uniq_counts_w_key;
  merge design (in=in1) uniq_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;
data multi_counts_w_key;
  merge design (in=in1) multi_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;

data uniqlocal_counts_w_key;
  merge design (in=in1) uniqlocal_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;
data multilocal_counts_w_key;
  merge design (in=in1) multilocal_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;


data uniq2_counts_w_key;
  merge design (in=in1) uniq2_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;
data multi2_counts_w_key;
  merge design (in=in1) multi2_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;


data uniq_bwa_counts_w_key;
  merge design (in=in1) uniq_bwa_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;
data multi_bwa_counts_w_key;
  merge design (in=in1) multi_bwa_all_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;


data genes;
   set tair.tair20_exon2gene;
   keep gene_id;
run;

proc sort data=genes nodup;
  by gene_id;
proc sort data=uniq_counts_w_key;
  by gene_id;
proc sort data=multi_counts_w_key;
  by gene_id;
proc sort data=uniqlocal_counts_w_key;
  by gene_id;
proc sort data=multilocal_counts_w_key;
  by gene_id;
proc sort data=uniq2_counts_w_key;
  by gene_id;
proc sort data=multi2_counts_w_key;
  by gene_id;
run;

data uniq_counts_w_key2;
  merge uniq_counts_w_key (in=in1) genes (in=in2);
  by gene_id;
  if in2 then flag_TE=0;
  else flag_TE=1;
run;

data multi_counts_w_key2;
  merge multi_counts_w_key (in=in1) genes (in=in2);
  by gene_id;
  if in2 then flag_TE=0;
  else flag_TE=1;
run;
data uniqlocal_counts_w_key2;
  merge uniqlocal_counts_w_key (in=in1) genes (in=in2);
  by gene_id;
  if in2 then flag_TE=0;
  else flag_TE=1;
run;

data multilocal_counts_w_key2;
  merge multilocal_counts_w_key (in=in1) genes (in=in2);
  by gene_id;
  if in2 then flag_TE=0;
  else flag_TE=1;
run;



data uniq2_counts_w_key2;
  merge uniq2_counts_w_key (in=in1) genes (in=in2);
  by gene_id;
  if in2 then flag_TE=0;
  else flag_TE=1;
run;

data multi2_counts_w_key2;
  merge multi2_counts_w_key (in=in1) genes (in=in2);
  by gene_id;
  if in2 then flag_TE=0;
  else flag_TE=1;
run;



/* Make permenant */

data arabRNA.counts_by_te_uniq;
  set uniq_counts_w_key2;
run;

data arabRNA.counts_by_te_multi;
  set multi_counts_w_key2;
run;

data arabRNA.counts_by_te_xscript_uniq;
  set uniq2_counts_w_key2;
run;

data arabRNA.counts_by_te_xscript_multi;
  set multi2_counts_w_key2;
run;

data arabRNA.counts_by_te_local_uniq;
  set uniqlocal_counts_w_key2;
run;

data arabRNA.counts_by_te_local_multi;
  set multilocal_counts_w_key2;
run;

data arabRNA.counts_by_te_bwa_uniq;
  set uniq_bwa_counts_w_key;
run;

data arabRNA.counts_by_te_bwa_multi;
  set multi_bwa_counts_w_key;
run;


