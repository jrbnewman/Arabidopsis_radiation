/* Rep-to-rep concordance -- do for mean log_q3_apn */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";


/* Generate gene expression matrix */

data apn_tall;
  set rs.arab_mean_q3_apn_by_gene;
  log_q3_q3_apn=log(q3_q3_apn+1);
  keep sample_number gene_id log_q3_q3_apn;
run;

data design;
   set rs.arab_design_file;
   length sample_name $20.;
   sample_name=catx("_","Sample",treatment,time,replicate);
   keep sample_number sample_name;
run;

proc sort data=apn_tall;
  by sample_number;
proc sort data=design nodup;
   by sample_number;
run;

data apn_w_key;
  merge design (in=in1) apn_tall (in=in2);
   by sample_number;
  if in1 and in2;
run;

proc sort data=apn_w_key;
   by gene_id sample_name;
proc transpose data=apn_w_key out=apn_sbys;
   by gene_id;
   id sample_name;
   var log_q3_q3_apn;
run;


/* Calculate means and differences for BA plots */

data apn_sbys_ba;
   set apn_sbys;
   mean_01gy_1h_12=(Sample_0_1gy_1_1+Sample_0_1gy_1_2)/2;
   diff_01gy_1h_12=(Sample_0_1gy_1_1-Sample_0_1gy_1_2);
   mean_01gy_1h_13=(Sample_0_1gy_1_1+Sample_0_1gy_1_3)/2;
   diff_01gy_1h_13=(Sample_0_1gy_1_1-Sample_0_1gy_1_3);
   mean_01gy_1h_23=(Sample_0_1gy_1_2+Sample_0_1gy_1_3)/2;
   diff_01gy_1h_23=(Sample_0_1gy_1_2-Sample_0_1gy_1_3);

   mean_01gy_3h_12=(Sample_0_1gy_3_1+Sample_0_1gy_3_2)/2;
   diff_01gy_3h_12=(Sample_0_1gy_3_1-Sample_0_1gy_3_2);
   mean_01gy_3h_13=(Sample_0_1gy_3_1+Sample_0_1gy_3_3)/2;
   diff_01gy_3h_13=(Sample_0_1gy_3_1-Sample_0_1gy_3_3);
   mean_01gy_3h_23=(Sample_0_1gy_3_2+Sample_0_1gy_3_3)/2;
   diff_01gy_3h_23=(Sample_0_1gy_3_2-Sample_0_1gy_3_3);

   mean_01gy_24h_12=(Sample_0_1gy_24_1+Sample_0_1gy_24_2)/2;
   diff_01gy_24h_12=(Sample_0_1gy_24_1-Sample_0_1gy_24_2);
   mean_01gy_24h_13=(Sample_0_1gy_24_1+Sample_0_1gy_24_3)/2;
   diff_01gy_24h_13=(Sample_0_1gy_24_1-Sample_0_1gy_24_3);
   mean_01gy_24h_23=(Sample_0_1gy_24_2+Sample_0_1gy_24_3)/2;
   diff_01gy_24h_23=(Sample_0_1gy_24_2-Sample_0_1gy_24_3);

   mean_01gy_72h_12=(Sample_0_1gy_72_1+Sample_0_1gy_72_2)/2;
   diff_01gy_72h_12=(Sample_0_1gy_72_1-Sample_0_1gy_72_2);
   mean_01gy_72h_13=(Sample_0_1gy_72_1+Sample_0_1gy_72_3)/2;
   diff_01gy_72h_13=(Sample_0_1gy_72_1-Sample_0_1gy_72_3);
   mean_01gy_72h_23=(Sample_0_1gy_72_2+Sample_0_1gy_72_3)/2;
   diff_01gy_72h_23=(Sample_0_1gy_72_2-Sample_0_1gy_72_3);

   drop _NAME_;
run;

data genes_on;
   set rs.arab_flag_gene_on_gt0;
   where flag_gene_on_01gy_apn0=1 and flag_gene_on_1gy_apn0=1 and flag_gene_on_mock_apn0=1;
   keep gene_id;
run;

proc sort data=apn_sbys;
  by gene_id;
proc sort data=genes_on;
  by gene_id;
run;

data apn_sbys2;
  merge genes_on (in=in1) apn_sbys (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Rep-to-rep concordance -- Pearson */

ods graphics on;
ods pdf file = "!PATCON/arabidopsis/rep_to_rep_concordance_and_ba_plots.pdf";

%macro plotQC(sample1,sample2,mean,diff);

data ba_plot;
  set apn_sbys2;
  &mean.=(&sample1. + &sample2.) /2 ;
  &diff.=&sample1. - &sample2. ;
run;

title "Concordance between &sample1. and &sample2.";
proc corr data=apn_sbys2 nomiss;
   var &sample1. &sample2.;
run;

title "Concordance between &sample1. and &sample2.";
proc sgplot data=ba_plot;
   scatter x=&sample1. y=&sample2.;
run;

title "Bland-Altmans plot for &sample1. and &sample2.";
proc sgplot data=ba_plot;
   scatter x=&mean. y=&diff.;
run;

%mend;


%plotQC(Sample_0_1gy_1_1,Sample_0_1gy_1_2,mean_01gy_1h_rep1_2,diff_01gy_1h_rep1_2);
%plotQC(Sample_0_1gy_1_1,Sample_0_1gy_1_3,mean_01gy_1h_rep1_3,diff_01gy_1h_rep1_3);
%plotQC(Sample_0_1gy_1_2,Sample_0_1gy_1_3,mean_01gy_1h_rep2_3,diff_01gy_1h_rep2_3);
%plotQC(Sample_0_1gy_3_1,Sample_0_1gy_3_2,mean_01gy_3h_rep1_2,diff_01gy_3h_rep1_2);
%plotQC(Sample_0_1gy_3_1,Sample_0_1gy_3_3,mean_01gy_3h_rep1_3,diff_01gy_3h_rep1_3);
%plotQC(Sample_0_1gy_3_2,Sample_0_1gy_3_3,mean_01gy_3h_rep2_3,diff_01gy_3h_rep2_3);
%plotQC(Sample_0_1gy_24_1,Sample_0_1gy_24_2,mean_01gy_24h_rep1_2,diff_01gy_24h_rep1_2);
%plotQC(Sample_0_1gy_24_1,Sample_0_1gy_24_3,mean_01gy_24h_rep1_3,diff_01gy_24h_rep1_3);
%plotQC(Sample_0_1gy_24_2,Sample_0_1gy_24_3,mean_01gy_24h_rep2_3,diff_01gy_24h_rep2_3);
%plotQC(Sample_0_1gy_72_1,Sample_0_1gy_72_2,mean_01gy_72h_rep1_2,diff_01gy_72h_rep1_2);
%plotQC(Sample_0_1gy_72_1,Sample_0_1gy_72_3,mean_01gy_72h_rep1_3,diff_01gy_72h_rep1_3);
%plotQC(Sample_0_1gy_72_2,Sample_0_1gy_72_3,mean_01gy_72h_rep2_3,diff_01gy_72h_rep2_3);

%plotQC(Sample_1gy_1_1,Sample_1gy_1_2,mean_1gy_1h_rep1_2,diff_1gy_1h_rep1_2);
%plotQC(Sample_1gy_1_1,Sample_1gy_1_3,mean_1gy_1h_rep1_3,diff_1gy_1h_rep1_3);
%plotQC(Sample_1gy_1_2,Sample_1gy_1_3,mean_1gy_1h_rep2_3,diff_1gy_1h_rep2_3);
%plotQC(Sample_1gy_3_1,Sample_1gy_3_2,mean_1gy_3h_rep1_2,diff_1gy_3h_rep1_2);
%plotQC(Sample_1gy_3_1,Sample_1gy_3_3,mean_1gy_3h_rep1_3,diff_1gy_3h_rep1_3);
%plotQC(Sample_1gy_3_2,Sample_1gy_3_3,mean_1gy_3h_rep2_3,diff_1gy_3h_rep2_3);
%plotQC(Sample_1gy_24_1,Sample_1gy_24_2,mean_1gy_24h_rep1_2,diff_1gy_24h_rep1_2);
%plotQC(Sample_1gy_24_1,Sample_1gy_24_3,mean_1gy_24h_rep1_3,diff_1gy_24h_rep1_3);
%plotQC(Sample_1gy_24_2,Sample_1gy_24_3,mean_1gy_24h_rep2_3,diff_1gy_24h_rep2_3);
%plotQC(Sample_1gy_72_1,Sample_1gy_72_2,mean_1gy_72h_rep1_2,diff_1gy_72h_rep1_2);
%plotQC(Sample_1gy_72_1,Sample_1gy_72_3,mean_1gy_72h_rep1_3,diff_1gy_72h_rep1_3);
%plotQC(Sample_1gy_72_2,Sample_1gy_72_3,mean_1gy_72h_rep2_3,diff_1gy_72h_rep2_3);

%plotQC(Sample_Mock_1_1,Sample_Mock_1_2,mean_Mock_1h_rep1_2,diff_Mock_1h_rep1_2);
%plotQC(Sample_Mock_1_1,Sample_Mock_1_3,mean_Mock_1h_rep1_3,diff_Mock_1h_rep1_3);
%plotQC(Sample_Mock_1_2,Sample_Mock_1_3,mean_Mock_1h_rep2_3,diff_Mock_1h_rep2_3);
%plotQC(Sample_Mock_3_1,Sample_Mock_3_2,mean_Mock_3h_rep1_2,diff_Mock_3h_rep1_2);
%plotQC(Sample_Mock_3_1,Sample_Mock_3_3,mean_Mock_3h_rep1_3,diff_Mock_3h_rep1_3);
%plotQC(Sample_Mock_3_2,Sample_Mock_3_3,mean_Mock_3h_rep2_3,diff_Mock_3h_rep2_3);
%plotQC(Sample_Mock_24_1,Sample_Mock_24_2,mean_Mock_24h_rep1_2,diff_Mock_24h_rep1_2);
%plotQC(Sample_Mock_24_1,Sample_Mock_24_3,mean_Mock_24h_rep1_3,diff_Mock_24h_rep1_3);
%plotQC(Sample_Mock_24_2,Sample_Mock_24_3,mean_Mock_24h_rep2_3,diff_Mock_24h_rep2_3);
%plotQC(Sample_Mock_72_1,Sample_Mock_72_2,mean_Mock_72h_rep1_2,diff_Mock_72h_rep1_2);
%plotQC(Sample_Mock_72_1,Sample_Mock_72_3,mean_Mock_72h_rep1_3,diff_Mock_72h_rep1_3);
%plotQC(Sample_Mock_72_2,Sample_Mock_72_3,mean_Mock_72h_rep2_3,diff_Mock_72h_rep2_3);


ods pdf close ;
ods graphics off;


/* Export data for JMP */

data apn_sbys3;
   set apn_sbys2;
   drop _NAME_;
run;


data design2;
   length columnname $20.;
   length Array $20.;
   length sampleID $20.;
   length condition $20.;
   set rs.arab_design_file;
   columnname=catx("_","Sample",treatment,time,replicate);
   Array=catx("_","Sample",treatment,time,replicate);
   sampleID=catx("_","Sample",treatment,time,replicate);
   condition=catx("_",treatment,time);
   keep columnname array sampleID condition treatment time;
run;

proc sort data=apn_sbys3;
  by gene_id;
proc sort data=design2 nodup;
  by sampleID;
run;

proc export data=apn_sbys3 outfile="!PATCON/arabidopsis/analysis_output/arabidopsis_jmp_all_exp_genes_data.csv"
    dbms=tab replace;
run;

proc export data=design2 outfile="!PATCON/arabidopsis/analysis_output/arabidopsis_design.csv"
    dbms=tab replace;
run;


/* Subset genes with at least one contrast at nominal P<0.05 */


data genes2keep;
   set rs.fdr_by_gene_v2;
   where (p_01gy_v_Mock_1h < 0.05 and p_01gy_v_Mock_1h ne .) 
      or (p_01gy_v_Mock_3h < 0.05 and p_01gy_v_Mock_3h ne .) 
      or (p_01gy_v_Mock_24h < 0.05 and p_01gy_v_Mock_24h ne .) 
      or (p_01gy_v_Mock_72h < 0.05 and p_01gy_v_Mock_72h ne .) 
      or (p_1gy_v_Mock_1h < 0.05 and p_1gy_v_Mock_1h ne .) 
      or (p_1gy_v_Mock_3h < 0.05 and p_1gy_v_Mock_3h ne .) 
      or (p_1gy_v_Mock_24h < 0.05 and p_1gy_v_Mock_24h ne .) 
      or (p_1gy_v_Mock_72h < 0.05 and p_1gy_v_Mock_72h ne .) 
      ;
   keep gene_id;
run;

proc sort data=apn_sbys3;
   by gene_id;
proc sort data=genes2keep;
  by gene_id;
run;

data apn_sbys_p05;
   merge genes2keep (in=in1) apn_sbys2 (in=in2);
   by gene_id;
   if in1 and in2;
run;


/* Rep-to-rep concordance -- Pearson */

ods graphics on;
ods pdf file = "!PATCON/arabidopsis/rep_to_rep_concordance_and_ba_plots_p05.pdf";

%macro plotQCsub(sample1,sample2,mean,diff);

data ba_plot;
  set apn_sbys_p05;
  &mean.=(&sample1. + &sample2.) /2 ;
  &diff.=&sample1. - &sample2. ;
run;


title "Concordance between &sample1. and &sample2.";
proc corr data=apn_sbys2 nomiss;
   var &sample1. &sample2.;
run;

title "Concordance between &sample1. and &sample2.";
proc sgplot data=ba_plot;
   scatter x=&sample1. y=&sample2.;
run;

title "Bland-Altman plot for &sample1. and &sample2.";
proc sgplot data=ba_plot;
   scatter x=&mean. y=&diff.;
run;

%mend;

%plotQCsub(Sample_0_1gy_1_1,Sample_0_1gy_1_2,mean_01gy_1h_rep1_2,diff_01gy_1h_rep1_2);
%plotQCsub(Sample_0_1gy_1_1,Sample_0_1gy_1_3,mean_01gy_1h_rep1_3,diff_01gy_1h_rep1_3);
%plotQCsub(Sample_0_1gy_1_2,Sample_0_1gy_1_3,mean_01gy_1h_rep2_3,diff_01gy_1h_rep2_3);
%plotQCsub(Sample_0_1gy_3_1,Sample_0_1gy_3_2,mean_01gy_3h_rep1_2,diff_01gy_3h_rep1_2);
%plotQCsub(Sample_0_1gy_3_1,Sample_0_1gy_3_3,mean_01gy_3h_rep1_3,diff_01gy_3h_rep1_3);
%plotQCsub(Sample_0_1gy_3_2,Sample_0_1gy_3_3,mean_01gy_3h_rep2_3,diff_01gy_3h_rep2_3);
%plotQCsub(Sample_0_1gy_24_1,Sample_0_1gy_24_2,mean_01gy_24h_rep1_2,diff_01gy_24h_rep1_2);
%plotQCsub(Sample_0_1gy_24_1,Sample_0_1gy_24_3,mean_01gy_24h_rep1_3,diff_01gy_24h_rep1_3);
%plotQCsub(Sample_0_1gy_24_2,Sample_0_1gy_24_3,mean_01gy_24h_rep2_3,diff_01gy_24h_rep2_3);
%plotQCsub(Sample_0_1gy_72_1,Sample_0_1gy_72_2,mean_01gy_72h_rep1_2,diff_01gy_72h_rep1_2);
%plotQCsub(Sample_0_1gy_72_1,Sample_0_1gy_72_3,mean_01gy_72h_rep1_3,diff_01gy_72h_rep1_3);
%plotQCsub(Sample_0_1gy_72_2,Sample_0_1gy_72_3,mean_01gy_72h_rep2_3,diff_01gy_72h_rep2_3);

%plotQCsub(Sample_1gy_1_1,Sample_1gy_1_2,mean_1gy_1h_rep1_2,diff_1gy_1h_rep1_2);
%plotQCsub(Sample_1gy_1_1,Sample_1gy_1_3,mean_1gy_1h_rep1_3,diff_1gy_1h_rep1_3);
%plotQCsub(Sample_1gy_1_2,Sample_1gy_1_3,mean_1gy_1h_rep2_3,diff_1gy_1h_rep2_3);
%plotQCsub(Sample_1gy_3_1,Sample_1gy_3_2,mean_1gy_3h_rep1_2,diff_1gy_3h_rep1_2);
%plotQCsub(Sample_1gy_3_1,Sample_1gy_3_3,mean_1gy_3h_rep1_3,diff_1gy_3h_rep1_3);
%plotQCsub(Sample_1gy_3_2,Sample_1gy_3_3,mean_1gy_3h_rep2_3,diff_1gy_3h_rep2_3);
%plotQCsub(Sample_1gy_24_1,Sample_1gy_24_2,mean_1gy_24h_rep1_2,diff_1gy_24h_rep1_2);
%plotQCsub(Sample_1gy_24_1,Sample_1gy_24_3,mean_1gy_24h_rep1_3,diff_1gy_24h_rep1_3);
%plotQCsub(Sample_1gy_24_2,Sample_1gy_24_3,mean_1gy_24h_rep2_3,diff_1gy_24h_rep2_3);
%plotQCsub(Sample_1gy_72_1,Sample_1gy_72_2,mean_1gy_72h_rep1_2,diff_1gy_72h_rep1_2);
%plotQCsub(Sample_1gy_72_1,Sample_1gy_72_3,mean_1gy_72h_rep1_3,diff_1gy_72h_rep1_3);
%plotQCsub(Sample_1gy_72_2,Sample_1gy_72_3,mean_1gy_72h_rep2_3,diff_1gy_72h_rep2_3);

%plotQCsub(Sample_Mock_1_1,Sample_Mock_1_2,mean_Mock_1h_rep1_2,diff_Mock_1h_rep1_2);
%plotQCsub(Sample_Mock_1_1,Sample_Mock_1_3,mean_Mock_1h_rep1_3,diff_Mock_1h_rep1_3);
%plotQCsub(Sample_Mock_1_2,Sample_Mock_1_3,mean_Mock_1h_rep2_3,diff_Mock_1h_rep2_3);
%plotQCsub(Sample_Mock_3_1,Sample_Mock_3_2,mean_Mock_3h_rep1_2,diff_Mock_3h_rep1_2);
%plotQCsub(Sample_Mock_3_1,Sample_Mock_3_3,mean_Mock_3h_rep1_3,diff_Mock_3h_rep1_3);
%plotQCsub(Sample_Mock_3_2,Sample_Mock_3_3,mean_Mock_3h_rep2_3,diff_Mock_3h_rep2_3);
%plotQCsub(Sample_Mock_24_1,Sample_Mock_24_2,mean_Mock_24h_rep1_2,diff_Mock_24h_rep1_2);
%plotQCsub(Sample_Mock_24_1,Sample_Mock_24_3,mean_Mock_24h_rep1_3,diff_Mock_24h_rep1_3);
%plotQCsub(Sample_Mock_24_2,Sample_Mock_24_3,mean_Mock_24h_rep2_3,diff_Mock_24h_rep2_3);
%plotQCsub(Sample_Mock_72_1,Sample_Mock_72_2,mean_Mock_72h_rep1_2,diff_Mock_72h_rep1_2);
%plotQCsub(Sample_Mock_72_1,Sample_Mock_72_3,mean_Mock_72h_rep1_3,diff_Mock_72h_rep1_3);
%plotQCsub(Sample_Mock_72_2,Sample_Mock_72_3,mean_Mock_72h_rep2_3,diff_Mock_72h_rep2_3);


ods pdf close ;
ods graphics off;


proc export data=apn_sbys_p05 outfile="!PATCON/arabidopsis/analysis_output/arabidopsis_jmp_any_p05_genes_data.csv"
    dbms=tab replace;
run;

