ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Upper-quartile normalizing the fusion counts. Then, I'm also going to line center the data (like with
   CEGS1), using UQ/mean/median of the log normalized APN per treatment*time */

data coverage_data;
   set rs.counts_by_isoform;
run;

data coverage_data_tpm0;
   set coverage_data;
   if tpm gt 0;
run;

/* quartiles */
proc sort data=coverage_data_tpm0;
   by sample_id;
run;

proc univariate data=coverage_data_tpm0 noprint;
  by sample_id;
  var tpm;
  output out=quartiles_tpm0 Q3=Q3;
  run;

/* Calculate the Q3 of Q3 for normalizations */
proc means data=quartiles_tpm0 noprint;
            var q3;
            output out=q3_q3 q3=q3;
            run;

data _null_;
   set q3_q3;
   call symput('q3q3',q3);
run;

%put &q3q3.;

/* Q3 of Q3 = 13.915*/

/* Now we need to merge the stats back into the original dataset */

    proc sort data=coverage_data;
        by sample_id;
        run;

    proc sort data=quartiles_tpm0;
        by sample_id;
        run;

/* merge and oops test */

    data coverage_data_w_q3 oops_tpm0;
        merge coverage_data (in=in1) quartiles_tpm0 (in=in2);
        by sample_id;
        if in1 and in2 then output coverage_data_w_q3;
        else output oops_tpm0; 
        run;

/* Calculate the adjustments per sample */
/* Divide the across-sample UQ by the sample UQ and then output */
/* Then we want to export this as plot the distributions */

data iso_counts_q3_norm;
        set coverage_data_w_q3;
        log_tpm=log(tpm + 1);

        * Q3 normalization ;
        q3_q3_tpm=(tpm/q3) * &q3q3.;
        q3_q3_ff=&q3q3./q3;
        log_q3_q3_tpm=log(q3_q3_tpm + 1);
run;

/* Line centering */

* centering by genotype, treatment and time;
proc sort data=iso_counts_q3_norm;
   by treatment time;
proc means data=iso_counts_q3_norm noprint;
   by treatment time;
   var log_q3_q3_tpm;
   output out=means
   mean(log_q3_q3_tpm)=mean_log_q3_tpm
   median(log_q3_q3_tpm)=median_log_q3_tpm
   q3(log_q3_q3_tpm)=q3_log_q3_tpm;
run;

proc sort data=means;
   by treatment time;
run;

data iso_counts_q3_norm2;
   merge iso_counts_q3_norm (in=in1) means (in=in2);
   by treatment time;
   if in1;
   drop _type_ _freq_;
run;

/* Calculate centered values */

data iso_counts_q3_norm_cent;
   set iso_counts_q3_norm2;
   mean_log_q3_center=log_q3_q3_tpm - mean_log_q3_tpm;
   median_log_q3_center=log_q3_q3_tpm - median_log_q3_tpm;
   q3_log_q3_center=log_q3_q3_tpm - q3_log_q3_tpm;
run;

data rs.arab_xscript_counts_q3_norm_cent;
   set iso_counts_q3_norm_cent;
run;

/* Export data so we can make plots */


    proc sort data=iso_counts_q3_norm_cent;
        by sample_id;
        run;

data counts_q3_norm_cent_gt0;
   length condition $20.;
   set iso_counts_q3_norm_cent;
   condition=catx("_",treatment,time);
   if tpm gt 0;
run;


/* Export CSV for plotting */

proc export data= counts_q3_norm_cent_gt0 outfile='!PATCON/arabidopsis/analysis_output/transcripts_q3_norm_tpm_centered_gt0.csv' dbms=csv replace;
    putnames=yes;
    run;

