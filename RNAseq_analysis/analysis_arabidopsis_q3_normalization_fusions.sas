ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Upper-quartile normalizing the fusion counts. Then, I'm also going to line center the data (like with
   CEGS1), using UQ/mean/median of the log normalized APN per treatment*time */

data coverage_data;
   set rs.arab_fusion_counts_w_key;
run;

data coverage_data_apn0;
   set coverage_data;
   if apn gt 0;
run;

/* quartiles */
proc sort data=coverage_data_apn0;
   by sample_number;
run;

proc univariate data=coverage_data_apn0 noprint;
  by sample_number;
  var apn;
  output out=quartiles_apn0 Q3=Q3;
  run;

/* Calculate the Q3 of Q3 for normalizations */
proc means data=quartiles_apn0 noprint;
            var q3;
            output out=q3_q3 q3=q3;
            run;

data _null_;
   set q3_q3;
   call symput('q3q3',q3);
run;

%put &q3q3.;

/* Q3 of Q3 = 120.83228417 */

/* Now we need to merge the stats back into the original dataset */

    proc sort data=coverage_data;
        by sample_number;
        run;

    proc sort data=quartiles_apn0;
        by sample_number;
        run;

/* merge and oops test */

    data coverage_data_w_q3 oops_apn0;
        merge coverage_data (in=in1) quartiles_apn0 (in=in2);
        by sample_number;
        if in1 and in2 then output coverage_data_w_q3;
        else output oops_apn0; 
        run;

/* Calculate the adjustments per sample */
/* Divide the across-sample UQ by the sample UQ and then output */
/* Then we want to export this as plot the distributions */

data fusion_counts_q3_norm;
        set coverage_data_w_q3;
        log_apn=log(apn + 1);

        * Q3 normalization ;
        q3_q3_apn=(apn/q3) * &q3q3.;
        q3_q3_ff=&q3q3./q3;
        log_q3_q3_apn=log(q3_q3_apn + 1);
run;

/* Line centering */

*Get design file, as I need treatment, genotype and timepoint;

data sample_key;
  set rs.arab_design_file;
  length condition $20.;
  condition=catx("_",treatment,time);
  keep sample_number condition treatment time;
run;

proc sort data=sample_key nodup;
   by sample_number;
proc sort data=fusion_counts_q3_norm;
   by sample_number;
run;

data norm_q3_w_key;
   merge sample_key (in=in1) fusion_counts_q3_norm (in=in2);
   by sample_number;
   if in1 and in2;
run;

* centering by genotype, treatment and time;
proc sort data=norm_q3_w_key;
   by treatment time;
proc means data=norm_q3_w_key noprint;
   by treatment time;
   var log_q3_q3_apn;
   output out=means
   mean(log_q3_q3_apn)=mean_log_q3_apn
   median(log_q3_q3_apn)=median_log_q3_apn
   q3(log_q3_q3_apn)=q3_log_q3_apn;
run;

proc sort data=means;
   by treatment time;
run;

data norm_q3_w_key2;
   merge norm_q3_w_key (in=in1) means (in=in2);
   by treatment time;
   if in1;
   drop _type_ _freq_;
run;

/* Calculate centered values */

data fusion_counts_q3_norm_cent;
   set norm_q3_w_key2;
   mean_log_q3_center=log_q3_q3_apn - mean_log_q3_apn;
   median_log_q3_center=log_q3_q3_apn - median_log_q3_apn;
   q3_log_q3_center=log_q3_q3_apn - q3_log_q3_apn;
run;

data rs.arab_fusion_counts_q3_norm_cent;
   set fusion_counts_q3_norm_cent;
run;

/* Export data so we can make plots */


    proc sort data=fusion_counts_q3_norm_cent;
        by sample_number;
        run;

data counts_q3_norm_cent_gt0;
   set fusion_counts_q3_norm_cent;
   if apn gt 0;
run;


/* Export CSV for plotting */

proc export data= counts_q3_norm_cent_gt0 outfile='!PATCON/arabidopsis/analysis_output/fusions_q3_norm_apn_centered_gt0.csv' dbms=csv replace;
    putnames=yes;
    run;

