/* Running model -- by fusion, time: time*treatment, plate as random */

ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";


/* Get list of "on" fusions */
* as we are now doing contrasts that involve mock-only or treatment-only comparisons, I am changing;
* this from "on in all treatments" to "on in any treatment";
* Can always drop fusions later, but before FDR;

data on_xs;
   set rs.arab_flag_transcript_on_gt0;
   if flag_xscript_on_mock_apn0=1 or flag_xscript_on_01gy_apn0=1 or flag_xscript_on_1gy_apn0=1 ;
   keep transcript_id;
run;


/* Merge in with count data */


proc sort data=on_xs;
   by transcript_id;
proc sort data=rs.arab_xscript_counts_q3_norm_cent;
   by transcript_id;
run;

data xs_counts_on;
   merge rs.arab_xscript_counts_q3_norm_cent (in=in1) on_xs (in=in2);
   by transcript_id;
   if in1 and in2;
run;

/* Merge in design file
   Want variables: treatment, time */

/* Model for contrasts */

*Cat together genotype, time and treatment, then find order for contrasts;

data xs_counts_on2;
  set xs_counts_on;
  length time_trt $20.;
  time_trt = catx('_', time, treatment);
run;

proc freq data=xs_counts_on2 noprint;
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

proc sort data=xs_counts_on2;
   by transcript_id time;
   run;

ods listing close;


%macro runModels(measure,outname);

proc mixed data=xs_counts_on2;
  by transcript_id;
  class time treatment time_trt;
  model &measure. = time_trt / htype=1;

  lsmeans time_trt;
                  /* Order of contrasts:     1hr              24hr             3h               72h  
                                             0.1gy 1gy Mock   0.1gy 1gy Mock   0.1gy 1gy Mock   0.1gy 1gy Mock */
  contrast '0.1gy-Mock: 1h'  time_trt        1     0    -1    0     0   0      0     0   0      0     0   0  ;
  contrast '1gy-Mock: 1h'    time_trt        0     1    -1    0     0   0      0     0   0      0     0   0  ;
  contrast '1gy-0.1gy: 1h'   time_trt        -1    1     0    0     0   0      0     0   0      0     0   0  ;

  contrast '0.1gy-Mock: 24h' time_trt        0     0     0    1     0   -1     0     0   0      0     0   0  ;
  contrast '1gy-Mock: 24h'   time_trt        0     0     0    0     1   -1     0     0   0      0     0   0  ;
  contrast '1gy-0.1gy: 24h'  time_trt        0     0     0    -1    1   0      0     0   0      0     0   0  ;
 
  contrast '0.1gy-Mock: 3h'  time_trt        0     0     0    0     0   0      1     0  -1      0     0   0     ;   
  contrast '1gy-Mock: 3h'    time_trt        0     0     0    0     0   0      0     1  -1      0     0   0     ;   
  contrast '1gy-0.1gy: 3h'   time_trt        0     0     0    0     0   0     -1     1   0      0     0   0     ;   

  contrast '0.1gy-Mock: 72h' time_trt        0     0     0    0     0   0      0     0   0      1     0  -1    ;
  contrast '1gy-Mock: 72h'   time_trt        0     0     0    0     0   0      0     0   0      0     1  -1  ;
  contrast '1gy-0.1gy: 72h'  time_trt        0     0     0    0     0   0      0     0   0     -1     1   0  ;

  contrast 'Mock: 3h-1h'     time_trt        0     0    -1    0     0   0      0     0   1      0     0   0  ;
  contrast 'Mock: 24h-1h'    time_trt        0     0    -1    0     0   1      0     0   0      0     0   0  ;
  contrast 'Mock: 72h-1h'    time_trt        0     0    -1    0     0   0      0     0   0      0     0   1  ;
  contrast 'Mock: 24h-3h'    time_trt        0     0     0    0     0   1      0     0  -1      0     0   0  ;
  contrast 'Mock: 72h-3h'    time_trt        0     0     0    0     0   0      0     0  -1      0     0   1  ;
  contrast 'Mock: 72h-24h'   time_trt        0     0     0    0     0  -1      0     0   0      0     0   1  ;

  contrast '0.1gy: 3h-1h'    time_trt       -1     0     0    0     0   0      1     0   0      0     0   0   ;
  contrast '0.1gy: 24h-1h'   time_trt       -1     0     0    1     0   0      0     0   0      0     0   0  ;
  contrast '0.1gy: 72h-1h'   time_trt       -1     0     0    0     0   0      0     0   0      1     0   0  ;
  contrast '0.1gy: 24h-3h'   time_trt        0     0     0    1     0   0     -1     0   0      0     0   0  ;
  contrast '0.1gy: 72h-3h'   time_trt        0     0     0    0     0   0     -1     0   0      1     0   0  ;
  contrast '0.1gy: 72h-24h'  time_trt        0     0     0   -1     0   0      0     0   0      1     0   0  ;

  contrast '1gy: 3h-1h'      time_trt        0    -1     0    0     0   0      0     1   0      0     0   0  ;
  contrast '1gy: 24h-1h'     time_trt        0    -1     0    0     1   0      0     0   0      0     0   0  ;
  contrast '1gy: 72h-1h'     time_trt        0    -1     0    0     0   0      0     0   0      0     1   0  ;
  contrast '1gy: 24h-3h'     time_trt        0     0     0    0     1   0      0    -1   0      0     0   0  ;
  contrast '1gy: 72h-3h'     time_trt        0     0     0    0     0   0      0    -1   0      0     1   0  ;
  contrast '1gy: 72h-24h'    time_trt        0     0     0    0    -1   0      0     0   0      0     1   0  ;

contrast '0.1gy 3h-1h = Mock 3h-1h' time_trt    -1  0   1     0      0  0      1     0   -1     0     0    0 ; 
contrast '0.1gy 24h-1h = Mock 24h-1h' time_trt  -1  0   1     1      0  -1     0     0    0     0     0    0 ;
contrast '0.1gy 72h-1h = Mock 72h-1h' time_trt  -1  0   1     0      0  0      0      0  0      1     0   -1 ;

contrast '1gy 3h-1h = Mock 3h-1h' time_trt       0  -1  1     0      0  0      1     0   -1     0     0    0 ; 
contrast '1gy 24h-1h = Mock 24h-1h' time_trt     0  -1  1     1      0  -1     0     0    0     0     0    0 ;
contrast '1gy 72h-1h = Mock 72h-1h' time_trt     0  -1  1     0      0  0      0      0  0      1     0   -1 ;

  estimate '0.1gy-Mock: 1h'  time_trt        1     0    -1    0     0   0      0     0   0      0     0   0  ;
  estimate '1gy-Mock: 1h'    time_trt        0     1    -1    0     0   0      0     0   0      0     0   0  ;
  estimate '1gy-0.1gy: 1h'   time_trt        -1    1     0    0     0   0      0     0   0      0     0   0  ;

  estimate '0.1gy-Mock: 24h' time_trt        0     0     0    1     0   -1     0     0   0      0     0   0  ;
  estimate '1gy-Mock: 24h'   time_trt        0     0     0    0     1   -1     0     0   0      0     0   0  ;
  estimate '1gy-0.1gy: 24h'  time_trt        0     0     0    -1    1   0      0     0   0      0     0   0  ;
 
  estimate '0.1gy-Mock: 3h'  time_trt        0     0     0    0     0   0      1     0  -1      0     0   0     ;   
  estimate '1gy-Mock: 3h'    time_trt        0     0     0    0     0   0      0     1  -1      0     0   0     ;   
  estimate '1gy-0.1gy: 3h'   time_trt        0     0     0    0     0   0     -1     1   0      0     0   0     ;   

  estimate '0.1gy-Mock: 72h' time_trt        0     0     0    0     0   0      0     0   0      1     0  -1    ;
  estimate '1gy-Mock: 72h'   time_trt        0     0     0    0     0   0      0     0   0      0     1  -1  ;
  estimate '1gy-0.1gy: 72h'  time_trt        0     0     0    0     0   0      0     0   0     -1     1   0  ;

  estimate 'Mock: 3h-1h'     time_trt        0     0    -1    0     0   0      0     0   1      0     0   0  ;
  estimate 'Mock: 24h-1h'    time_trt        0     0    -1    0     0   1      0     0   0      0     0   0  ;
  estimate 'Mock: 72h-1h'    time_trt        0     0    -1    0     0   0      0     0   0      0     0   1  ;
  estimate 'Mock: 24h-3h'    time_trt        0     0     0    0     0   1      0     0  -1      0     0   0  ;
  estimate 'Mock: 72h-3h'    time_trt        0     0     0    0     0   0      0     0  -1      0     0   1  ;
  estimate 'Mock: 72h-24h'   time_trt        0     0     0    0     0  -1      0     0   0      0     0   1  ;

  estimate '0.1gy: 3h-1h'    time_trt       -1     0     0    0     0   0      1     0   0      0     0   0   ;
  estimate '0.1gy: 24h-1h'   time_trt       -1     0     0    1     0   0      0     0   0      0     0   0  ;
  estimate '0.1gy: 72h-1h'   time_trt       -1     0     0    0     0   0      0     0   0      1     0   0  ;
  estimate '0.1gy: 24h-3h'   time_trt        0     0     0    1     0   0     -1     0   0      0     0   0  ;
  estimate '0.1gy: 72h-3h'   time_trt        0     0     0    0     0   0     -1     0   0      1     0   0  ;
  estimate '0.1gy: 72h-24h'  time_trt        0     0     0   -1     0   0      0     0   0      1     0   0  ;

  estimate '1gy: 3h-1h'      time_trt        0    -1     0    0     0   0      0     1   0      0     0   0  ;
  estimate '1gy: 24h-1h'     time_trt        0    -1     0    0     1   0      0     0   0      0     0   0  ;
  estimate '1gy: 72h-1h'     time_trt        0    -1     0    0     0   0      0     0   0      0     1   0  ;
  estimate '1gy: 24h-3h'     time_trt        0     0     0    0     1   0      0    -1   0      0     0   0  ;
  estimate '1gy: 72h-3h'     time_trt        0     0     0    0     0   0      0    -1   0      0     1   0  ;
  estimate '1gy: 72h-24h'    time_trt        0     0     0    0    -1   0      0     0   0      0     1   0  ;

estimate '0.1gy 3h-1h = Mock 3h-1h' time_trt    -1  0   1     0      0  0      1     0   -1     0     0    0 ; 
estimate '0.1gy 24h-1h = Mock 24h-1h' time_trt  -1  0   1     1      0  -1     0     0    0     0     0    0 ;
estimate '0.1gy 72h-1h = Mock 72h-1h' time_trt  -1  0   1     0      0  0      0      0  0      1     0   -1 ;

estimate '1gy 3h-1h = Mock 3h-1h' time_trt       0  -1  1     0      0  0      1     0   -1     0     0    0 ; 
estimate '1gy 24h-1h = Mock 24h-1h' time_trt     0  -1  1     1      0  -1     0     0    0     0     0    0 ;
estimate '1gy 72h-1h = Mock 72h-1h' time_trt     0  -1  1     0      0  0      0      0  0      1     0   -1 ;


ods output tests1=tests1 lsmeans=lsmeans contrasts=contrasts estimates=estimates ;
run;
quit;


/* Make permenant */


data rs.arab_xs_cntrs_tests1_&outname. ;
  set tests1 ;
  run ;

data rs.arab_xs_cntrs_lsmeans_&outname. ;
  set lsmeans ;
  run ;

data rs.arab_xs_cntrs_constr_&outname. ;
  set contrasts ;
  run ;

data rs.arab_xs_cntrs_estim_&outname. ;
  set estimates ;
  run ;


%mend;

%runModels(log_q3_q3_tpm,non);
*%runModels(mean_log_q3_center,mean);
*%runModels(median_log_q3_center,med);
*%runModels(q3_log_q3_center,q3);
   


