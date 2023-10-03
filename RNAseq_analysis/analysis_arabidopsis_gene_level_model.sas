/* Calculating gene-level coverage, flagging on/off at APN>0 and re-running models */

ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";


/* Going to calculate the mean coverage across all fusions for a gene, skipping   
   multigene fusions. So we won't be able to get some genes, but this is okay for now */

data fus2gene;
  set tair.tair20_fusion_si_info_unique;
  where flag_multigene=0;
  keep fusion_id gene_id;
run;

/* Merge with normalized data -- any gene with a normalized APN of 0 are removed */

data q3_counts;
  set rs.arab_fusion_counts_q3_norm_cent;
  keep sample_number fusion_id q3_q3_apn;
run;

proc sort data=fus2gene nodup;
  by fusion_id gene_id;
proc sort data=q3_counts;
  by fusion_id;
run;

data q3_counts2keep;
  merge fus2gene (in=in1) q3_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=q3_counts2keep;
  by sample_number gene_id;
proc means data=q3_counts2keep noprint;
  by sample_number gene_id;
  var q3_q3_apn;
  output out=mean_q3_apn_by_gene_sample mean=;
run;


data rs.arab_mean_q3_apn_by_gene;
  set mean_q3_apn_by_gene_sample;
run;



/* flagging on/off at APN>0 */


data design;
  set rs.arab_design_file;
  keep sample_number treatment time replicate;
run;

proc sort data=design nodups;
  by sample_number;
proc sort data=mean_q3_apn_by_gene_sample;
  by sample_number;
run;

data gene_counts_w_key;
  merge design (in=in1) mean_q3_apn_by_gene_sample (in=in2);
  by sample_number;
  if in1 and in2;
  if q3_q3_apn > 0 then flag_apn_gt0=1;
  else flag_apn_gt0=0;
run;

/* Flag on per treatment */

proc sort data=gene_counts_w_key;
   by treatment gene_id;
proc means data=gene_counts_w_key noprint;
   by treatment gene_id;
   var flag_apn_gt0;
   output out=mean_on_by_trt mean=;
run;

data flag_gene_on_by_trt;
  set mean_on_by_trt;
  /* For APN > 0 */
  if flag_apn_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_apn_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;
  keep gene_id treatment flag_gene_on_gt0;
run;

proc sort data=flag_gene_on_by_trt;
   by gene_id treatment;
proc transpose data=flag_gene_on_by_trt out=flag_gene_trt_sbys_apn0;
   by gene_id;
   id treatment;
   var flag_gene_on_gt0;
run;

data treat_apn0;
   set flag_gene_trt_sbys_apn0;
   drop _NAME_;
   rename _0_1gy=flag_gene_on_01gy_apn0
          _1gy=flag_gene_on_1gy_apn0
          Mock=flag_gene_on_mock_apn0;
run;

/* Flag on per time */

proc sort data=gene_counts_w_key;
   by time gene_id;
proc means data=gene_counts_w_key noprint;
   by time gene_id;
   var flag_apn_gt0;
   output out=mean_on_by_time mean=;
run;

data flag_gene_on_by_time;
  set mean_on_by_time;

  /* For APN > 0 */
  if flag_apn_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_apn_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  keep gene_id time flag_gene_on_gt0 ;
run;

proc sort data=flag_gene_on_by_time;
   by gene_id time;
proc transpose data=flag_gene_on_by_time out=flag_gene_time_sbys_apn0;
   by gene_id;
   id time;
   var flag_gene_on_gt0;
run;


data time_apn0;
   set flag_gene_time_sbys_apn0;
   drop _NAME_;
   rename _1=flag_gene_on_1hr_apn0
          _3=flag_gene_on_3hr_apn0
          _24=flag_gene_on_34hr_apn0
          _72=flag_gene_on_72hr_apn0;
run;

/* Flag on per treatment-by-time */

proc sort data=gene_counts_w_key;
   by treatment time gene_id;
proc means data=gene_counts_w_key noprint;
   by treatment time gene_id;
   var flag_apn_gt0 ;
   output out=mean_on_by_trttime mean=;
run;

data flag_gene_on_by_trttime;
  length condition $10.;
  set mean_on_by_trttime;
  condition=catx("_",treatment,time);

  /* For APN > 0 */
  if flag_apn_gt0 ge 0.5 then flag_gene_on_gt0=1;
  else if flag_apn_gt0 = 0 then flag_gene_on_gt0=0;
  else flag_gene_on_gt0=.;

  keep gene_id condition flag_gene_on_gt0;
run;

proc sort data=flag_gene_on_by_trttime;
   by gene_id condition;
proc transpose data=flag_gene_on_by_trttime out=flag_gene_trttime_sbys_apn0;
   by gene_id;
   id condition;
   var flag_gene_on_gt0;
run;


data trttime_apn0;
   set flag_gene_trttime_sbys_apn0;
   drop _NAME_;
   rename _0_1gy_1=flag_gene_on_01gy_1hr_apn0
          _0_1gy_3=flag_gene_on_01gy_3hr_apn0
          _0_1gy_24=flag_gene_on_01gy_24hr_apn0
          _0_1gy_72=flag_gene_on_01gy_72hr_apn0
          _1gy_1=flag_gene_on_1gy_1hr_apn0
          _1gy_3=flag_gene_on_1gy_3hr_apn0
          _1gy_24=flag_gene_on_1gy_24hr_apn0
          _1gy_72=flag_gene_on_1gy_72hr_apn0
          Mock_1=flag_gene_on_mock_1hr_apn0
          Mock_3=flag_gene_on_mock_3hr_apn0
          Mock_24=flag_gene_on_mock_24hr_apn0
          Mock_72=flag_gene_on_mock_72hr_apn0;
run;


/* Merge and make permenant */

proc sort data=treat_apn0;
   by gene_id;
proc sort data=time_apn0;
   by gene_id;
proc sort data=trttime_apn0;
   by gene_id;
run;

data rs.arab_flag_gene_on_gt0;
  merge treat_apn0 (in=in1) time_apn0 (in=in2) trttime_apn0 (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;

/*  Re-runimodels on gene */

data on_genes;
   set rs.arab_flag_gene_on_gt0;
   if flag_gene_on_mock_apn0=1 or flag_gene_on_01gy_apn0=1 or flag_gene_on_1gy_apn0=1 ;
   keep gene_id;
run;

/* Merge in with count data */

proc sort data=on_genes;
   by gene_id;
proc sort data=mean_q3_apn_by_gene_sample;
   by gene_id;
run;

data gene_counts_on;
   merge mean_q3_apn_by_gene_sample (in=in1) on_genes (in=in2);
   by gene_id;
   if in1 and in2;
run;

/* Merge in design file
   Want variables: treatment, time */
   
data sample_key;
  set rs.arab_design_file;
  keep sample_number treatment time;
run;

proc sort data=sample_key nodup;
   by sample_number;
proc sort data=gene_counts_on;
   by sample_number;
run;

data genes_on_w_key;
   merge sample_key (in=in1) gene_counts_on (in=in2);
   by sample_number;
   if in1 and in2;
   log_q3_q3_apn=log(q3_q3_apn+1);
   run;

/* Run model */

* Need to sort by fusion and time!;


data genes_on_w_key2;
  set genes_on_w_key;
  length time_trt $20.;
  time_trt = catx('_', time, treatment);
run;

proc freq data=genes_on_w_key2 noprint;
  tables time_trt / out=order;
run;


proc sort data=genes_on_w_key2;
  by gene_id time_trt;
run;

ods listing close;
proc glimmix data=genes_on_w_key2;
  by gene_id;
  class time treatment ;
  model log_q3_q3_apn = time|treatment / htype=1;
  output out=resid resid=resid pred=pred student=stu;
ods output tests1=anova;
run;
quit;

/* Flag residuals */

proc univariate data = resid normal noprint;
  by gene_id;
  var Resid;
  output out = normtest probn=pnorm;
  run;

data flag_resids;
  set normtest;
  if pnorm = . then flag_fail_norm = .;
        else if pnorm le 0.05 then flag_fail_norm = 1;
        else flag_fail_norm = 0;
  run;

proc freq data = flag_resids noprint;
  tables flag_fail_norm / out=genes_flag_fail_norm;
  run;

/* Make permenant */

data rs.arab_resid_main_trt_tm_genes;
  set flag_resids;
  run;
data rs.arab_anova_main_trt_tm_genes ;
  set anova ;
  run ;    
data rs.arab_norm_main_trt_tm_genes ;
  set genes_flag_fail_norm ;
  run ;


/* Contrasts */



proc mixed data=genes_on_w_key2;
  by gene_id;
  class time treatment time_trt;
  model log_q3_q3_apn = time_trt / htype=1;

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


data rs.arab_gene_cntrs_tests1 ;
  set tests1 ;
  run ;

data rs.arab_gene_cntrs_lsmeans ;
  set lsmeans ;
  run ;

data rs.arab_gene_cntrs_constr ;
  set contrasts ;
  run ;

data rs.arab_gene_cntrs_estim ;
  set estimates ;
  run ;



   


