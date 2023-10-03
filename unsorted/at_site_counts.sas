data meth_data;
  set wgbsA.methylation_data_cg_chg_chh;
  keep chr stop_pos total_C treatment units rep site_type;
run;

proc sort data=meth_data;
   by site_type chr stop_pos treatment units;
proc means data=meth_data noprint;
   by site_type chr stop_pos treatment units;
   var total_C;
   output out=meth_data_sum sum=;
run;

data gc_data;
  set wgbsA.methylation_data_gc;
  keep chr stop_pos total_C treatment units rep ;
  rename total_C = total_C_GC;
run;

proc sort data=gc_data;
   by  chr stop_pos treatment units;
proc means data=gc_data noprint;
   by  chr stop_pos treatment units;
   var total_C_GC;
   output out=gc_data_sum sum=;
run;

proc sort data=meth_data;
  by chr stop_pos treatment units rep;
proc sort data=gc_data;
  by chr stop_pos treatment units rep;
proc sort data=meth_data_sum;
  by chr stop_pos treatment units ;
proc sort data=gc_data_sum;
  by chr stop_pos treatment units ;
run;

data meth2gc_byrep;
  merge meth_data (in=in1) gc_data (in=in2);
  by chr stop_pos treatment units rep;
  if in1 then flag_endo=1; else flag_endo=0;
  if in2 then flag_gc=1; else flag_gc=0;
run;

data meth2gc_sum;
  merge meth_data_sum (in=in1) gc_data_sum (in=in2);
  by chr stop_pos treatment units;
  if in1 then flag_endo=1; else flag_endo=0;
  if in2 then flag_gc=1; else flag_gc=0;
run;


data meth2gc_byrep2;
  set meth2gc_byrep;
  if total_C_GC = . then total_C_highest=total_C;
  else if total_C = . then total_C_highest=total_C_GC;
  else if total_C_GC >= total_C then total_C_highest=total_C_GC;
  else if total_C_GC <= total_C then total_C_highest=total_C;
run;

data meth2gc_sum2;
  set meth2gc_sum;
  if total_C_GC = . then total_C_highest=total_C;
  else if total_C = . then total_C_highest=total_C_GC;
  else if total_C_GC >= total_C then total_C_highest=total_C_GC;
  else if total_C_GC <= total_C then total_C_highest=total_C;
run;

proc sort data=meth2gc_byrep2;
  by treatment  units rep;
proc means data=meth2gc_byrep2 noprint;
  where chr in ("1","2","3","4","5") and total_C_highest >= 10;
  by treatment units  rep;
  var total_C_highest;
  output out=mean_depth_rep mean=;
run;

proc sort data=meth2gc_sum2;
  by treatment  units ;
proc means data=meth2gc_sum2 noprint;
  where chr in ("1","2","3","4","5") and total_C_highest >= 10;
  by treatment units  ;
  var total_C_highest;
  output out=mean_depth_sum mean=;
run;


proc freq data=meth2gc_byrep2 noprint;
  where total_C >= 10 ;
  tables treatment*units*rep*flag_gc*flag_endo*site_type / out=ctabs_count_10x_rep ;
run;

proc freq data=meth2gc_sum2 noprint;
  where total_C >= 10;
  tables treatment*units*flag_gc*flag_endo*site_type / out=ctabs_count_10x_sum;
run;


proc freq data=meth2gc_byrep2 noprint;
  where total_C_GC >= 10 ;
  tables treatment*units*rep*flag_gc*flag_endo*site_type / out=ctabs_count_10x_rep_GC;
run;

proc freq data=meth2gc_sum2 noprint;
  where total_C_GC >= 10;
  tables treatment*units*flag_gc*flag_endo*site_type / out=ctabs_count_10x_sum_GC;
run;

proc print data=ctabs_count_10x_rep;
proc print data=ctabs_count_10x_sum;
proc print data=ctabs_count_10x_rep_gc;
proc print data=ctabs_count_10x_sum_gc;
run;


/*
ctabs_count_10x_rep

                                                 flag_    site_
treatment    units             rep    flag_gc     endo    type       COUNT

  01Gy        0U                 1       0         1       CG       2678814
  01Gy        0U                 1       0         1       CHG      2755727
  01Gy        0U                 1       0         1       CHH     12612514
  01Gy        0U                 1       1         1       CG        534053
  01Gy        0U                 1       1         1       CHG       754704
  01Gy        0U                 1       1         1       CHH      2835525
  01Gy        0U                 2       0         1       CG       2528511
  01Gy        0U                 2       0         1       CHG      2596568
  01Gy        0U                 2       0         1       CHH     12613990
  01Gy        0U                 2       1         1       CG        492017
  01Gy        0U                 2       1         1       CHG       692190
  01Gy        0U                 2       1         1       CHH      2697408
  0Gy         0U                 1       0         1       CG       2240455
  0Gy         0U                 1       0         1       CHG      2221660
  0Gy         0U                 1       0         1       CHH      8647684
  0Gy         0U                 1       1         1       CG        491253
  0Gy         0U                 1       1         1       CHG       675074
  0Gy         0U                 1       1         1       CHH      2228011
  0Gy         0U                 2       0         1       CG       1894404
  0Gy         0U                 2       0         1       CHG      1909141
  0Gy         0U                 2       0         1       CHH      8405465
  0Gy         0U                 2       1         1       CG        390827
  0Gy         0U                 2       1         1       CHG       535246
  0Gy         0U                 2       1         1       CHH      1944737
  1Gy         0U                 1       0         1       CG       1382288
  1Gy         0U                 1       0         1       CHG      1368344
  1Gy         0U                 1       0         1       CHH      6181487
  1Gy         0U                 1       1         1       CG        286239
  1Gy         0U                 1       1         1       CHG       382492
  1Gy         0U                 1       1         1       CHH      1402140
  1Gy         0U                 2       0         1       CG       2079264
  1Gy         0U                 2       0         1       CHG      2090119
  1Gy         0U                 2       0         1       CHH      9110857
  1Gy         0U                 2       1         1       CG        858835
  1Gy         0U                 2       1         1       CHG      1181788
  1Gy         0U                 2       1         1       CHH      4216505




ctabs_count_10x_sum


                                  flag_    site_
 treatment    units    flag_gc     endo    type       COUNT

   01Gy        0U         0         1       CG       3774394
   01Gy        0U         0         1       CHG      3918647
   01Gy        0U         0         1       CHH     19644385
   01Gy        0U         1         1       CG        716419
   01Gy        0U         1         1       CHG      1028863
   01Gy        0U         1         1       CHH      4058259
   0Gy         0U         0         1       CG       3270988
   0Gy         0U         0         1       CHG      3345322
   0Gy         0U         0         1       CHH     14907165
   0Gy         0U         1         1       CG        657251
   0Gy         0U         1         1       CHG       932258
   0Gy         0U         1         1       CHH      3407452
   1Gy         0U         0         1       CG       2978421
   1Gy         0U         0         1       CHG      3050624
   1Gy         0U         0         1       CHH     14251584
   1Gy         0U         1         1       CG        590021
   1Gy         0U         1         1       CHG       832270
   1Gy         0U         1         1       CHH      3125834



ctabs_count_10x_rep_gc


                                                  flag_    site_
 treatment    units             rep    flag_gc     endo    type       COUNT

   01Gy       0U                  1       1         0                     2
   01Gy       0U                  1       1         1       CG       426114
   01Gy       0U                  1       1         1       CHG      611764
   01Gy       0U                  1       1         1       CHH     2411613
   01Gy       0U                  2       1         0                     1
   01Gy       0U                  2       1         1       CG       410023
   01Gy       0U                  2       1         1       CHG      589546
   01Gy       0U                  2       1         1       CHH     2324583
   01Gy       100U                1       1         0               3483007
   01Gy       100U                2       1         0               3776572
   0Gy        0U                  1       1         0                     1
   0Gy        0U                  1       1         1       CG       419367
   0Gy        0U                  1       1         1       CHG      595860
   0Gy        0U                  1       1         1       CHH     2173724
   0Gy        0U                  2       1         0                     2
   0Gy        0U                  2       1         1       CG       329875
   0Gy        0U                  2       1         1       CHG      466242
   0Gy        0U                  2       1         1       CHH     1708711
   0Gy        100U                1       1         0               2080359
   0Gy        100U                2       1         0               2045400
   1Gy        0U                  1       1         0                     1
   1Gy        0U                  1       1         1       CG       285345
   1Gy        0U                  1       1         1       CHG      405017
   1Gy        0U                  1       1         1       CHH     1515218
   1Gy        0U                  2       1         0                     1
   1Gy        0U                  2       1         1       CG       642438
   1Gy        0U                  2       1         1       CHG      890687
   1Gy        0U                  2       1         1       CHH     3240329
   1Gy        100U                1       1         0               4402496
   1Gy        100U                2       1         0               4289662

                                 The SAS System



ctabs_count_10x_sum_gc




*/





