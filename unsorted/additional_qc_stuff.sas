/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

/* Additional data QC :

(1) Rep-to-rep concordance for each condition:
    RNA: { 0, 10, 100 } cGy * { 1, 3, 24, 72 } hr
    METH: { 0, 10, 100 } cGy * { CG 0U, HCG 0U, CHG 0U, HCHG 0U, CHH 0U, HCHH 0U, GC 0U, GCH 0U, GC 100U, GCH 100U } site*units
        unwindowed, 10bp window, 100bp window, 1000bp window

(2) Distribution of methylation:
    METH: { 0, 10, 100 } cGy * { CG 0U, HCG 0U, CHG 0U, HCHG 0U, CHH 0U, HCHH 0U, GC 0U, GCH 0U, GC 100U, GCH 100U } site*units
        unwindowed, 10bp window, 100bp window, 1000bp window

(3) Bisulfite conversion ratio: HCHH 100U, CHH 100U, HCHH 0U, CHH 0U on Pt sites only 

(4) Gene and TE methylation/accessibilty plots (+/- 2kb)

(5) DE vs non DE gene methylation/accessibility plots (+/2 2kb)
    -> do by time 
*/




/*

(1) Rep-to-rep concordance for each condition:
    RNA: { 0, 10, 100 } cGy * { 1, 3, 24, 72 } hr
    METH: { 0, 10, 100 } cGy * { CG 0U, HCG 0U, CHG 0U, HCHG 0U, CHH 0U, HCHH 0U, GC 0U, GCH 0U, GC 100U, GCH 100U } site*units
        unwindowed, 10bp window, 100bp window, 1000bp window
*/



/* Bin sites into 100bp windows */

data meth_data_cg_chg_chh;
  set arabMAP.methylation_data_cg_chg_chh;
  length sample_id $32.;
  sample_id=catx("_",treatment,units,rep);
  chr_bin=int(stop_pos/100) + 1;
  if total_C < 100 then flag_lt100_reads=1;
  else flag_lt100_reads=0;
  keep chr stop_pos perc_methyl sample_id site_type chr_bin flag_lt100_reads;
run;

data meth_data_gc;
  set arabMAP.methylation_data_gc;
  where  flag_normalized=1;
  length sample_id $32.;
  sample_id=catx("_",treatment,units,rep);
  chr_bin=int(stop_pos/100) + 1;
  if total_C < 100 then flag_lt100_reads=1;
  else flag_lt100_reads=0;
  keep chr stop_pos perc_methyl_norm sample_id site_type chr_bin flag_lt100_reads;
  rename perc_methyl_norm=perc_methyl;
run;

data gc_sites;
  set arabMAP.methylation_data_gc;
  keep chr stop_pos;
run;

data endo_sites;
  set arabMAP.methylation_data_cg_chg_chh;
  keep chr stop_pos;
run;

proc sort data=meth_data_cg_chg_chh;
  by chr stop_pos;
proc sort data=meth_data_gc;
  by chr stop_pos;
proc sort data=gc_sites nodup;
  by chr stop_pos;
proc sort data=endo_sites nodup;
  by chr stop_pos;
run;


data meth_data_cg_chg_chh2;
  merge meth_data_cg_chg_chh (in=in1) gc_sites (in=in2);
  by chr stop_pos;
  if in2 then delete;
run;

data meth_data_gc2;
  merge meth_data_gc (in=in1) endo_sites (in=in2);
  by chr stop_pos;
  if in2 then delete;
run;


proc sort data=meth_data_cg_chg_chh;
  by sample_id site_type chr chr_bin;
proc means data=meth_data_cg_chg_chh noprint;
  by sample_id site_type chr chr_bin;
  var perc_methyl;
  output out=mean_methyl_by_bin_endo mean=;
run;

proc sort data=meth_data_cg_chg_chh2;
  by sample_id site_type chr chr_bin;
proc means data=meth_data_cg_chg_chh2 noprint;
  by sample_id site_type chr chr_bin;
  var perc_methyl;
  output out=mean_methyl_by_bin_endo2 mean=;
run;

proc sort data=meth_data_gc;
  by sample_id site_type chr chr_bin;
proc means data=meth_data_gc noprint;
  by sample_id site_type chr chr_bin;
  var perc_methyl;
  output out=mean_methyl_by_bin_gc mean=;
run;

proc sort data=meth_data_gc2;
  by sample_id site_type chr chr_bin;
proc means data=meth_data_gc2 noprint;
  by sample_id site_type chr chr_bin;
  var perc_methyl;
  output out=mean_methyl_by_bin_gc2 mean=;
run;



proc means data=meth_data_cg_chg_chh noprint;
  by sample_id site_type chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl;
  output out=mean_methyl_by_bin_gt100_endo mean=;
run;


proc means data=meth_data_cg_chg_chh2 noprint;
  by sample_id site_type chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl;
  output out=mean_methyl_by_bin_gt100_endo2 mean=;
run;

proc means data=meth_data_gc noprint;
  by sample_id site_type chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl;
  output out=mean_methyl_by_bin_gt100_gc mean=;
run;


proc means data=meth_data_gc2 noprint;
  by sample_id site_type chr chr_bin;
  where flag_lt100_reads=0;
  var perc_methyl;
  output out=mean_methyl_by_bin_gt100_gc2 mean=;
run;




/* Calc correlations */
ods graphics / ANTIALIASMAX=50000000 MAXOBS=50000000 ;



/* Endogenous sites */ 
%macro corrPlot(sample1,sample2);

data sample1_all;
  set mean_methyl_by_bin_endo;
  where sample_id = "&sample1.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample1_10;
  set mean_methyl_by_bin_gt100_endo;
  where sample_id = "&sample1.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample2_all;
  set mean_methyl_by_bin_endo;
  where sample_id = "&sample2.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample2.;
run;

data sample2_10;
  set mean_methyl_by_bin_gt100_endo;
  where sample_id = "&sample2." ;
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample2.;
run;


data sample1_all_2;
  set mean_methyl_by_bin_endo2;
  where sample_id = "&sample1.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample1_10_2;
  set mean_methyl_by_bin_gt100_endo2;
  where sample_id = "&sample1.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample2_all_2;
  set mean_methyl_by_bin_endo2;
  where sample_id = "&sample2.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample2.;
run;

data sample2_10_2;
  set mean_methyl_by_bin_gt100_endo2;
  where sample_id = "&sample2." ;
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample2.;
run;

run;


proc sort data=sample1_all;
  by site_type chr chr_bin ;
  proc sort data=sample1_10;
  by site_type chr chr_bin ;
proc sort data=sample2_all;
  by site_type chr chr_bin ;
proc sort data=sample2_10;
  by site_type chr chr_bin ;
proc sort data=sample1_all_2;
  by site_type chr chr_bin ;
  proc sort data=sample1_10_2;
  by site_type chr chr_bin ;
proc sort data=sample2_all_2;
  by site_type chr chr_bin;
proc sort data=sample2_10_2;
  by site_type chr chr_bin ;
run;

data sample_1v2_all;
  merge sample1_all (in=in1) sample2_all (in=in2);
  by site_type chr chr_bin ;
run;

data sample_1v2_10;
  merge sample1_10 (in=in1) sample2_10 (in=in2);
  by site_type chr chr_bin ;
run;

data sample_1v2_all_2;
  merge sample1_all_2 (in=in1) sample2_all_2 (in=in2);
  by site_type chr chr_bin ;
run;

data sample_1v2_10_2;
  merge sample1_10_2 (in=in1) sample2_10_2 (in=in2);
  by site_type chr chr_bin ;
run;



title "Correlation endogenous &sample1. &sample2. (all sites) ";
ods text="Correlation endogenous &sample1. &sample2. (all sites) ";

proc corr data=sample_1v2_all pearson ;
   by site_type;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;
W
proc sgplot data=sample_1v2_all;
   by site_type;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;


title "Correlation endogenous &sample1. &sample2. (100X reads) ";
ods text="Correlation endogenous &sample1. &sample2. (100X reads) ";

proc corr data=sample_1v2_10 pearson ;
   by site_type;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_10;
   by site_type;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;









title "Correlation non-GC endogenous &sample1. &sample2. (all sites) ";
ods text="Correlation non-GC endogenous &sample1. &sample2. (all sites) ";

proc corr data=sample_1v2_all_2 pearson;
   by site_type;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_all_2;
   by site_type;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;


title "Correlation non-GC endogenous &sample1. &sample2. (100X reads) ";
ods text="Correlation non-GC endogenous  &sample1. &sample2. (100X reads) ";

proc corr data=sample_1v2_10_2 pearson ;
   by site_type;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_10_2;
   by site_type;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;

%mend;

%corrPlot(0Gy_0U_1,0Gy_0U_2);
%corrPlot(01Gy_0U_1,01Gy_0U_2);
%corrPlot(1Gy_0U_1,1Gy_0U_2);




/* GC sites */ 
%macro corrPlot(sample1,sample2);


data sample1_all_gc;
  set mean_methyl_by_bin_gc;
  where sample_id = "&sample1.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample1_10_gc;
  set mean_methyl_by_bin_gt100_gc;
  where sample_id = "&sample1.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample1.;
run;

data sample2_all_gc;
  set mean_methyl_by_bin_gc;
  where sample_id = "&sample2.";
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample2.;
run;

data sample2_10_gc;
  set mean_methyl_by_bin_gt100_gc;
  where sample_id = "&sample2." ;
  keep chr chr_bin perc_methyl site_type;
  rename perc_methyl=perc_meth_&sample2.;
run;


proc sort data=sample1_all_gc;
  by site_type chr chr_bin ;
  proc sort data=sample1_10_gc;
  by site_type chr chr_bin ;
proc sort data=sample2_all_gc;
  by site_type chr chr_bin ;
proc sort data=sample2_10_gc;
  by site_type chr chr_bin ;
run;


data sample_1v2_all_gc;
  merge sample1_all_gc (in=in1) sample2_all_gc (in=in2);
  by site_type chr chr_bin ;
run;

data sample_1v2_10_gc;
  merge sample1_10_gc (in=in1) sample2_10_gc (in=in2);
  by site_type chr chr_bin ;
run;


title "Correlation GC &sample1. &sample2. (all sites) ";
ods text="Correlation GC &sample1. &sample2. (all sites) ";

proc corr data=sample_1v2_all_gc pearson;
   by site_type;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_all_gc;
   by site_type;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;


title "Correlation GC &sample1. &sample2. (100X reads) ";
ods text="Correlation GC &sample1. &sample2. (100X reads) ";

proc corr data=sample_1v2_10_gc pearson;
   by site_type;
  var perc_meth_&sample1. perc_meth_&sample2.;
run;

proc sgplot data=sample_1v2_10_gc;
   by site_type;
   scatter x=perc_meth_&sample1. y=perc_meth_&sample2.;
run;

%mend;

%corrPlot(0Gy_0U_1,0Gy_0U_2);
%corrPlot(01Gy_0U_1,01Gy_0U_2);
%corrPlot(1Gy_0U_1,1Gy_0U_2);

%corrPlot(0Gy_100U_1,0Gy_100U_2);
%corrPlot(01Gy_100U_1,01Gy_100U_2);
%corrPlot(1Gy_100U_1,1Gy_100U_2);

%corrPlot(0Gy_0U_1,0Gy_100U_1);
%corrPlot(01Gy_0U_1,01Gy_100U_1);
%corrPlot(1Gy_0U_1,1Gy_100U_1);
%corrPlot(0Gy_0U_2,0Gy_100U_2);
%corrPlot(01Gy_0U_2,01Gy_100U_2);
%corrPlot(1Gy_0U_2,1Gy_100U_2);







%corrPlot(01Gy_0U_1,01Gy_100U_1);
%corrPlot(01Gy_0U_1,01Gy_100U_2);
%corrPlot(01Gy_0U_1,1Gy_100U_1);
%corrPlot(01Gy_0U_1,1Gy_100U_2);
%corrPlot(01Gy_0U_1,0Gy_100U_1);
%corrPlot(01Gy_0U_1,0Gy_100U_2);

%corrPlot(01Gy_0U_2,01Gy_100U_1);
%corrPlot(01Gy_0U_2,01Gy_100U_2);
%corrPlot(01Gy_0U_2,1Gy_100U_1);
%corrPlot(01Gy_0U_2,1Gy_100U_2);
%corrPlot(01Gy_0U_2,0Gy_100U_1);
%corrPlot(01Gy_0U_2,0Gy_100U_2);






/*
(2) Distribution of methylation:
    METH: { 0, 10, 100 } cGy * { CG 0U, HCG 0U, CHG 0U, HCHG 0U, CHH 0U, HCHH 0U, GC 0U, GCH 0U, GC 100U, GCH 100U } site*units
        
*/

%macro methDistrib(inData);

proc sort data=&inData.;
  by sample_id site_type;
run;


proc means data=&inData. noprint;
  by sample_id site_type;
  var perc_methyl;
  output out=distrib_meth mean=mean stddev=sd min=min 
         p5=p5 p10=p10 p25=p25 median=median p75=p75 p90=p90 p95=p95 max=max;
run;

proc print data=distrib_meth;
run;

%mend;

%methDistrib(meth_data_cg_chg_chh);
%methDistrib(meth_data_cg_chg_chh2);
%methDistrib(meth_data_gc);
%methDistrib(mean_methyl_by_bin_endo);
%methDistrib(mean_methyl_by_bin_endo2);
%methDistrib(mean_methyl_by_bin_gc);



/*
                       site_
   Obs    sample_id    type     _TYPE_     _FREQ_             mean              sd             min              p5             p10

     1    01Gy_0U_1     CG         0       5395457    0.2954217998    0.3982049017               0               0               0
     2    01Gy_0U_1     CHG        0       5918139    0.1225198894    0.2281658457               0               0               0
     3    01Gy_0U_1     CHH        0      29822894    0.0617413646     0.128418806               0               0               0
     4    01Gy_0U_2     CG         0       5386982    0.2927996144     0.397248944               0               0               0
     5    01Gy_0U_2     CHG        0       5911491     0.120110648    0.2279616387               0               0               0
     6    01Gy_0U_2     CHH        0      29895274    0.0589234309    0.1252768352               0               0               0
     7    0Gy_0U_1      CG         0       5290045    0.2971661009     0.406668082               0               0               0
     8    0Gy_0U_1      CHG        0       5793848    0.1234464357    0.2416618073               0               0               0
     9    0Gy_0U_1      CHH        0      28312150    0.0626678721    0.1506402556               0               0               0
    10    0Gy_0U_2      CG         0       5308130     0.296521366    0.4059851484               0               0               0
    11    0Gy_0U_2      CHG        0       5820939    0.1237781441    0.2428255171               0               0               0
    12    0Gy_0U_2      CHH        0      28857942    0.0630722975    0.1493140926               0               0               0
    13    1Gy_0U_1      CG         0       5202513    0.2934684765    0.4090489894               0               0               0
    14    1Gy_0U_1      CHG        0       5706510    0.1196544535    0.2466606699               0               0               0
    15    1Gy_0U_1      CHH        0      28317320    0.0569889272    0.1484919529               0               0               0
    16    1Gy_0U_2      CG         0       5325151    0.3038011795    0.4092040459               0               0               0


   Obs             p25          median             p75             p90             p95             max

     1               0            0.04       0.7777778      0.95238096               1               1
     2               0               0      0.11111111             0.5       0.7111111               1
     3               0               0      0.07692308      0.18181819      0.33333334               1
     4               0     0.033333335       0.7647059      0.95652175               1               1
     5               0               0      0.11111111             0.5       0.7058824               1
     6               0               0     0.071428575      0.16666667             0.3               1
     7               0               0             0.8      0.96153843               1               1
     8               0               0             0.1      0.53333336            0.75               1
     9               0               0     0.050847456             0.2      0.36666667               1
    10               0               0       0.7878788               1               1               1
    11               0               0      0.10526316        0.516129            0.75               1
    12               0               0     0.055555556             0.2       0.3529412               1
    13               0               0      0.77272725               1               1               1
    14               0               0      0.09090909             0.5            0.75               1
    15               0               0               0             0.2      0.33333334               1
    16               0               0             0.8               1               1               1

                                    Correlation non-GC endogenous 0Gy_0U_1 0Gy_0U_2 (100X reads)                     11:27 Thursday, March 24,

            sample_     site_
     Obs       id       type     _TYPE_     _FREQ_             mean              sd             min              p5             p10

      17    1Gy_0U_2     CHG        0       5837982     0.135427231    0.2575292686               0               0               0
      18    1Gy_0U_2     CHH        0      29209026    0.1054795774    0.2227479085               0               0               0


     Obs             p25          median             p75             p90             p95             max

      17               0               0      0.11111111      0.59090906       0.7878788               1
      18               0               0      0.09090909           0.375       0.6666667               1


                     site_
 Obs    sample_id    type     _TYPE_     _FREQ_             mean              sd             min              p5             p10

   1    01Gy_0U_1     CG         0       4558197    0.2906846907     0.395947358               0               0               0
   2    01Gy_0U_1     CHG        0       4706899    0.1261895448    0.2338202141               0               0               0
   3    01Gy_0U_1     CHH        0      24952437    0.0630058894    0.1313610864               0               0               0
   4    01Gy_0U_2     CG         0       4551744    0.2879187267     0.394932977               0               0               0
   5    01Gy_0U_2     CHG        0       4702813    0.1234822293    0.2332743379               0               0               0
   6    01Gy_0U_2     CHH        0      25027249    0.0599300431    0.1276986504               0               0               0
   7    0Gy_0U_1      CG         0       4461195    0.2927740227    0.4048802321               0               0               0
   8    0Gy_0U_1      CHG        0       4595833    0.1277331507    0.2483632566               0               0               0
   9    0Gy_0U_1      CHH        0      23568931    0.0645321627    0.1550067818               0               0               0
  10    0Gy_0U_2      CG         0       4480550    0.2918845072    0.4039293897               0               0               0
  11    0Gy_0U_2      CHG        0       4624488    0.1276189273    0.2485875491               0               0               0
  12    0Gy_0U_2      CHH        0      24079577    0.0647947828    0.1530384737               0               0               0
  13    1Gy_0U_1      CG         0       4389524    0.2887377504    0.4069537405               0               0               0
  14    1Gy_0U_1      CHG        0       4532647    0.1233536246    0.2523104704               0               0               0
  15    1Gy_0U_1      CHH        0      23629068    0.0583312925    0.1515740521               0               0               0
  16    1Gy_0U_2      CG         0       4496371    0.3002746769     0.407639908               0               0               0


 Obs             p25          median             p75             p90             p95             max

   1               0     0.037037037       0.7647059            0.95               1               1
   2               0               0      0.11111111             0.5      0.72727275               1
   3               0               0      0.07692308             0.2      0.33333334               1
   4               0     0.030303031            0.75      0.95238096               1               1
   5               0               0      0.11111111             0.5       0.7222222               1
   6               0               0     0.071428575      0.18181819          0.3125               1
   7               0               0             0.8      0.96153843               1               1
   8               0               0             0.1       0.5652174       0.7619048               1
   9               0               0     0.050847456             0.2      0.39506173               1
  10               0               0       0.7777778               1               1               1
  11               0               0      0.11111111      0.55263156       0.7586207               1
  12               0               0     0.055555556             0.2           0.375               1
  13               0               0            0.75               1               1               1
  14               0               0      0.09090909      0.53571427       0.7692308               1
  15               0               0               0             0.2      0.33333334               1
  16               0               0             0.8               1               1               1

          sample_     site_
   Obs       id       type     _TYPE_     _FREQ_             mean              sd             min              p5             p10

    17    1Gy_0U_2     CHG        0       4639522    0.1411946616    0.2649249847               0               0               0
    18    1Gy_0U_2     CHH        0      24415887     0.111977697    0.2315206922               0               0               0


   Obs             p25          median             p75             p90             p95             max

    17               0               0           0.125      0.61904764             0.8               1
    18               0               0             0.1      0.41666666       0.6666667               1



                       site_
 Obs     sample_id     type     _TYPE_     _FREQ_            mean              sd             min              p5             p10

   1    01Gy_0U_1       GC         0      5803543    0.0991440156    0.2033247259               0               0               0
   2    01Gy_0U_2       GC         0      5803543    0.1032281087    0.2169142545               0               0               0
   3    01Gy_100U_1     GC         0      5988381    0.4407246306    0.2060837635               0    0.1534358345    0.2148101683
   4    01Gy_100U_2     GC         0      5988381    0.4429558975    0.1822222517               0    0.1871002991    0.2361860047
   5    0Gy_0U_1        GC         0      4996963    0.1051640283    0.2156169022               0               0               0
   6    0Gy_0U_2        GC         0      4996963    0.1142699769    0.2416084869               0               0               0
   7    0Gy_100U_1      GC         0      4712204    0.2921458301    0.2186737622               0               0    0.0765855095
   8    0Gy_100U_2      GC         0      4712204    0.2973473216    0.2323860545               0               0               0
   9    1Gy_0U_1        GC         0      4539526    0.1194853053    0.2579054283               0               0               0
  10    1Gy_0U_2        GC         0      4539526    0.1081512735    0.2140737141               0               0               0
  11    1Gy_100U_1      GC         0      4703751    0.3451559751    0.2131640271               0    0.0598921328    0.1282298642
  12    1Gy_100U_2      GC         0      4703751    0.3405215299    0.2248408372               0               0    0.1026143932


 Obs             p25          median             p75             p90             p95             max

   1               0               0    0.0907304166    0.2690942227    0.6605040011    0.9687392017
   2               0               0     0.093934985    0.2850608606    0.6888970798    1.0333456197
   3    0.3068716689    0.4156126069    0.5370254207    0.7080807376    0.8587700193    1.0621211064
   4    0.3191966867    0.4165516762    0.5395734417    0.6937372822    0.8178797433    0.9712321951
   5               0               0    0.0949283298    0.3497359519     0.733537094    0.9492832981
   6               0               0    0.0956450974    0.3892153942    0.8126475263    1.0564417842
   7    0.1472069181    0.2440233547     0.374931841    0.5827980053    0.8152315628    0.9941848326
   8    0.1537699595    0.2498761842     0.381881293    0.6122769712    0.8382363296    1.0058835955
   9               0               0    0.0967176853      0.42786649    0.8596947068    1.1053217659
  10               0               0    0.0984932474    0.3873348126    0.7304027895    0.9130034869
  11    0.2022423665     0.302835399    0.4488045248    0.6472425053    0.8215001028    0.9708637579
  12    0.1881263875      0.30326886    0.4351248861    0.6344240879    0.8490087059    1.0309391428




                         site_
     Obs    sample_id    type     _TYPE_     _FREQ_            mean              sd             min              p5             p10

       1    01Gy_0U_1     CG         0       974981    0.2847499836    0.3735625185               0               0               0
       2    01Gy_0U_1     CHG        0      1051587    0.1193478135    0.1995983176               0               0               0
       3    01Gy_0U_1     CHH        0      1176944    0.0618824614     0.072158687               0    0.0127664108    0.0179166375
       4    01Gy_0U_2     CG         0       974578    0.2820394599    0.3719611541               0               0               0
       5    01Gy_0U_2     CHG        0      1051409    0.1167809099    0.1978760272               0               0               0
       6    01Gy_0U_2     CHH        0      1176587    0.0590930636    0.0702407687               0     0.011901802    0.0166956119
       7    0Gy_0U_1      CG         0       975036    0.2832119932    0.3812016649               0               0               0
       8    0Gy_0U_1      CHG        0      1048086     0.119640449    0.2143932801               0               0               0
       9    0Gy_0U_1      CHH        0      1177296    0.0617308069    0.0886554828               0    0.0041911765    0.0086206897
      10    0Gy_0U_2      CG         0       975661    0.2838027913    0.3787429271               0               0               0
      11    0Gy_0U_2      CHG        0      1050321    0.1205511831    0.2119137731               0               0               0
      12    0Gy_0U_2      CHH        0      1177801     0.062818967    0.0869716975               0    0.0053763442    0.0092933008
      13    1Gy_0U_1      CG         0       974858    0.2791081121     0.379289036               0               0               0
      14    1Gy_0U_1      CHG        0      1048835    0.1152138718    0.2112418859               0               0               0
      15    1Gy_0U_1      CHH        0      1177001    0.0564550169    0.0832035205               0    0.0022321429    0.0067340067
      16    1Gy_0U_2      CG         0       975628    0.2935107045    0.3795521679               0               0               0


     Obs             p25          median             p75             p90             p95             max

       1        0.015625    0.0439749225    0.7194527025       0.9081755     0.938690475               1
       2     0.013888889    0.0335097011    0.0833333333    0.4760299131     0.621190486               1
       3    0.0268473008    0.0383454886    0.0599102531    0.1442640716    0.2004794887               1
       4    0.0143073597    0.0420634941       0.7070707         0.90625      0.93877315               1
       5          0.0125    0.0323695271            0.08      0.46753917      0.61538464               1
       6    0.0252287419    0.0362741074    0.0563897508    0.1392266716    0.1952711663               1
       7               0    0.0324074076    0.7289858108    0.9169146875       0.9462963               1
       8               0    0.0222222233    0.0777777783    0.5069744925    0.6547052286               1
       9    0.0174305165    0.0303884721    0.0598639461    0.1654437123    0.2332796787               1
      10               0    0.0357142875    0.7191544525      0.91666665    0.9477195443               1
      11               0    0.0243819151     0.083333335             0.5    0.6496071278               1
      12    0.0179012348    0.0315555556    0.0642192897    0.1647294173    0.2309615388               1
      13               0         0.03125     0.704283644      0.91696952     0.951960795               1
      14               0            0.02            0.08    0.4957000533     0.646311868               1
      15      0.01472532     0.027777778           0.055    0.1509896767    0.2189631406               1
      16     0.004347826    0.0403689106       0.7357143    0.9181818233    0.9473593079               1


      s      s
      a      i
      m      t
      p      e  _    _                                                                                        m
      l      _  T    F                                                                                        e
      e      t  Y    R          m                                                                             d
 O    _      y  P    E          e                         m                         p            p            i            p            p            p            m
 b    i      p  E    Q          a            s            i            p            1            2            a            7            9            9            a
 s    d      e  _    _          n            d            n            5            0            5            n            5            0            5            x

17 1Gy_0U_2 CHG 0 1050971 0.1383379071 0.2290246345            0            0            0            0        0.025        0.125  0.547058825 0.6833333333            1
18 1Gy_0U_2 CHH 0 1178187 0.1062766794 0.1664929675            0  0.005714286 0.0096172142 0.0183666375 0.0333333335 0.1123282279 0.3178571467 0.5017451931            1



                                            Correlation non-GC endogenous 0Gy_0U_1 0Gy_0U_2 (100X reads)                     11:27 Thursday, March 24, 2022 486

                                site_
            Obs    sample_id    type     _TYPE_     _FREQ_            mean              sd             min              p5             p10

              1    01Gy_0U_1     CG         0       968287    0.2848304358    0.3745948596               0               0               0
              2    01Gy_0U_1     CHG        0      1036907     0.121684784    0.2058556512               0               0               0
              3    01Gy_0U_1     CHH        0      1176508    0.0628241135     0.074415954               0    0.0107339634    0.0162627554
              4    01Gy_0U_2     CG         0       968108    0.2819384405     0.372971906               0               0               0
              5    01Gy_0U_2     CHG        0      1036971    0.1188981249    0.2040475073               0               0               0
              6    01Gy_0U_2     CHH        0      1176179     0.059841972    0.0723835587               0    0.0099256079    0.0150793657
              7    0Gy_0U_1      CG         0       964938    0.2838342357    0.3827349527               0               0               0
              8    0Gy_0U_1      CHG        0      1029572    0.1223908206    0.2212712432               0               0               0
              9    0Gy_0U_1      CHH        0      1176493    0.0629318368    0.0914401436               0               0    0.0069204153
             10    0Gy_0U_2      CG         0       966857    0.2842131657    0.3802647832               0               0               0
             11    0Gy_0U_2      CHG        0      1033164    0.1230147112    0.2186426641               0               0               0
             12    0Gy_0U_2      CHH        0      1177202    0.0639419968    0.0897363358               0    0.0034782608    0.0077191197
             13    1Gy_0U_1      CG         0       964986    0.2794901834    0.3810141281               0               0               0
             14    1Gy_0U_1      CHG        0      1031017     0.117596203    0.2182786488               0               0               0
             15    1Gy_0U_1      CHH        0      1176332    0.0573581118    0.0858622841               0               0    0.0050029868
             16    1Gy_0U_2      CG         0       967410    0.2952775743    0.3819343696               0               0               0


            Obs             p25          median             p75             p90             p95             max

              1     0.013333334    0.0446428583       0.7199074    0.9095238183      0.94285715               1
              2    0.0104166667    0.0333333333    0.0866935475     0.493467648      0.64033969               1
              3    0.0260629392    0.0386904771    0.0621164345     0.148135103    0.2076849712               1
              4    0.0119047625    0.0426136363       0.7058824      0.90789473      0.94285715               1
              5    0.0085470089         0.03125     0.083333335      0.48333334       0.6341912               1
              6    0.0243708235    0.0364441629    0.0582750588    0.1427690729    0.2018988667               1
              7               0    0.0317686576      0.73243245      0.91953441      0.95108225               1
              8               0    0.0208333333    0.0833333333           0.525     0.672630866               1
              9    0.0162037039    0.0304814801          0.0625    0.1704752599    0.2425943837               1
             10               0    0.0357142875    0.7208423186       0.9189189      0.95283664               1
             11               0    0.0227272725      0.08888889     0.515378785       0.6666667               1
             12     0.016784937    0.0317307692    0.0666666667    0.1694444467    0.2398303776               1
             13               0     0.029670331     0.705086585      0.92105263      0.95833335               1
             14               0    0.0167032971    0.0833333333      0.50176055       0.6666667               1
             15    0.0134920636     0.027763254    0.0573135839     0.155250256    0.2270601988               1
             16               0            0.04      0.74404765    0.9213564133       0.9528302               1

                                            Correlation non-GC endogenous 0Gy_0U_1 0Gy_0U_2 (100X reads)                     11:27 Thursday, March 24, 2022 487

     s      s
     a      i
     m      t
     p      e  _    _                                                                                        m
     l      _  T    F                                                                                        e
     e      t  Y    R          m                                                                             d
O    _      y  P    E          e                         m                         p            p            i            p            p            p            m
b    i      p  E    Q          a            s            i            p            1            2            a            7            9            9            a
s    d      e  _    _          n            d            n            5            0            5            n            5            0            5            x

7 1Gy_0U_2 CHG 0 1034838 0.1428934042 0.2387168574            0            0            0            0  0.023809525  0.133333335   0.57256633    0.7094709            1
8 1Gy_0U_2 CHH 0 1177875 0.1091695337 0.1719676262            0 0.0039345713 0.0081453638  0.017434079 0.0338666247 0.1160291712 0.3297325115 0.5177057435            1



                             site_
       Obs     sample_id     type     _TYPE_     _FREQ_            mean              sd             min              p5             p10

         1    01Gy_0U_1       GC         0      1061479    0.0906234225    0.1325174059               0               0               0
         2    01Gy_0U_2       GC         0      1061479    0.0939635535    0.1405346026               0               0               0
         3    01Gy_100U_1     GC         0      1062837    0.4262890039    0.1294973953               0    0.2320954927    0.2732116303
         4    01Gy_100U_2     GC         0      1062837    0.4290702691    0.1163155554               0    0.2547321223    0.2914860391
         5    0Gy_0U_1        GC         0      1037206    0.0899761619    0.1430255769               0               0               0
         6    0Gy_0U_2        GC         0      1037206    0.0969886051    0.1572070711               0               0               0
         7    0Gy_100U_1      GC         0      1015867    0.2733895773    0.1509753843               0    0.0666997169    0.1095957119
         8    0Gy_100U_2      GC         0      1015867    0.2793111967    0.1528734121               0    0.0716027424    0.1138280254
         9    1Gy_0U_1        GC         0      1031730    0.0987729486    0.1665384078               0               0               0
        10    1Gy_0U_2        GC         0      1031730    0.0910560148      0.14543672               0               0               0
        11    1Gy_100U_1      GC         0      1026996    0.3277059984    0.1450395948               0    0.1246417497    0.1651829449
        12    1Gy_100U_2      GC         0      1026996    0.3242416632    0.1516098259               0     0.111198582    0.1544242727


       Obs             p25          median             p75             p90             p95             max

         1    0.0184937006    0.0390625108    0.0908193002     0.272382334    0.3914333915    0.9687392017
         2    0.0182171337    0.0392467476    0.0920444652    0.2863176407    0.4133382479    1.0333456197
         3    0.3404203848    0.4173387677    0.5019256879    0.5891945271     0.650211586    1.0463225648
         4    0.3519475968    0.4211850244     0.497095391    0.5753921134    0.6295334276    0.9712321951
         5    0.0052185883     0.030814822    0.0952458133    0.2920871686    0.4120666596    0.9492832981
         6    0.0063763398    0.0318816991    0.0978467711    0.3174131333    0.4529342086    1.0564417842
         7    0.1727590808    0.2507380635    0.3499363849    0.4675545271    0.5512070498    0.9941848326
         8    0.1784183234    0.2553358442    0.3557494117    0.4773516163    0.5658095225    1.0058835955
         9               0    0.0292682085    0.1015535695     0.333382917    0.4777466236    1.1053217659
        10    0.0043381232    0.0295652012    0.0977591598    0.3016391123     0.425534418    0.9130034869
        11    0.2309364031    0.3099311036    0.4066483592    0.5152759154    0.5894313055    0.9708637579
        12    0.2246270862    0.3069313356    0.4054433809    0.5173504373    0.5977066271    1.0309391428


*/



/*

(3) Bisulfite conversion ratio: HCHH 100U, CHH 100U, HCHH 0U, CHH 0U on Pt sites only 
*/



%macro importBED(trt,unit,rep,siteType);

     data WORK.input_data    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 "!HOME/concannon/DTRA/radiation_bed_files/arabidopsis/by_rep/at_rad_&siteType._all_&trt._&unit._&rep..bed"
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat chr $2. ;
        informat start_pos best32. ;
        informat stop_pos best32. ;
        informat perc_methyl best32. ;
        informat total_C best32. ;
        informat methyl_C best32. ;
        informat VAR7 $1. ;
        informat VAR8 best32. ;
        informat VAR9 $1. ;
        informat VAR10 best32. ;
        informat VAR11 best32. ;
        informat VAR12 $1. ;
        informat VAR13 best32. ;
        informat VAR14 best32. ;
        format chr $2. ;
        format start_pos best12. ;
        format stop_pos best12. ;
        format perc_methyl best12. ;
        format total_C best12. ;
        format methyl_C best12. ;
        format VAR7 $1. ;
        format VAR8 best12. ;
        format VAR9 $1. ;
        format VAR10 best12. ;
        format VAR11 best12. ;
        format VAR12 $1. ;
        format VAR13 best12. ;
        format VAR14 best12. ;
     input
               chr $
               start_pos
               stop_pos
               perc_methyl
               total_C
               methyl_C
               VAR7 $
               VAR8
               VAR9 $
               VAR10
               VAR11
               VAR12 $
               VAR13
               VAR14
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;



data &siteType._&trt._&unit._&rep.;
  length condition $15.;
  length dose $4.;
  length units $4.;
  length site_type $3.;
  length rep $4.;
  set input_data;
  dose="&trt.";
  rep="&rep.";
  units="&unit.";
  condition=catx("_",dose,units);
  site_type=upcase("&siteType.");
  keep chr  start_pos stop_pos perc_methyl rep total_C methyl_C  condition dose units site_type;
run;

%mend;

%importBED(0Gy,100U,1,chh);
%importBED(0Gy,100U,2,chh);
%importBED(01Gy,100U,1,chh);
%importBED(01Gy,100U,2,chh);
%importBED(1Gy,100U,1,chh);
%importBED(1Gy,100U,2,chh);


data chh_100u;
   set chh_: ;
    where chr="Pt";
run;

data gc_sites;
   set arabmap.methylation_data_gc;
   where chr="Pt";
   keep chr start_pos stop_pos;
run;

proc sort data=gc_sites nodup;
  by chr start_pos stop_pos;
proc sort data=chh_100u;
  by chr start_pos stop_pos;
run;

data hchh_100u;
  merge chh_100u (in=in1) gc_sites (in=in2);
  by chr start_pos stop_pos;
  if in2 then delete;
run;


proc sort data=chh_100u;
   by dose units rep;
proc means data=chh_100u noprint;
   by dose units rep;
   var perc_methyl;
   output out=mean_meth_CHH_Pt_by_rep mean=mean_methyl_all stddev=sd_methyl_all;
run;

proc print data=mean_meth_CHH_Pt_by_rep;
run;

proc sort data=hchh_100u;
   by dose units rep;
proc means data=hchh_100u noprint;
   by dose units rep;
   var perc_methyl;
   output out=mean_meth_HCHH_Pt_by_rep mean=mean_methyl_all stddev=sd_methyl_all;
run;

proc print data=mean_meth_HCHH_Pt_by_rep;
run;



/*
                                              mean_methyl_      sd_methyl_
  dose    units    rep    _TYPE_    _FREQ_        all                  all

  01Gy    100U      1        0       26743    0.1603825312    0.2849672177
  01Gy    100U      2        0       26725    0.1667935361    0.2926086694
  0Gy     100U      1        0       28119     0.126473861    0.2482502224
  0Gy     100U      2        0       27097    0.1190796803    0.2300784374
  1Gy     100U      1        0       27593    0.1196787223    0.2187123496
  1Gy     100U      2        0       27501     0.108110657    0.2087951266



                                                    mean_methyl_      sd_methyl_
 Obs    dose    units    rep    _TYPE_    _FREQ_        all                  all

  1     01Gy    100U      1        0       22782    0.0452758381    0.0566418001
  2     01Gy    100U      2        0       22762    0.0493652411    0.0689724972
  3     0Gy     100U      1        0       23887    0.0284024473    0.0484747028
  4     0Gy     100U      2        0       23067    0.0287800745    0.0451179667
  5     1Gy     100U      1        0       23523    0.0350333003    0.0513646279
  6     1Gy     100U      2        0       23436    0.0276555258    0.0432934304


*/





data chh_100u;
   set arabmap.methylation_data_cg_chg_chh ;
    where chr="Pt" and site_type="CHH" and units="0U";
run;

data gc_sites;
   set arabmap.methylation_data_gc;
   where chr="Pt";
   keep chr start_pos stop_pos;
run;

proc sort data=gc_sites nodup;
  by chr start_pos stop_pos;
proc sort data=chh_100u;
  by chr start_pos stop_pos;
run;

data hchh_100u;
  merge chh_100u (in=in1) gc_sites (in=in2);
  by chr start_pos stop_pos;
  if in2 then delete;
run;


proc sort data=chh_100u;
   by treatment units rep;
proc means data=chh_100u noprint;
   by treatment units rep;
   var perc_methyl;
   output out=mean_meth_CHH_Pt_by_rep mean=mean_methyl_all stddev=sd_methyl_all;
run;

proc print data=mean_meth_CHH_Pt_by_rep;
run;

proc sort data=hchh_100u;
   by treatment units rep;
proc means data=hchh_100u noprint;
   by treatment units rep;
   var perc_methyl;
   output out=mean_meth_HCHH_Pt_by_rep mean=mean_methyl_all stddev=sd_methyl_all;
run;

proc print data=mean_meth_HCHH_Pt_by_rep;
run;


/*

                                                                mean_methyl_      sd_methyl_
bs    treatment    units             rep    _TYPE_    _FREQ_        all                  all

1       01Gy        0U                 1       0       26617    0.0340701666    0.0514181153
2       01Gy        0U                 2       0       26173    0.0319882793    0.0443313611
3       0Gy         0U                 1       0       27822    0.0271857812    0.0480980497
4       0Gy         0U                 2       0       27406    0.0292746323    0.0430725983
5       1Gy         0U                 1       0       27075    0.0244003449    0.0397991234
6       1Gy         0U                 2       0       26838    0.0467796456    0.0948773482


                                                                 mean_methyl_      sd_methyl_
Obs    treatment    units             rep    _TYPE_    _FREQ_        all                  all

 1       01Gy        0U                 1       0       22716    0.0340420658    0.0505953282
 2       01Gy        0U                 2       0       22336    0.0319650193    0.0449660834
 3       0Gy         0U                 1       0       23699    0.0272799752    0.0474223711
 4       0Gy         0U                 2       0       23347    0.0292460636    0.0404076607
 5       1Gy         0U                 1       0       23112    0.0245740443    0.0410213716
 6       1Gy         0U                 2       0       22874    0.0490378368    0.0985669606



*/




/*

(4) Gene and TE methylation/accessibilty plots (+/- 2kb)
*/




/*

(5) DE vs non DE gene methylation/accessibility plots (+/2 2kb)
    -> do by time 
*/

proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/Arabidopsis_thaliana.TAIR10.37.gtf"
  out=gtf dbms=tab replace;
  guessingrows=max;
  getnames=no;
run;


data gene;
  set gtf;
  where VAR3 = "gene";
  length gene_id $15.;
  gene_id=compress(tranwrd(tranwrd(scan(VAR9,2," "), ";", ""), '"', ''));
  keep gene_id VAR1 VAR4 VAR5 VAR7; 
  rename VAR7=strand VAR1=chr VAR4=gene_start VAR5=gene_stop;
run;


/* prep methylation and accessibility counts */



data meth_data;
  set arabMAP.methylation_data_cg_chg_chh;
  keep site_type chr stop_pos treatment rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by site_type chr stop_pos treatment  rep ;
proc means data=meth_data noprint;
  by site_type chr stop_pos treatment   ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;

data meth_data3;
  set meth_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=meth_data3 out=meth_sbys10;
  where total_C >= 10;
  by site_type chr stop_pos;
  id treatment ;
  var perc_methyl2;
run;

data meth_sbys10_2;
  set meth_sbys10;
  if _01Gy ne . and _0Gy ne . then _01Gy_common=_01Gy; else _01Gy_common=.;
  if _1Gy ne . and _0Gy ne . then _1Gy_common=_1Gy; else _1Gy_common=.;
  if (_1Gy ne . or _01Gy ne .) and _0Gy ne . then _0Gy_common=_0Gy; else _0Gy_common=.;

  if _01Gy_common ne . and _1Gy_common ne .  then do;
    _01Gy_common_all = _01Gy_common;
    _1Gy_common_all = _1Gy_common;
    _0Gy_common_all = _0Gy_common;
    end;
  else do;
   _01Gy_common_all = .;
   _1Gy_common_all = .;
   _0Gy_common_all = . ;
    end;
  rename stop_pos=pos;
run;

proc sort data=gene;
  by chr gene_id;
run;

data gene_2;
  set gene;
  by chr gene_id;
  gene_length=gene_stop-gene_start;
  do pos = gene_start to gene_stop;
  output;
  end;
run;

proc sort data=gene_2;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data meth_sbys10_gene;
  merge meth_sbys10_2 (in=in1) gene_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data meth_sbys10_gene2;
  set meth_sbys10_gene;
  if strand="+" then pos_bin=int((pos-gene_start)/gene_length * 2000) + 2000 ;
  else if strand="-" then pos_bin=int((gene_stop-pos)/gene_length * 2000) + 2000 ;
run;

data promoter;
  set gene;
  if strand="+" then do;
    promoter_start=gene_start-2000;
    promoter_stop=gene_start-1;
    end;
  else do;
    promoter_start=gene_stop+2000;
    promoter_stop=gene_stop+1;
    end;
run;

data downstream;
  set gene;
  if strand="-" then do;
    downstream_start=gene_start-1;
    downstream_stop=gene_start-2000;
    end;
  else do;
    downstream_start=gene_stop+1;
    downstream_stop=gene_stop+2000;
    end;
run;

proc sort data=promoter;
  by chr gene_id ;
proc sort data=downstream;
  by chr gene_id ;
run;


data promoter_2;
  set promoter;
  by chr gene_id;
  if strand="+" then do pos = promoter_start to promoter_stop;
  output;
  end;
  else if strand="-" then do pos = promoter_stop to promoter_start;
  output;
  end;
run;


data downstream_2;
  set downstream;
  by chr gene_id;
  if strand="+" then do pos = downstream_start to downstream_stop;
  output;
  end;
  else if strand="-" then do pos = downstream_stop to downstream_start;
  output;
  end;
run;



proc sort data=promoter_2;
  by chr pos;
proc sort data=downstream_2;
  by chr pos;
run;

data meth_sbys10_promoter;
  merge meth_sbys10_2 (in=in1) promoter_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_downstream;
  merge meth_sbys10_2 (in=in1) downstream_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_promoter2;
  set meth_sbys10_promoter;
  if strand="+" then pos_bin=(pos-promoter_start);
  else if strand="-" then pos_bin=(promoter_start-pos) ;
run;

data meth_sbys10_downstream2;
  set meth_sbys10_downstream;
  if strand="+" then pos_bin=(pos-downstream_start) + 4000;
  else if strand="-" then pos_bin=(downstream_start-pos) + 4000;
run;



data meth_sbys10_genestack;
   set meth_sbys10_gene2 (in=in1) meth_sbys10_promoter2 (in=in2) meth_sbys10_downstream2 (in=in3);
run;

proc sort data=meth_sbys10_genestack;
  by gene_id chr pos_bin;
run;



data de_results;
  set arabRNA.arab_results_by_gene;
  keep gene_id flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
  flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
  flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
  mean_cpm_: 
  fdr_: 
  flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_3h_fdr05
  flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_72h_fdr05
  flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_3h_fdr05
  flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_72h_fdr05 ;
run;

proc contents data=de_results;run;quit;


data de_results2;
  set de_results;
  log2fc_01gy_Mock_1h = log2(mean_cpm_01gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_01gy_Mock_3h = log2(mean_cpm_01gy_3h ) - log2(mean_cpm_Mock_3h );
  log2fc_01gy_Mock_24h = log2(mean_cpm_01gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_01gy_Mock_72h = log2(mean_cpm_01gy_72h ) - log2(mean_cpm_Mock_72h );
  log2fc_1gy_Mock_1h = log2(mean_cpm_1gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_1gy_Mock_3h = log2(mean_cpm_1gy_3h) - log2(mean_cpm_Mock_3h );
  log2fc_1gy_Mock_24h = log2(mean_cpm_1gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_1gy_Mock_72h = log2(mean_cpm_1gy_72h ) - log2(mean_cpm_Mock_72h );
run;

data up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
     dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1
     up_1_1_fc1  up_1_3_fc1  up_1_24_fc1  up_1_72_fc1
     dn_1_1_fc1  dn_1_3_fc1  dn_1_24_fc1  dn_1_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3_fc1;
if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24_fc1;
if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3_fc1;
if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24_fc1;
if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;
keep gene_id;
run;

data gene_on_01_1 gene_on_01_3 gene_on_01_24 gene_on_01_72
     gene_on_1_1 gene_on_1_3 gene_on_1_24 gene_on_1_72;
  set de_results2;
  if flag_01gy_1hr_on = 1 or flag_Mock_1hr_on=1 then output gene_on_01_1;
  if flag_01gy_3hr_on = 1 or flag_Mock_3hr_on=1 then output gene_on_01_3;
  if flag_01gy_24hr_on = 1 or flag_Mock_24hr_on=1 then output gene_on_01_24;
  if flag_01gy_72hr_on = 1 or flag_Mock_72hr_on=1 then output gene_on_01_72;
  if flag_1gy_1hr_on = 1 or flag_Mock_1hr_on=1 then output gene_on_1_1;
  if flag_1gy_3hr_on = 1 or flag_Mock_3hr_on=1 then output gene_on_1_3;
  if flag_1gy_24hr_on = 1 or flag_Mock_24hr_on=1 then output gene_on_1_24;
  if flag_1gy_72hr_on = 1 or flag_Mock_72hr_on=1 then output gene_on_1_72;
  keep gene_id;
run;


%macro exportGBDE(siteType, geneList, upList, downList);

proc sort data=&geneList.;
  by gene_id;
proc sort data=&upList.;
  by gene_id;
proc sort data=&downList.;
  by gene_id;
run;

data genes2plot;
  merge &geneList. (in=in1) &upList. (in=in2) &downList.  (in=in3);
  by gene_id;
  length DE_type $12.;
  if in2 then DE_type = "Upreg";
  else if in3 then DE_type = "Downreg";
  else DE_type = "not_DE";
  if in1 then output;
run;



proc sort data=genes2plot;
   by gene_id;
proc sort data=meth_sbys10_genestack;
  by gene_id;
run;

data meth_data;
  merge genes2plot (in=in1) meth_sbys10_genestack (in=in2);
  by gene_id;
  if in1 and in2;
run;

data meth_data2;
  set meth_data;
  where site_type="&siteType.";
  grouped_pos=int(pos_bin/10)*10;
  diff_meth_10cgy=_01Gy_common - _0Gy_common ;
  diff_meth_100cgy=_1Gy_common - _0Gy_common;
run;


proc sort data=meth_data2;
  by DE_type pos_bin grouped_pos;
run;


proc means data=meth_data2 noprint;
  by DE_type pos_bin grouped_pos  ;
  var  diff_meth_10cgy diff_meth_100cgy
   ;
  output out=meth_data3
  mean(diff_meth_10cgy)=mean_diff_meth_10cgy
  mean(diff_meth_100cgy)=mean_diff_meth_100cgy;
run;


proc sort data=meth_data3 ;
  by DE_type grouped_pos;
run;


proc means data=meth_data3 noprint;
  by DE_type grouped_pos  ;
  var  mean_diff_meth_10cgy mean_diff_meth_100cgy  
       ;
  output out=meth_data4
  mean(mean_diff_meth_10cgy)=mean_diff_meth_10cgy
  mean(mean_diff_meth_100cgy)=mean_diff_meth_100cgy;
run;

proc sort data=meth_data4;
  by grouped_pos DE_type;
proc transpose data=meth_data4 out=meth_data_10cgy;
  by grouped_pos ;
  id DE_type;
  var mean_diff_meth_10cgy;
run;

proc transpose data=meth_data4 out=meth_data_100cgy;
  by grouped_pos ;
  id DE_type;
  var mean_diff_meth_100cgy;
run;



data meth_data_export3;
  set meth_data_10cgy;
  drop _NAME_;
  rename grouped_pos=pos;
run;

data meth_data_export4;
  set meth_data_100cgy;
  drop _NAME_ ;
   rename grouped_pos=pos;
run;


proc export data=meth_data_export3
   outfile="!HOME/concannon/DTRA/genebody_plot_2kb_&geneList._&siteType._10cGy_diff_meth.csv"
   dbms=csv replace;
run;
proc export data=meth_data_export4
   outfile="!HOME/concannon/DTRA/genebody_plot_2kb_&geneList._&siteType._100cGy_diff_meth.csv"
   dbms=csv replace;
run;

%mend;





%exportGBDE(CG, gene_on_01_1, up_01_1_fc1, dn_01_1_fc1);    %exportGBDE(CHG, gene_on_01_1, up_01_1_fc1, dn_01_1_fc1); %exportGBDE(CHH, gene_on_01_1, up_01_1_fc1, dn_01_1_fc1);
%exportGBDE(CG, gene_on_01_3, up_01_3_fc1, dn_01_3_fc1);    %exportGBDE(CHG, gene_on_01_3, up_01_3_fc1, dn_01_3_fc1);    %exportGBDE(CHH, gene_on_01_3, up_01_3_fc1, dn_01_3_fc1);
%exportGBDE(CG, gene_on_01_24, up_01_24_fc1, dn_01_24_fc1); %exportGBDE(CHG, gene_on_01_24, up_01_24_fc1, dn_01_24_fc1); %exportGBDE(CHH, gene_on_01_24, up_01_24_fc1, dn_01_24_fc1);
%exportGBDE(CG, gene_on_01_72, up_01_72_fc1, dn_01_72_fc1); %exportGBDE(CHG, gene_on_01_72, up_01_72_fc1, dn_01_72_fc1); %exportGBDE(CHH, gene_on_01_72, up_01_72_fc1, dn_01_72_fc1);

%exportGBDE(CG, gene_on_1_1, up_1_1_fc1, dn_1_1_fc1);    %exportGBDE(CHG, gene_on_1_1, up_1_1_fc1, dn_1_1_fc1); %exportGBDE(CHH, gene_on_1_1, up_1_1_fc1, dn_1_1_fc1);
%exportGBDE(CG, gene_on_1_3, up_1_3_fc1, dn_1_3_fc1);    %exportGBDE(CHG, gene_on_1_3, up_1_3_fc1, dn_1_3_fc1);    %exportGBDE(CHH, gene_on_1_3, up_1_3_fc1, dn_1_3_fc1);
%exportGBDE(CG, gene_on_1_24, up_1_24_fc1, dn_1_24_fc1); %exportGBDE(CHG, gene_on_1_24, up_1_24_fc1, dn_1_24_fc1); %exportGBDE(CHH, gene_on_1_24, up_1_24_fc1, dn_1_24_fc1);
%exportGBDE(CG, gene_on_1_72, up_1_72_fc1, dn_1_72_fc1); %exportGBDE(CHG, gene_on_1_72, up_1_72_fc1, dn_1_72_fc1); %exportGBDE(CHH, gene_on_1_72, up_1_72_fc1, dn_1_72_fc1);






%macro exportGB(siteType,geneList);

proc sort data=&geneList.;
   by gene_id;
proc sort data=meth_sbys10_genestack;
  by gene_id;
run;

data meth_data;
  merge &geneList. (in=in1) meth_sbys10_genestack (in=in2);
  by gene_id;
  if in1 and in2;
run;

data meth_data2;
  set meth_data;
  where site_type="&siteType.";
  grouped_pos=int(pos_bin/10)*10;
run;


proc sort data=meth_data2;
  by pos_bin grouped_pos;
run;


proc means data=meth_data2 noprint;
  by pos_bin grouped_pos  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
   ;
  output out=meth_data3
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common;
run;


proc sort data=meth_data3 ;
  by grouped_pos;
run;


proc means data=meth_data3 noprint;
  by grouped_pos  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       ;
  output out=meth_data4
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common;
run;


data meth_data_export3;
  set meth_data3;
  drop grouped_pos _TYPE_ _FREQ_;
   rename pos_bin=pos;
run;

data meth_data_export4;
  set meth_data4;
  drop _TYPE_ _FREQ_ mean_01gy_;
   rename grouped_pos=pos;
run;


proc export data=meth_data_export3
   outfile="!HOME/concannon/DTRA/genebody_plot_2kb_&geneList._&siteType._full.csv"
   dbms=csv replace;
run;
proc export data=meth_data_export4
   outfile="!HOME/concannon/DTRA/genebody_plot_2kb_&geneList._&siteType._10bp.csv"
   dbms=csv replace;
run;

%mend;

%exportGB(CG,up_01_1_fc1); %exportGB(CG,up_01_3_fc1); %exportGB(CG,up_01_24_fc1); %exportGB(CG,up_01_72_fc1);
%exportGB(CG,dn_01_1_fc1); %exportGB(CG,dn_01_3_fc1); %exportGB(CG,dn_01_24_fc1); %exportGB(CG,dn_01_72_fc1);
%exportGB(CG,up_1_1_fc1); %exportGB(CG,up_1_3_fc1); %exportGB(CG,up_1_24_fc1); %exportGB(CG,up_1_72_fc1);
%exportGB(CG,dn_1_1_fc1); %exportGB(CG,dn_1_3_fc1); %exportGB(CG,dn_1_24_fc1); %exportGB(CG,dn_1_72_fc1);

   
%exportGB(CHG,up_01_1_fc1); %exportGB(CHG,up_01_3_fc1); %exportGB(CHG,up_01_24_fc1); %exportGB(CHG,up_01_72_fc1);
%exportGB(CHG,dn_01_1_fc1); %exportGB(CHG,dn_01_3_fc1); %exportGB(CHG,dn_01_24_fc1); %exportGB(CHG,dn_01_72_fc1);
%exportGB(CHG,up_1_1_fc1); %exportGB(CHG,up_1_3_fc1); %exportGB(CHG,up_1_24_fc1); %exportGB(CHG,up_1_72_fc1);
%exportGB(CHG,dn_1_1_fc1); %exportGB(CHG,dn_1_3_fc1); %exportGB(CHG,dn_1_24_fc1); %exportGB(CHG,dn_1_72_fc1);

%exportGB(CHH,up_01_1_fc1); %exportGB(CHH,up_01_3_fc1); %exportGB(CHH,up_01_24_fc1); %exportGB(CHH,up_01_72_fc1);
%exportGB(CHH,dn_01_1_fc1); %exportGB(CHH,dn_01_3_fc1); %exportGB(CHH,dn_01_24_fc1); %exportGB(CHH,dn_01_72_fc1);
%exportGB(CHH,up_1_1_fc1); %exportGB(CHH,up_1_3_fc1); %exportGB(CHH,up_1_24_fc1); %exportGB(CHH,up_1_72_fc1);
%exportGB(CHH,dn_1_1_fc1); %exportGB(CHH,dn_1_3_fc1); %exportGB(CHH,dn_1_24_fc1); %exportGB(CHH,dn_1_72_fc1);

   
        
           
           


/* get DEG lists */



proc sort data=mean_diff_tss ;
  by grouped_pos2;
run;


proc means data=mean_diff_tss noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_tss_2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_tss_1;
  set mean_diff_tss;
  drop _TYPE_ _FREQ_;
  keep distance_To_tss mean_: ;
  rename distance_to_tss=pos;
run;

data mean_diff_tss_3;
  set mean_diff_tss_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;


/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, var3, outName);

data export;
  retain pos &var1. &var2. &var3.;
  set &inData.;
  keep pos &var1. &var2. &var3.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/&siteType._methylation_&outName..csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_tss_1, mean_0Gy_common, mean_01Gy_common, mean_1Gy_common, TSS_1kb);
%exportLine( mean_diff_tss_3, mean_0Gy_common, mean_01Gy_common, mean_1Gy_common, TSS_1kb_binned);

%exportLine( mean_diff_tss_1, mean_0Gy_common_all, mean_01Gy_common_all, mean_1Gy_common_all, TSS_1kb_common);
%exportLine( mean_diff_tss_3, mean_0Gy_common_all, mean_01Gy_common_all, mean_1Gy_common_all, TSS_1kb_common_binned);


%mend;

%methDataGen(CG);
%methDataGen(CHG);
%methDataGen(CHH);









