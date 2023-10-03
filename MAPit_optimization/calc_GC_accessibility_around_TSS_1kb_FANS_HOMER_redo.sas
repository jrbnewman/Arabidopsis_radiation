/* For 22C and 4C, show and plot the average level of GC accessibility
   (calculated as 100U - 0U ) +/- 1kb around TSS for all genes 

doing this now for transcripts. I think what I want to do is pull in the transcript TSS sites
and then compare them against the HOMER annotations to see which looks "better"
let Mingqi decide which he likes better
idc
*/

libname cold '!PATCON/DTRA/arabidopsis_wgbs_cold/sas_data';
libname coldloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


* GC methylation data;
data gc_data;
  set coldloc.methylation_data;
  where site_type="GC" ;
  perc_methyl2=perc_methyl*100;
run;

proc sort data=gc_data;
  by chr pos_end condition;
proc transpose data=gc_data out=gc_sbys;
  by chr pos_end;
  id condition;
  var perc_methyl2;
run;


proc transpose data=gc_data out=gc_sbys10;
  where total_C >= 10;
  by chr pos_end;
  id condition;
  var perc_methyl2;
run;

data gc_sbys_2;
  set gc_sbys;
  if _22C_0U ne . and _22C_100U ne . then _22C_100U_0U=_22C_100U-_22C_0U;
  else _22C_100U_0U=.;
  if _4C_0U ne . and _4C_100U ne . then _4C_100U_0U=_4C_100U-_4C_0U;
  else _4C_100U_0U=.;
  if _22C_0U ne . and _22C_100U ne . and _4C_0U ne . and _4C_100U ne . then do;
        _22C_100U_0U_common=_22C_100U - _22C_0U;
        _4C_100U_0U_common=_4C_100U - _4C_0U;
  end;
  rename pos_end=pos;
run;


data gc_sbys10_2;
  set gc_sbys10;
  if _22C_0U ne . and _22C_100U ne . then _22C_100U_0U=_22C_100U-_22C_0U;
  else _22C_100U_0U=.;
  if _4C_0U ne . and _4C_100U ne . then _4C_100U_0U=_4C_100U-_4C_0U;
  else _4C_100U_0U=.;
  if _22C_0U ne . and _22C_100U ne . and _4C_0U ne . and _4C_100U ne . then do;
        _22C_100U_0U_common=_22C_100U - _22C_0U;
        _4C_100U_0U_common=_4C_100U - _4C_0U;
  end;
  rename pos_end=pos;
run;


/* Get TSSs and regions from Ensembl TAIR10/Araport11 GTF annotation */

data annot;
  set cold.homer_annotations_sites;
run;


data site2promoter;
  set annot;
  where siteID ? "GC";
  length gene_id $20.;
  if abs(distance_to_tss) > 1000 then delete;
  if count(nearest_promoterID, "-T1") > 0 then gene_ID=compress(upcase(tranwrd(Nearest_PromoterID,"-T1","")));
  else gene_ID=compress(upcase(scan(Nearest_PromoterID,1,".")));
  keep siteID gene_ID chr start strand distance_to_tss nearest_promoterID;
  rename start=pos nearest_promoterID=transcript_id;
run;

/* Get strand for each transcript so I know what sites should be up/downstream of TSS */


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
  keep gene_id VAR7; 
  rename VAR7=strand;
run;

proc sort data=gene nodup;
  by gene_id;
proc sort data=site2promoter;
    by gene_id;
run;

data site2promoter2 no_strand no_site;
  merge site2promoter (in=in1) gene (in=in2);
  by gene_id;
  if in1 and in2 then output site2promoter2;
  else if in1 then output no_strand;
  else output no_site;
run;

/* distance to TSS is relative to strand of transcript, so I don't need to flip anything!!! */

proc sort data=site2promoter2 nodup;
  by chr pos;
proc sort data=gc_sbys10_2;
  by chr pos;
run;

data tss_pos_gc_10H;
  merge gc_sbys10_2 (in=in1) site2promoter2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

/* Also going to group every 10 positions */


data tss_pos_gc_10_2H;
  set tss_pos_gc_10H;
  grouped_pos=int(distance_to_tss/10) * 10;
run;

/* For each, calculate the mean and SD */

proc sort data=tss_pos_gc_10_2H;
  by distance_to_tss grouped_pos;
run;


proc means data=tss_pos_gc_10_2H noprint;
  by distance_to_tss grouped_pos;
  var  _22C_100U _4C_100U   _22C_0U _4C_0U _22C_100U_0U _4C_100U_0U _22C_100U_0U_common _4C_100U_0U_common FANS_0p5U FANS_1p5U FANS_25U FANS_5U;
  output out=mean_diff_tss_10X_1000H
  mean(_22C_100U)=mean_22C_100U  stddev(_22C_100U)=sd_22C_100U

  mean(_22C_0U)=mean_22C_0U  stddev(_22C_0U)=sd_22C_0U
  mean(_4C_100U)=mean_4C_100U  stddev(_4C_100U)=sd_4C_100U
  mean(_4C_0U)=mean_4C_0U  stddev(_4C_0U)=sd_4C_0U

  mean(_22C_100U_0U)=mean_22C_100U_0U  stddev(_22C_100U_0U)=sd_22C_100U_0U
  mean(_4C_100U_0U)=mean_4C_100U_0U  stddev(_4C_100U_0U)=sd_4C_100U_0U
  mean(_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(FANS_0p5U)=mean_FANS_0p5U  stddev(FANS_0p5U)=sd_FANS_0p5U
  mean(FANS_1p5U)=mean_FANS_1p5U   stddev(FANS_1p5U)=sd_FANS_1p5U
  mean(FANS_25U)=mean_FANS_25U   stddev(FANS_25U)=sd_FANS_25U
  mean(FANS_5U)=mean_FANS_5U  stddev(FANS_5U)=sd_FANS_5U;
run;


proc sort data=mean_diff_tss_10X_1000H;
  by grouped_pos;
run;


proc means data=mean_diff_tss_10X_1000H noprint;
  by grouped_pos;
  var mean_22C_100U  mean_4C_100U   mean_22C_0U mean_4C_0U mean_22C_100U_0U mean_4C_100U_0U mean_22C_100U_0U_common mean_4C_100U_0U_common mean_FANS_0p5U mean_FANS_1p5U mean_FANS_25U mean_FANS_5U;;
  output out=mean_diff_tss_10X_100H
  mean(mean_22C_100U)=mean_22C_100U stddev(mean_22C_100U)=sd_22C_100U

  mean(mean_22C_0U)=mean_22C_0U stddev(mean_22C_0U)=sd_22C_0U
  mean(mean_4C_100U)=mean_4C_100U stddev(mean_4C_100U)=sd_4C_100U
  mean(mean_4C_0U)=mean_4C_0U stddev(mean_4C_0U)=sd_4C_0U

  mean(mean_22C_100U_0U)=mean_22C_100U_0U  stddev(mean_22C_100U_0U)=sd_22C_100U_0U
  mean(mean_4C_100U_0U)=mean_4C_100U_0U  stddev(mean_4C_100U_0U)=sd_4C_100U_0U
  mean(mean_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(mean_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(mean_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(mean_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(mean_FANS_0p5U)=mean_FANS_0p5U   stddev(mean_FANS_0p5U)=sd_FANS_0p5U
  mean(mean_FANS_1p5U)=mean_FANS_1p5U   stddev(mean_FANS_1p5U)=sd_FANS_1p5U
  mean(mean_FANS_25U)=mean_FANS_25U   stddev(mean_FANS_25U)=sd_FANS_25U
  mean(mean_FANS_5U)=mean_FANS_5U   stddev(mean_FANS_5U)=sd_FANS_5U;
run;


/* T-test for 22C vs 4C */

data model_22c;
  set mean_diff_tss_10X_100H;
  length temperature $3.;
  temperature="22C";
  keep grouped_pos mean_22C_100U_0U temperature;
  rename mean_22C_100U_0U=mean_access;
run;

data model_4c;
  set mean_diff_tss_10X_100H;
  length temperature $3.;
  temperature="4C";
  keep grouped_pos mean_4C_100U_0U temperature;
  rename mean_4C_100U_0U=mean_access;
run;

data model_data;
  set model_22c model_4c;
run;

proc sort data=model_data;
   by temperature grouped_pos ;

proc ttest data=model_data;
   class temperature;
   var mean_access;
run;


proc ttest data=model_data;
   where grouped_pos >= -250 and grouped_pos <= 200;
   class temperature;
   var mean_access;
run;

proc ttest data=mean_diff_tss_10X_100H;
   paired mean_22C_100U_0U*mean_4C_100U_0U;
run;

proc ttest data=mean_diff_tss_10X_100H;
   where grouped_pos >= -250 and grouped_pos <= 200;
   paired mean_22C_100U_0U*mean_4C_100U_0U;
run;


/* Whole region (+/- 1kb), unpaired ttest:


             temperature       N        Mean     Std Dev     Std Err     Minimum     Maximum

           22C             201     42.4977      2.5223      0.1779     40.1485     50.8378
           4C              201     39.9211      2.5155      0.1774     37.0569     48.1499
           Diff (1-2)               2.5766      2.5189      0.2513

    temperature    Method               Mean       95% CL Mean        Std Dev      95% CL Std Dev

    22C                              42.4977     42.1469  42.8485      2.5223      2.2975   2.7963
    4C                               39.9211     39.5713  40.2710      2.5155      2.2913   2.7887
    Diff (1-2)     Pooled             2.5766      2.0826   3.0706      2.5189      2.3558   2.7064
    Diff (1-2)     Satterthwaite      2.5766      2.0826   3.0706

                     Method           Variances        DF    t Value    Pr > |t|

                     Pooled           Equal           400      10.25      <.0001
                     Satterthwaite    Unequal         400      10.25      <.0001

                                        Equality of Variances

                          Method      Num DF    Den DF    F Value    Pr > F

                          Folded F       200       200       1.01    0.9695

zoomed region (+/- 250bb), unpaired ttest:

                                      The SAS System                                         22:13 Sunda

                                   The TTEST Procedure

                                  Variable:  mean_access

      temperature      N        Mean     Std Dev     Std Err     Minimum     Maximum

      22C             46     46.2657      2.8224      0.4161     41.9853     50.8378
      4C              46     43.7162      2.7161      0.4005     39.5157     48.1499
      Diff (1-2)              2.5495      2.7698      0.5775

mperature    Method               Mean       95% CL Mean        Std Dev      95% CL Std Dev

C                              46.2657     45.4275  47.1038      2.8224      2.3410   3.5548
                               43.7162     42.9096  44.5228      2.7161      2.2528   3.4210
ff (1-2)     Pooled             2.5495      1.4021   3.6968      2.7698      2.4175   3.2431
ff (1-2)     Satterthwaite      2.5495      1.4021   3.6968

               Method           Variances        DF    t Value    Pr > |t|

               Pooled           Equal            90       4.41      <.0001
               Satterthwaite    Unequal      89.868       4.41      <.0001

                                  Equality of Variances

                    Method      Num DF    Den DF    F Value    Pr > F

                    Folded F        45        45       1.08    0.7980


whole region (+/- 1kb), paired ttest:


       Difference:  mean_22C_100U_0U - mean_4C_100U_0U

  N        Mean     Std Dev     Std Err     Minimum     Maximum

201      2.5766      0.2828      0.0199      1.8687      3.5965

    Mean       95% CL Mean        Std Dev      95% CL Std Dev

  2.5766      2.5373   2.6159      0.2828      0.2576   0.3135

                    DF    t Value    Pr > |t|

                   200     129.18      <.0001


zoomed region (+/- 250bb), paired ttest:



        Difference:  mean_22C_100U_0U - mean_4C_100U_0U

  N        Mean     Std Dev     Std Err     Minimum     Maximum

 46      2.5495      0.2433      0.0359      2.0047      3.1233

     Mean       95% CL Mean        Std Dev      95% CL Std Dev

   2.5495      2.4772   2.6217      0.2433      0.2018   0.3065

                     DF    t Value    Pr > |t|

                     45      71.07      <.0001
*/


data mean_diff_tss_10X_1000_2H;
  set mean_diff_tss_10X_1000H;
  drop grouped_pos;
  rename distance_to_tss=pos;
run;


data mean_diff_tss_10X_100_2H;
  set mean_diff_tss_10X_100H;
  rename grouped_pos=pos;
run;


data fans_1000;
   set mean_diff_tss_10X_1000_2H;
   keep pos mean_fans_0p5u mean_fans_1p5u mean_fans_5u mean_fans_25u ;
run;

data gc_1000;
   set mean_diff_tss_10X_1000_2H;
   keep pos mean_22c_100U mean_22c_0U pos mean_4c_100U mean_4c_0U;
run;


data fans_100;
   set mean_diff_tss_10X_100_2H;
   keep pos mean_fans_0p5u mean_fans_1p5u mean_fans_5u mean_fans_25u ;
run;

data gc_100;
   set mean_diff_tss_10X_100_2H;
   keep pos mean_22c_100U mean_22c_0U pos mean_4c_100U mean_4c_0U;
run;



/* Export and plot in python */

proc export data=fans_1000
   outfile="!PATCON/DTRA/arabidopsis_wgbs_cold/analysis_output/FANS_only_GC_1kb_TSS_10X_sites_nobin_HOMER.csv"
   dbms=csv
   replace;
run;

proc export data=gc_1000
   outfile="!PATCON/DTRA/arabidopsis_wgbs_cold/analysis_output/GC_22C_4C_100U_0U_1kb_TSS_10X_sites_nobin_HOMER.csv"
   dbms=csv
   replace;
run;




proc export data=fans_100
   outfile="!PATCON/DTRA/arabidopsis_wgbs_cold/analysis_output/FANS_only_GC_1kb_TSS_10X_sites_binned_HOMER.csv"
   dbms=csv
   replace;
run;

proc export data=gc_100
   outfile="!PATCON/DTRA/arabidopsis_wgbs_cold/analysis_output/GC_22C_4C_100U_0U_1kb_TSS_10X_sites_binned_HOMER.csv"
   dbms=csv
   replace;
run;


