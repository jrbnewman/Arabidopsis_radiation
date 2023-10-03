libname wgbs '!PATCON/arabidopsis_wgbs/sas_data';
libname wgbsloc '/TB14/TB14/sandbox/wgbs_sandbox/sas_data';

/* Output filtered bedGraph-like files 

For CG, CHG, CHH sites:
(1) Flag sites <10 coverage in any sample
(2) Create by-sample BEDs, avg-over-condition BEDs, difference BEDs (from perc methyl)
    for each
        -> Flag-filtered and unfiltered
(3) Convert to bigWig

For GC sites:
(1) Flag sites <10 coverage in any sample
(2) Flag sites where 0U > 100U
(3) Create by-sample BEDs, condition average BEDs
    Diffs between 100U and 0U by dose
    Diffs between doses by units
    Diff of diffs
        -> Flag-filtered and unfiltered
(4) Convert all to bigWig

*/

%macro makeBedGraph(site);

data &site.;
  set wgbsloc.methylation_data;
  where units="0U" and site_type="&site.";
  keep sample_id chr pos perc_methyl total_C methyl_C;
run;

proc sort data=&site.;
  by chr pos sample_id;
run;

proc transpose data=&site. out=&site._sbys_perc;
  by chr pos;
  id sample_id;
  var perc_methyl;
run;

proc transpose data=&site. out=&site._sbys_total;
  by chr pos;
  id sample_id;
  var total_C;
run;

proc transpose data=&site. out=&site._sbys_meth;
  by chr pos;
  id sample_id;
  var methyl_C;
run;

data &site._sbys_perc2;
  set &site._sbys_perc;
  keep chr pos _01Gy_0U_1 _01Gy_0U_2
           _0Gy_0U_1 _0Gy_0U_2
           _1Gy_0U_1 _1Gy_0U_2;
  rename _01Gy_0U_1=perc_01Gy_0U_1
         _01Gy_0U_2=perc_01Gy_0U_2
         _0Gy_0U_1=perc_0Gy_0U_1
         _0Gy_0U_2=perc_0Gy_0U_2
         _1Gy_0U_1=perc_1Gy_0U_1
         _1Gy_0U_2=perc_1Gy_0U_2;
run;

data &site._sbys_total2;
  set &site._sbys_total;
  keep chr pos _01Gy_0U_1 _01Gy_0U_2
           _0Gy_0U_1 _0Gy_0U_2
           _1Gy_0U_1 _1Gy_0U_2;
  rename _01Gy_0U_1=total_01Gy_0U_1
         _01Gy_0U_2=total_01Gy_0U_2
         _0Gy_0U_1=total_0Gy_0U_1
         _0Gy_0U_2=total_0Gy_0U_2
         _1Gy_0U_1=total_1Gy_0U_1
         _1Gy_0U_2=total_1Gy_0U_2;
run;

data &site._sbys_meth2;
  set &site._sbys_meth;
  keep chr pos _01Gy_0U_1 _01Gy_0U_2
           _0Gy_0U_1 _0Gy_0U_2
           _1Gy_0U_1 _1Gy_0U_2;
  rename _01Gy_0U_1=meth_01Gy_0U_1
         _01Gy_0U_2=meth_01Gy_0U_2
         _0Gy_0U_1=meth_0Gy_0U_1
         _0Gy_0U_2=meth_0Gy_0U_2
         _1Gy_0U_1=meth_1Gy_0U_1
         _1Gy_0U_2=meth_1Gy_0U_2;
run;

proc sort data=&site._sbys_perc2;
  by chr pos;
proc sort data=&site._sbys_total2;
  by chr pos;
proc sort data=&site._sbys_meth2;
  by chr pos;
run;

data &site._sites;
  merge &site._sbys_perc2 (in=in1)
        &site._sbys_total2 (in=in2)
        &site._sbys_meth2 (in=in3);
  by chr pos;
  if in1 and in2 and in3;
run;

data &site._flag_sites;
  set &site._sites;
  if total_01Gy_0U_1 = . or total_01Gy_0U_2 = .
  or total_0Gy_0U_1 = . or total_0Gy_0U_2 = .
  or total_1Gy_0U_1 = . or total_1Gy_0U_2 = .
  then flag_any_missing=1;
  else flag_any_missing=0;

  if total_01Gy_0U_1 < 10 or total_01Gy_0U_2 < 10
  or total_0Gy_0U_1 < 10 or total_0Gy_0U_2 < 10
  or total_1Gy_0U_1 < 10 or total_1Gy_0U_2 < 10
  then flag_any_total_lt10=1;
  else flag_any_total_lt10=0;
  if meth_01Gy_0U_1 >= 1 or meth_01Gy_0U_2 >= 1
  or meth_0Gy_0U_1 >= 1 or meth_0Gy_0U_2 >= 1
  or meth_1Gy_0U_1 >= 1 or meth_1Gy_0U_2 >= 1
  then flag_any_meth_ge1=1;
  else flag_any_meth_ge1=0;
run;

data &site._calc_diff_mean;
  retain chr pos pos2;
  set &site._flag_sites;
  mean_perc_01Gy=mean(perc_01Gy_0U_1,perc_01Gy_0U_2);
  mean_perc_0Gy=mean(perc_0Gy_0U_1,perc_0Gy_0U_2);
  mean_perc_1Gy=mean(perc_1Gy_0U_1,perc_1Gy_0U_2);
  mean_diff_01Gy_0Gy=mean_perc_01Gy-mean_perc_0Gy;
  mean_diff_1Gy_0Gy=mean_perc_1Gy-mean_perc_0Gy;
  pos2=pos+1;
run;

/* Export bedGraph-like text files (bedGraphs without header) */

data sample1;
  set &site._calc_diff_mean;
  if perc_01Gy_0U_1 ^= .;
  keep chr pos pos2 perc_01Gy_0U_1;
run;

data sample1_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if perc_01Gy_0U_1 ^= .;
  keep chr pos pos2 perc_01Gy_0U_1;
run;


data sample2;
  set &site._calc_diff_mean;
  if perc_01Gy_0U_2 ^= .;
  keep chr pos pos2 perc_01Gy_0U_2;
run;

data sample2_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if perc_01Gy_0U_2 ^= .;
  keep chr pos pos2 perc_01Gy_0U_2;
run;


data sample3;
  set &site._calc_diff_mean;
  if perc_0Gy_0U_1 ^= .;
  keep chr pos pos2 perc_0Gy_0U_1;
run;

data sample3_filter;
  set &site._calc_diff_mean;
  if perc_0Gy_0U_1 ^= .;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  keep chr pos pos2 perc_0Gy_0U_1;
run;


data sample4;
  set &site._calc_diff_mean;
  if perc_0Gy_0U_2 ^= .;
  keep chr pos pos2 perc_0Gy_0U_2;
run;

data sample4_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if perc_0Gy_0U_2 ^= .;
  keep chr pos pos2 perc_0Gy_0U_2;
run;


data sample5;
  set &site._calc_diff_mean;
  if perc_1Gy_0U_1 ^= .;
  keep chr pos pos2 perc_1Gy_0U_1;
run;

data sample5_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if perc_1Gy_0U_1 ^= .;
  keep chr pos pos2 perc_1Gy_0U_1;
run;


data sample6;
  set &site._calc_diff_mean;
  if perc_1Gy_0U_2 ^= .;
  keep chr pos pos2 perc_1Gy_0U_2;
run;

data sample6_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if perc_1Gy_0U_2 ^= .;
  keep chr pos pos2 perc_1Gy_0U_2;
run;

data group1_mean;
  set &site._calc_diff_mean;
  if mean_perc_01Gy ^= .;
  keep chr pos pos2 mean_perc_01Gy;
run;

data group1_mean_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if mean_perc_01Gy ^= .;
  keep chr pos pos2 mean_perc_01Gy;
run;


data group2_mean;
  set &site._calc_diff_mean;
  if mean_perc_0Gy ^= .;
  keep chr pos pos2 mean_perc_0Gy;
run;

data group2_mean_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if mean_perc_0Gy ^= .;
  keep chr pos pos2 mean_perc_0Gy;
run;


data group3_mean;
  set &site._calc_diff_mean;
  if mean_perc_1Gy ^= .;
  keep chr pos pos2 mean_perc_1Gy;
run;

data group3_mean_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if mean_perc_1Gy ^= .;
  keep chr pos pos2 mean_perc_1Gy;
run;

data diff1;
  set &site._calc_diff_mean;
  if mean_diff_01Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_01Gy_0Gy;
run;

data diff1_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if mean_diff_01Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_01Gy_0Gy;
run;

data diff2;
  set &site._calc_diff_mean;
  if mean_diff_1Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_1Gy_0Gy;
run;

data diff2_filter;
  set &site._calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1;
  if mean_diff_1Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_1Gy_0Gy;
run;

proc export data=sample1
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_rep1_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample1_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_rep1_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;


proc export data=sample2
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_rep2_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample2_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_rep2_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;


proc export data=sample3
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._0Gy_0U_rep1_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample3_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._0Gy_0U_rep1_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;


proc export data=sample4
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._0Gy_0U_rep2_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample4_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._0Gy_0U_rep2_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample5
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_rep1_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample5_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_rep1_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;


proc export data=sample6
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_rep2_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=sample6_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_rep2_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=group1_mean
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_avg_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=group1_mean_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_avg_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;


proc export data=group2_mean
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._0Gy_0U_avg_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=group2_mean_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._0Gy_0U_avg_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;


proc export data=group3_mean
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_avg_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=group3_mean_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_avg_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;



proc export data=diff1
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_sub_0Gy_0U_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=diff1_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._01Gy_0U_sub_0Gy_0U_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;




proc export data=diff2
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_sub_0Gy_0U_all.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

proc export data=diff2_filter
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&site._1Gy_0U_sub_0Gy_0U_10reads.bedGraph"
     dbms=tab replace;
     putnames=no;
run;

%mend;


*%makeBedGraph(CG);
*%makeBedGraph(CHG);
*%makeBedGraph(CHH);


/* GCH sites:
   Per sample
   Per plant 100U-0U diff
   Group average
   Per group 100U-0U diff
   01 v 0 0U, 100U
   1 v 0 0U, 100U
   01 v 0, 100U-0U
   1 v 0, 100U-0U
*/




data gc;
  set wgbsloc.methylation_data;
  where site_type="GC";
  keep sample_id chr pos perc_methyl total_C methyl_C;
run;

proc sort data=gc;
  by chr pos sample_id;
run;

proc transpose data=gc out=gch_sbys_perc;
  by chr pos;
  id sample_id;
  var perc_methyl;
run;

proc transpose data=gc out=gch_sbys_total;
  by chr pos;
  id sample_id;
  var total_C;
run;

proc transpose data=gc out=gch_sbys_meth;
  by chr pos;
  id sample_id;
  var methyl_C;
run;

data gch_sbys_perc2;
  set gch_sbys_perc;
  keep chr pos _01Gy_0U_1 _01Gy_0U_2
           _0Gy_0U_1 _0Gy_0U_2
           _1Gy_0U_1 _1Gy_0U_2
           _01Gy_100U_1 _01Gy_100U_2
           _0Gy_100U_1 _0Gy_100U_2
           _1Gy_100U_1 _1Gy_100U_2;
  rename _01Gy_0U_1=perc_01Gy_0U_1
         _01Gy_0U_2=perc_01Gy_0U_2
         _0Gy_0U_1=perc_0Gy_0U_1
         _0Gy_0U_2=perc_0Gy_0U_2
         _1Gy_0U_1=perc_1Gy_0U_1
         _1Gy_0U_2=perc_1Gy_0U_2
         _01Gy_100U_1=perc_01Gy_100U_1
         _01Gy_100U_2=perc_01Gy_100U_2
         _0Gy_100U_1=perc_0Gy_100U_1
         _0Gy_100U_2=perc_0Gy_100U_2
         _1Gy_100U_1=perc_1Gy_100U_1
         _1Gy_100U_2=perc_1Gy_100U_2;
run;

data gch_sbys_total2;
  set gch_sbys_total;
  keep chr pos _01Gy_0U_1 _01Gy_0U_2
           _0Gy_0U_1 _0Gy_0U_2
           _1Gy_0U_1 _1Gy_0U_2
           _01Gy_100U_1 _01Gy_100U_2
           _0Gy_100U_1 _0Gy_100U_2
           _1Gy_100U_1 _1Gy_100U_2;
  rename _01Gy_0U_1=total_01Gy_0U_1
         _01Gy_0U_2=total_01Gy_0U_2
         _0Gy_0U_1=total_0Gy_0U_1
         _0Gy_0U_2=total_0Gy_0U_2
         _1Gy_0U_1=total_1Gy_0U_1
         _1Gy_0U_2=total_1Gy_0U_2
         _01Gy_100U_1=total_01Gy_100U_1
         _01Gy_100U_2=total_01Gy_100U_2
         _0Gy_100U_1=total_0Gy_100U_1
         _0Gy_100U_2=total_0Gy_100U_2
         _1Gy_100U_1=total_1Gy_100U_1
         _1Gy_100U_2=total_1Gy_100U_2;
run;

data gch_sbys_meth2;
  set gch_sbys_meth;
  keep chr pos _01Gy_0U_1 _01Gy_0U_2
           _0Gy_0U_1 _0Gy_0U_2
           _1Gy_0U_1 _1Gy_0U_2
           _01Gy_100U_1 _01Gy_100U_2
           _0Gy_100U_1 _0Gy_100U_2
           _1Gy_100U_1 _1Gy_100U_2;
  rename _01Gy_0U_1=meth_01Gy_0U_1
         _01Gy_0U_2=meth_01Gy_0U_2
         _0Gy_0U_1=meth_0Gy_0U_1
         _0Gy_0U_2=meth_0Gy_0U_2
         _1Gy_0U_1=meth_1Gy_0U_1
         _1Gy_0U_2=meth_1Gy_0U_2
         _01Gy_100U_1=meth_01Gy_100U_1
         _01Gy_100U_2=meth_01Gy_100U_2
         _0Gy_100U_1=meth_0Gy_100U_1
         _0Gy_100U_2=meth_0Gy_100U_2
         _1Gy_100U_1=meth_1Gy_100U_1
         _1Gy_100U_2=meth_1Gy_100U_2;
run;

proc sort data=gch_sbys_perc2;
  by chr pos;
proc sort data=gch_sbys_total2;
  by chr pos;
proc sort data=gch_sbys_meth2;
  by chr pos;
run;

data gch_sites;
  merge gch_sbys_perc2 (in=in1)
        gch_sbys_total2 (in=in2)
        gch_sbys_meth2 (in=in3);
  by chr pos;
  if in1 and in2 and in3;
run;

data GCH_flag_sites;
  set GCH_sites;
  if total_01Gy_0U_1 = . or total_01Gy_0U_2 = .
  or total_0Gy_0U_1 = . or total_0Gy_0U_2 = .
  or total_1Gy_0U_1 = . or total_1Gy_0U_2 = .
  or total_01Gy_100U_1 = . or total_01Gy_100U_2 = .
  or total_0Gy_100U_1 = . or total_0Gy_100U_2 = .
  or total_1Gy_100U_1 = . or total_1Gy_100U_2 = .
  then flag_any_missing=1;
  else flag_any_missing=0;

  if total_01Gy_0U_1 < 10 or total_01Gy_0U_2 < 10
  or total_0Gy_0U_1 < 10 or total_0Gy_0U_2 < 10
  or total_1Gy_0U_1 < 10 or total_1Gy_0U_2 < 10
  or total_01Gy_100U_1 < 10 or total_01Gy_100U_2 < 10
  or total_0Gy_100U_1 < 10 or total_0Gy_100U_2 < 10
  or total_1Gy_100U_1 < 10 or total_1Gy_100U_2 < 10
  then flag_any_total_lt10=1;
  else flag_any_total_lt10=0;

  if meth_01Gy_0U_1 >= 1 or meth_01Gy_0U_2 >= 1
  or meth_0Gy_0U_1 >= 1 or meth_0Gy_0U_2 >= 1
  or meth_1Gy_0U_1 >= 1 or meth_1Gy_0U_2 >= 1
  or meth_01Gy_100U_1 >= 1 or meth_01Gy_100U_2 >= 1
  or meth_0Gy_100U_1 >= 1 or meth_0Gy_100U_2 >= 1
  or meth_1Gy_100U_1 >= 1 or meth_1Gy_100U_2 >= 1
  then flag_any_meth_ge1=1;
  else flag_any_meth_ge1=0;

  if mean(perc_01Gy_100U_1-perc_01Gy_0U_1,
          perc_01Gy_100U_2-perc_01Gy_0U_2,
          perc_0Gy_100U_1-perc_0Gy_0U_1,
          perc_0Gy_100U_2-perc_0Gy_0U_2,
          perc_1Gy_100U_1-perc_1Gy_0U_1,
          perc_1Gy_100U_2-perc_1Gy_0U_2) < 0
   then flag_0U_gt_100U=1;
   else flag_0U_gt_100U=0;
run;



data GCH_calc_diff_mean;
  retain chr pos pos2;
  set GCH_flag_sites;
  /* 0U */
  mean_perc_01Gy_0U=mean(perc_01Gy_0U_1,perc_01Gy_0U_2);
  mean_perc_0Gy_0U=mean(perc_0Gy_0U_1,perc_0Gy_0U_2);
  mean_perc_1Gy_0U=mean(perc_1Gy_0U_1,perc_1Gy_0U_2);
  mean_diff_01Gy_0Gy_0U=mean_perc_01Gy_0U-mean_perc_0Gy_0U;
  mean_diff_1Gy_0Gy_0U=mean_perc_1Gy_0U-mean_perc_0Gy_0U;

  /* 100U */
  mean_perc_01Gy_100U=mean(perc_01Gy_100U_1,perc_01Gy_100U_2);
  mean_perc_0Gy_100U=mean(perc_0Gy_100U_1,perc_0Gy_100U_2);
  mean_perc_1Gy_100U=mean(perc_1Gy_100U_1,perc_1Gy_100U_2);
  mean_diff_01Gy_0Gy_100U=mean_perc_01Gy_100U-mean_perc_0Gy_100U;
  mean_diff_1Gy_0Gy_100U=mean_perc_1Gy_100U-mean_perc_0Gy_100U;

  /* 100U - 0U */
  diff_perc_01Gy_0U_100U_1=perc_01Gy_100U_1-perc_01Gy_0U_1;
  diff_perc_01Gy_0U_100U_2=perc_01Gy_100U_2-perc_01Gy_0U_2;
  diff_perc_0Gy_0U_100U_1=perc_0Gy_100U_1-perc_0Gy_0U_1;
  diff_perc_0Gy_0U_100U_2=perc_0Gy_100U_2-perc_0Gy_0U_2;
  diff_perc_1Gy_0U_100U_1=perc_1Gy_100U_1-perc_1Gy_0U_1;
  diff_perc_1Gy_0U_100U_2=perc_1Gy_100U_2-perc_1Gy_0U_2;
  mean_diff_perc_01Gy_0U_100U=mean(diff_perc_01Gy_0U_100U_1,diff_perc_01Gy_0U_100U_2);
  mean_diff_perc_0Gy_0U_100U=mean(diff_perc_0Gy_0U_100U_1,diff_perc_0Gy_0U_100U_2);
  mean_diff_perc_1Gy_0U_100U=mean(diff_perc_1Gy_0U_100U_1,diff_perc_1Gy_0U_100U_2);
  mean_diff_perc_01Gy_0Gy=mean_diff_perc_01Gy_0U_100U-mean_diff_perc_0Gy_0U_100U;
  mean_diff_perc_1Gy_0Gy=mean_diff_perc_1Gy_0U_100U-mean_diff_perc_0Gy_0U_100U;

  pos2=pos+1;
run;


/* Export bedGraph-like text files (bedGraphs without header) */

/* By sample */

data sample1_0U;
  set GCH_calc_diff_mean;
  if perc_01Gy_0U_1 ^= . ;
  keep chr pos pos2 perc_01Gy_0U_1;
run;

data sample1_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_01Gy_0U_1 ^= .;
  keep chr pos pos2 perc_01Gy_0U_1;
run;


data sample2_0U;
  set GCH_calc_diff_mean;
  if perc_01Gy_0U_2 ^= .;
  keep chr pos pos2 perc_01Gy_0U_2;
run;

data sample2_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_01Gy_0U_2 ^= .;
  keep chr pos pos2 perc_01Gy_0U_2;
run;


data sample3_0U;
  set GCH_calc_diff_mean;
  if perc_0Gy_0U_1 ^= .;
  keep chr pos pos2 perc_0Gy_0U_1;
run;

data sample3_0U_filter;
  set GCH_calc_diff_mean;
  if perc_0Gy_0U_1 ^= .;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  keep chr pos pos2 perc_0Gy_0U_1;
run;


data sample4_0U;
  set GCH_calc_diff_mean;
  if perc_0Gy_0U_2 ^= .;
  keep chr pos pos2 perc_0Gy_0U_2;
run;

data sample4_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_0Gy_0U_2 ^= .;
  keep chr pos pos2 perc_0Gy_0U_2;
run;


data sample5_0U;
  set GCH_calc_diff_mean;
  if perc_1Gy_0U_1 ^= .;
  keep chr pos pos2 perc_1Gy_0U_1;
run;

data sample5_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_1Gy_0U_1 ^= .;
  keep chr pos pos2 perc_1Gy_0U_1;
run;


data sample6_0U;
  set GCH_calc_diff_mean;
  if perc_1Gy_0U_2 ^= .;
  keep chr pos pos2 perc_1Gy_0U_2;
run;

data sample6_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_1Gy_0U_2 ^= .;
  keep chr pos pos2 perc_1Gy_0U_2;
run;


data sample1_100U;
  set GCH_calc_diff_mean;
  if perc_01Gy_100U_1 ^= .;
  keep chr pos pos2 perc_01Gy_0U_1;
run;

data sample1_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_01Gy_100U_1 ^= .;
  keep chr pos pos2 perc_01Gy_100U_1;
run;


data sample2_100U;
  set GCH_calc_diff_mean;
  if perc_01Gy_100U_2 ^= .;
  keep chr pos pos2 perc_01Gy_100U_2;
run;

data sample2_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_01Gy_100U_2 ^= .;
  keep chr pos pos2 perc_01Gy_100U_2;
run;


data sample3_100U;
  set GCH_calc_diff_mean;
  if perc_0Gy_100U_1 ^= .;
  keep chr pos pos2 perc_0Gy_100U_1;
run;

data sample3_100U_filter;
  set GCH_calc_diff_mean;
  if perc_0Gy_100U_1 ^= .;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  keep chr pos pos2 perc_0Gy_100U_1;
run;


data sample4_100U;
  set GCH_calc_diff_mean;
  if perc_0Gy_100U_2 ^= .;
  keep chr pos pos2 perc_0Gy_100U_2;
run;

data sample4_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_0Gy_100U_2 ^= .;
  keep chr pos pos2 perc_0Gy_100U_2;
run;


data sample5_100U;
  set GCH_calc_diff_mean;
  if perc_1Gy_100U_1 ^= .;
  keep chr pos pos2 perc_1Gy_100U_1;
run;

data sample5_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_1Gy_100U_1 ^= .;
  keep chr pos pos2 perc_1Gy_100U_1;
run;


data sample6_100U;
  set GCH_calc_diff_mean;
  if perc_1Gy_100U_2 ^= .;
  keep chr pos pos2 perc_1Gy_100U_2;
run;

data sample6_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if perc_1Gy_100U_2 ^= .;
  keep chr pos pos2 perc_1Gy_100U_2;
run;


/* by sample ( unit difference) */

data sample1_100U_0U;
  set GCH_calc_diff_mean;
  if diff_perc_01Gy_0U_100U_1 ^= . ;
  keep chr pos pos2 diff_perc_01Gy_0U_100U_1;
run;

data sample1_100U_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if diff_perc_01Gy_0U_100U_1 ^= . ;
  keep chr pos pos2 diff_perc_01Gy_0U_100U_1;
run;

data sample2_100U_0U;
  set GCH_calc_diff_mean;
  if diff_perc_01Gy_0U_100U_2 ^= . ;
  keep chr pos pos2 diff_perc_01Gy_0U_100U_2;
run;

data sample2_100U_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if diff_perc_01Gy_0U_100U_2 ^= . ;
  keep chr pos pos2 diff_perc_01Gy_0U_100U_2;
run;


data sample3_100U_0U;
  set GCH_calc_diff_mean;
  if diff_perc_0Gy_0U_100U_1 ^= . ;
  keep chr pos pos2 diff_perc_0Gy_0U_100U_1;
run;

data sample3_100U_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if diff_perc_0Gy_0U_100U_1 ^= . ;
  keep chr pos pos2 diff_perc_0Gy_0U_100U_1;
run;

data sample4_100U_0U;
  set GCH_calc_diff_mean;
  if diff_perc_0Gy_0U_100U_2 ^= . ;
  keep chr pos pos2 diff_perc_0Gy_0U_100U_2;
run;

data sample4_100U_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if diff_perc_0Gy_0U_100U_2 ^= . ;
  keep chr pos pos2 diff_perc_0Gy_0U_100U_2;
run;


data sample5_100U_0U;
  set GCH_calc_diff_mean;
  if diff_perc_1Gy_0U_100U_1 ^= . ;
  keep chr pos pos2 diff_perc_1Gy_0U_100U_1;
run;

data sample5_100U_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if diff_perc_1Gy_0U_100U_1 ^= . ;
  keep chr pos pos2 diff_perc_1Gy_0U_100U_1;
run;

data sample6_100U_0U;
  set GCH_calc_diff_mean;
  if diff_perc_1Gy_0U_100U_2 ^= . ;
  keep chr pos pos2 diff_perc_1Gy_0U_100U_2;
run;

data sample6_100U_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if diff_perc_1Gy_0U_100U_2 ^= . ;
  keep chr pos pos2 diff_perc_1Gy_0U_100U_2;
run;



/* By condition */

data group1_0U_mean;
  set GCH_calc_diff_mean;
  if mean_perc_01Gy_0U ^= .;
  keep chr pos pos2 mean_perc_01Gy_0U;
run;

data group1_0U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_perc_01Gy_0U ^= .;
  keep chr pos pos2 mean_perc_01Gy_0U;
run;


data group2_0U_mean;
  set GCH_calc_diff_mean;
  if mean_perc_0Gy_0U ^= .;
  keep chr pos pos2 mean_perc_0Gy_0U;
run;

data group2_0U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_perc_0Gy_0U ^= .;
  keep chr pos pos2 mean_perc_0Gy_0U;
run;


data group3_0U_mean;
  set GCH_calc_diff_mean;
  if mean_perc_1Gy_0U ^= .;
  keep chr pos pos2 mean_perc_1Gy_0U;
run;

data group3_0U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_perc_1Gy_0U ^= .;
  keep chr pos pos2 mean_perc_1Gy_0U;
run;


data group1_100U_mean;
  set GCH_calc_diff_mean;
  if mean_perc_01Gy_100U ^= .;
  keep chr pos pos2 mean_perc_01Gy_100U;
run;

data group1_100U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_perc_01Gy_100U ^= .;
  keep chr pos pos2 mean_perc_01Gy_100U;
run;


data group2_100U_mean;
  set GCH_calc_diff_mean;
  if mean_perc_0Gy_100U ^= .;
  keep chr pos pos2 mean_perc_0Gy_100U;
run;

data group2_100U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_perc_0Gy_100U ^= .;
  keep chr pos pos2 mean_perc_0Gy_100U;
run;


data group3_100U_mean;
  set GCH_calc_diff_mean;
  if mean_perc_1Gy_100U ^= .;
  keep chr pos pos2 mean_perc_1Gy_100U;
run;

data group3_100U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_perc_1Gy_100U ^= .;
  keep chr pos pos2 mean_perc_1Gy_100U;
run;


/* By condition unit diff */


data group1_100U_0U_mean;
  set GCH_calc_diff_mean;
  if mean_diff_perc_01Gy_0U_100U ^= .;
  keep chr pos pos2 mean_diff_perc_01Gy_0U_100U;
run;

data group1_100U_0U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_perc_01Gy_0U_100U ^= .;
  keep chr pos pos2 mean_diff_perc_01Gy_0U_100U;
run;

data group2_100U_0U_mean;
  set GCH_calc_diff_mean;
  if mean_diff_perc_0Gy_0U_100U ^= .;
  keep chr pos pos2 mean_diff_perc_0Gy_0U_100U;
run;

data group2_100U_0U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_perc_0Gy_0U_100U ^= .;
  keep chr pos pos2 mean_diff_perc_0Gy_0U_100U;
run;

data group3_100U_0U_mean;
  set GCH_calc_diff_mean;
  if mean_diff_perc_1Gy_0U_100U ^= .;
  keep chr pos pos2 mean_diff_perc_1Gy_0U_100U;
run;

data group3_100U_0U_mean_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_perc_1Gy_0U_100U ^= .;
  keep chr pos pos2 mean_diff_perc_1Gy_0U_100U;
run;

/* By unit, dose diff */

data diff_01Gy_0Gy_0U;
  set GCH_calc_diff_mean;
  if mean_diff_01Gy_0Gy_0U ^= .;
  keep chr pos pos2 mean_diff_01Gy_0Gy_0U;
run;

data diff_01Gy_0Gy_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_01Gy_0Gy_0U ^= .;
  keep chr pos pos2 mean_diff_01Gy_0Gy_0U;
run;

data diff_01Gy_0Gy_100U;
  set GCH_calc_diff_mean;
  if mean_diff_01Gy_0Gy_100U ^= .;
  keep chr pos pos2 mean_diff_01Gy_0Gy_100U;
run;

data diff_01Gy_0Gy_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_01Gy_0Gy_100U ^= .;
  keep chr pos pos2 mean_diff_01Gy_0Gy_100U;
run;

data diff_1Gy_0Gy_0U;
  set GCH_calc_diff_mean;
  if mean_diff_1Gy_0Gy_0U ^= .;
  keep chr pos pos2 mean_diff_1Gy_0Gy_0U;
run;

data diff_1Gy_0Gy_0U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_1Gy_0Gy_0U ^= .;
  keep chr pos pos2 mean_diff_1Gy_0Gy_0U;
run;

data diff_1Gy_0Gy_100U;
  set GCH_calc_diff_mean;
  if mean_diff_1Gy_0Gy_100U ^= .;
  keep chr pos pos2 mean_diff_1Gy_0Gy_100U;
run;

data diff_1Gy_0Gy_100U_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_1Gy_0Gy_100U ^= .;
  keep chr pos pos2 mean_diff_1Gy_0Gy_100U;
run;


/* Difference of differences */

data diff_01Gy_0Gy;
  set GCH_calc_diff_mean;
  if mean_diff_perc_01Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_perc_01Gy_0Gy;
run;

data diff_01Gy_0Gy_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_perc_01Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_perc_01Gy_0Gy;
run;


data diff_1Gy_0Gy;
  set GCH_calc_diff_mean;
  if mean_diff_perc_1Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_perc_1Gy_0Gy;
run;

data diff_1Gy_0Gy_filter;
  set GCH_calc_diff_mean;
  where flag_any_total_lt10=0 and flag_any_meth_ge1=1 and flag_0U_gt_100U=0;
  if mean_diff_perc_1Gy_0Gy ^= .;
  keep chr pos pos2 mean_diff_perc_1Gy_0Gy;
run;


/* Export */

%macro exportDAR(datain,outname);

proc export data=&datain.
     outfile="!PATCON/arabidopsis_wgbs/analysis_output/bedGraphs/&outname..bedGraph"
     dbms=tab replace;
     putnames=no;
run;

%mend;




%exportDAR(sample1_0U,GC_01Gy_0U_rep1_all);
%exportDAR(sample1_0U_filter,GC_01Gy_0U_rep1_10reads);
%exportDAR(sample2_0U,GC_01Gy_0U_rep2_all);
%exportDAR(sample2_0U_filter,GC_01Gy_0U_rep2_10reads);
%exportDAR(sample3_0U,GC_0Gy_0U_rep1_all);
%exportDAR(sample3_0U_filter,GC_0Gy_0U_rep1_10reads);
%exportDAR(sample4_0U,GC_0Gy_0U_rep2_all);
%exportDAR(sample4_0U_filter,GC_0Gy_0U_rep2_10reads);
%exportDAR(sample5_0U,GC_1Gy_0U_rep1_all);
%exportDAR(sample5_0U_filter,GC_1Gy_0U_rep1_10reads);
%exportDAR(sample6_0U,GC_1Gy_0U_rep2_all);
%exportDAR(sample6_0U_filter,GC_1Gy_0U_rep2_10reads);


%exportDAR(sample1_100U,GC_01Gy_100U_rep1_all);
%exportDAR(sample1_100U_filter,GC_01Gy_100U_rep1_10reads);
%exportDAR(sample2_100U,GC_01Gy_100U_rep2_all);
%exportDAR(sample2_100U_filter,GC_01Gy_100U_rep2_10reads);
%exportDAR(sample3_100U,GC_0Gy_100U_rep1_all);
%exportDAR(sample3_100U_filter,GC_0Gy_100U_rep1_10reads);
%exportDAR(sample4_100U,GC_0Gy_100U_rep2_all);
%exportDAR(sample4_100U_filter,GC_0Gy_100U_rep2_10reads);
%exportDAR(sample5_100U,GC_1Gy_100U_rep1_all);
%exportDAR(sample5_100U_filter,GC_1Gy_100U_rep1_10reads);
%exportDAR(sample6_100U,GC_1Gy_100U_rep2_all);
%exportDAR(sample6_100U_filter,GC_1Gy_100U_rep2_10reads);

%exportDAR(sample1_100U_0U,GC_01Gy_100U_sub_0U_rep1_all);
%exportDAR(sample1_100U_0U_filter,GC_01Gy_100U_sub_0U_rep1_10reads);
%exportDAR(sample2_100U_0U,GC_01Gy_100U_sub_0U_rep2_all);
%exportDAR(sample2_100U_0U_filter,GC_01Gy_100U_sub_0U_rep2_10reads);
%exportDAR(sample3_100U_0U,GC_0Gy_100U_sub_0U_rep1_all);
%exportDAR(sample3_100U_0U_filter,GC_0Gy_100U_sub_0U_rep1_10reads);
%exportDAR(sample4_100U_0U,GC_0Gy_100U_sub_0U_rep2_all);
%exportDAR(sample4_100U_0U_filter,GC_0Gy_100U_sub_0U_rep2_10reads);
%exportDAR(sample5_100U_0U,GC_1Gy_100U_sub_0U_rep1_all);
%exportDAR(sample5_100U_0U_filter,GC_1Gy_100U_sub_0U_rep1_10reads);
%exportDAR(sample6_100U_0U,GC_1Gy_100U_sub_0U_rep2_all);
%exportDAR(sample6_100U_0U_filter,GC_1Gy_100U_sub_0U_rep2_10reads);

%exportDAR(group1_0U_mean,GC_01Gy_0U_avg_all);
%exportDAR(group1_0U_mean_filter,GC_01Gy_0U_avg_10reads);
%exportDAR(group2_0U_mean,GC_0Gy_0U_avg_all);
%exportDAR(group2_0U_mean_filter,GC_0Gy_0U_avg_10reads);
%exportDAR(group3_0U_mean,GC_1Gy_0U_avg_all);
%exportDAR(group3_0U_mean_filter,GC_1Gy_0U_avg_10reads);

%exportDAR(group1_100U_mean,GC_01Gy_100U_avg_all);
%exportDAR(group1_100U_mean_filter,GC_01Gy_100U_avg_10reads);
%exportDAR(group2_100U_mean,GC_0Gy_100U_avg_all);
%exportDAR(group2_100U_mean_filter,GC_0Gy_100U_avg_10reads);
%exportDAR(group3_100U_mean,GC_1Gy_100U_avg_all);
%exportDAR(group3_100U_mean_filter,GC_1Gy_100U_avg_10reads);

%exportDAR(group1_100U_0U_mean,GC_01Gy_100U_sub_0U_avg_all);
%exportDAR(group1_100U_0U_mean_filter,GC_01Gy_100U_sub_0U_avg_10reads);
%exportDAR(group2_100U_0U_mean,GC_0Gy_100U_sub_0U_avg_all);
%exportDAR(group2_100U_0U_mean_filter,GC_0Gy_100U_sub_0U_avg_10reads);
%exportDAR(group3_100U_0U_mean,GC_1Gy_100U_sub_0U_avg_all);
%exportDAR(group3_100U_0U_mean_filter,GC_1Gy_100U_sub_0U_avg_10reads);

%exportDAR(diff_01Gy_0Gy_0U,GC_01Gy_sub_0Gy_0U_avg_all);
%exportDAR(diff_01Gy_0Gy_0U_filter,GC_01Gy_sub_0Gy_0U_avg_10reads);
%exportDAR(diff_01Gy_0Gy_100U,GC_01Gy_sub_0Gy_100U_avg_all);
%exportDAR(diff_01Gy_0Gy_100U_filter,GC_01Gy_sub_0Gy_100U_avg_10reads);

%exportDAR(diff_1Gy_0Gy_0U,GC_1Gy_sub_0Gy_0U_avg_all);
%exportDAR(diff_1Gy_0Gy_0U_filter,GC_1Gy_sub_0Gy_0U_avg_10reads);
%exportDAR(diff_1Gy_0Gy_100U,GC_1Gy_sub_0Gy_100U_avg_all);
%exportDAR(diff_1Gy_0Gy_100U_filter,GC_1Gy_sub_0Gy_100U_avg_10reads);

%exportDAR(diff_01Gy_0Gy,GC_01Gy_100U_0U_sub_0Gy_100U_0U_avg_all);
%exportDAR(diff_01Gy_0Gy_filter,GC_01Gy_100U_0U_sub_0Gy_100U_0U_avg_10reads);
%exportDAR(diff_1Gy_0Gy,GC_1Gy_100U_0U_sub_0Gy_100U_0U_avg_all);
%exportDAR(diff_1Gy_0Gy_filter,GC_1Gy_100U_0U_sub_0Gy_100U_0U_avg_10reads);




