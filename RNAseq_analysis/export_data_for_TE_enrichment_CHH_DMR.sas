/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/output/TE_fragments.bed"
   out=te_bed dbms=tab replace;
   guessingrows=max;
   getnames=no;
run;

data te2fragment;
  set te_bed;
  length geneID $20.;
  length teID $20.;
  geneID = compress(scan(VAR4,2,"/"));
  teID = compress(scan(VAR4,1,"/"));
  keep geneID teID;
run;



data dmr_01_chh_te dmr_1_chh_te;
     set arabMAP.results_by_dmr_annot_te;
     where site_type="CHH";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output dmr_01_chh_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output dmr_01_chh_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output dmr_1_chh_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output dmr_1_chh_te;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output dmr_01_chh_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output dmr_1_chh_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output dmr_01_chh_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output dmr_1_chh_te;

     keep comparison site_type chr region_Start region_stop geneID;
run;



data dmr_01_chh_te_all dmr_1_chh_te_all;
     set arabMAP.results_by_dmr_annot_te;
     length feature $32.;
     where site_type="CHH";
     feature = scan(annotation, 1, " ");
     if feature = "Intergenic" then geneID="";
     if feature = "" then geneID="";
     /* FET */
     if comparison="0Gy_vs_01G" then output dmr_01_chh_te_all;
     if comparison="0Gy_vs_1Gy" then output dmr_1_chh_te_all;
     keep comparison site_type chr region_Start region_stop geneID;
run;

proc sort data=dmr_01_chh_te_all nodup;
  by comparison site_type chr region_start region_stop;
proc sort data=dmr_1_chh_te_all nodup;
  by comparison site_type chr region_start region_stop;
proc sort data=dmr_01_chh_te nodup;
  by comparison site_type chr region_start region_stop;
proc sort data=dmr_1_chh_te nodup;
  by comparison site_type chr region_start region_stop;
run;

data dmr_01_chh_flag;
  merge dmr_01_chh_te_all (in=in1) dmr_01_chh_te (in=in2);
  by comparison site_type chr region_start region_stop;
  if in2 then flag_dmr=1; else flag_dmr=0;
  if in1 then output;
run;

data dmr_1_chh_flag;
  merge dmr_1_chh_te_all (in=in1) dmr_1_chh_te (in=in2);
  by comparison site_type chr region_start region_stop;
  if in2 then flag_dmr=1; else flag_dmr=0;
  if in1 then output;
run;

proc sort data=dmr_01_chh_flag;
  by geneID;
proc sort data=dmr_1_chh_flag;
  by geneID;
proc sort data=te2fragment;
  by geneID;
run;

data dmr_01_chh_flag2;
  merge dmr_01_chh_flag (in=in1) te2fragment (in=in2);
  by geneID;
  if in1;
run;

data dmr_1_chh_flag2;
  merge dmr_1_chh_flag (in=in1) te2fragment (in=in2);
  by geneID;
  if in1;
run;


proc export data=dmr_01_chh_flag2 
outfile="!HOME/concannon/DTRA/DMR_10cgy_CHH_TE_annot.csv"
dbms=csv replace;
run;

proc export data=dmr_1_chh_flag2 
outfile="!HOME/concannon/DTRA/DMR_100cgy_CHH_TE_annot.csv"
dbms=csv replace;
run;


