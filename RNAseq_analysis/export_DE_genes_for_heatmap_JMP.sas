ods listing; ods html close;
libname arabrna '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname tair '!HOME/concannon/useful_arabidopsis_data/TAIR10/sas_data';

proc datasets lib=work kill noprint;
run;
quit;



/* get data from summary table */

data de_genes_fc;
   set arabrna.arab_results_by_gene;
   if (flag_01gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >=1 )
   or (flag_01gy_v_Mock_3h_fdr05=1 and abs(log2_fc_01gy_v_Mock_3h) >=1 )
   or (flag_01gy_v_Mock_24h_fdr05=1 and abs(log2_fc_01gy_v_Mock_24h) >=1 )
   or (flag_01gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >=1 )
   or (flag_1gy_v_Mock_1h_fdr05=1 and abs(log2_fc_1gy_v_Mock_lh) >=1 )
   or (flag_1gy_v_Mock_3h_fdr05=1 and abs(log2_fc_1gy_v_Mock_3h) >=1 )
   or (flag_1gy_v_Mock_24h_fdr05=1 and abs(log2_fc_1gy_v_Mock_24h) >=1 )
   or (flag_1gy_v_Mock_72h_fdr05=1 and abs(log2_fc_1gy_v_Mock_72h) >=1 );
   keep gene_id symbol  log2_fc_: ;
run;



data de_genes_fc2;
   set arabrna.arab_results_by_gene;

       if log2_fc_01gy_v_Mock_lh=.
       or log2_fc_01gy_v_Mock_3h=.
     or log2_fc_01gy_v_Mock_24h=.
      or log2_fc_01gy_v_Mock_72h=.
      or log2_fc_1gy_v_Mock_lh=.
      or log2_fc_1gy_v_Mock_3h=.
      or log2_fc_1gy_v_Mock_24h=.
      or log2_fc_1gy_v_Mock_72h=. then delete;



   if (flag_01gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >=1 )
   or (flag_01gy_v_Mock_3h_fdr05=1 and abs(log2_fc_01gy_v_Mock_3h) >=1 )
   or (flag_01gy_v_Mock_24h_fdr05=1 and abs(log2_fc_01gy_v_Mock_24h) >=1 )
   or (flag_01gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >=1 )
   or (flag_1gy_v_Mock_1h_fdr05=1 and abs(log2_fc_1gy_v_Mock_lh) >=1 )
   or (flag_1gy_v_Mock_3h_fdr05=1 and abs(log2_fc_1gy_v_Mock_3h) >=1 )
   or (flag_1gy_v_Mock_24h_fdr05=1 and abs(log2_fc_1gy_v_Mock_24h) >=1 )
   or (flag_1gy_v_Mock_72h_fdr05=1 and abs(log2_fc_1gy_v_Mock_72h) >=1 );


   keep gene_id symbol  log2_fc_: ;
run;


data de_genes_mean;
   set arabrna.arab_results_by_gene;
   if (flag_01gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >=1 )
   or (flag_01gy_v_Mock_3h_fdr05=1 and abs(log2_fc_01gy_v_Mock_3h) >=1 )
   or (flag_01gy_v_Mock_24h_fdr05=1 and abs(log2_fc_01gy_v_Mock_24h) >=1 )
   or (flag_01gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >=1 )
   or (flag_1gy_v_Mock_1h_fdr05=1 and abs(log2_fc_1gy_v_Mock_lh) >=1 )
   or (flag_1gy_v_Mock_3h_fdr05=1 and abs(log2_fc_1gy_v_Mock_3h) >=1 )
   or (flag_1gy_v_Mock_24h_fdr05=1 and abs(log2_fc_1gy_v_Mock_24h) >=1 )
   or (flag_1gy_v_Mock_72h_fdr05=1 and abs(log2_fc_1gy_v_Mock_72h) >=1 );
   keep gene_id symbol mean_cpm_: ;
run;


proc export data=de_genes_fc outfile="!HOME/concannon/DTRA/arabidopsis_FC1_any_DE_genes_JMP.csv"
   dbms=csv replace;
run;

proc export data=de_genes_fc2 outfile="!HOME/concannon/DTRA/arabidopsis_FC1_any_DE_genes_nomiss_JMP.csv"
   dbms=csv replace;
run;

proc export data=de_genes_mean outfile="!HOME/concannon/DTRA/arabidopsis_means_FC1_any_DE_genes_JMP.csv"
   dbms=csv replace;
run;


data de_genes_counts;
   set arabrna.arab_results_by_gene;
   if log2_fc_01gy_v_Mock_lh = . then sign_01_1=.;
   else if log2_fc_01gy_v_Mock_lh >= 1 then sign_01_1=1;
   else if log2_fc_01gy_v_Mock_lh <= -1 then sign_01_1=-1;
   else sign_01_1=0;

   if log2_fc_01gy_v_Mock_3h = . then sign_01_3=.;
   else if log2_fc_01gy_v_Mock_3h >= 1 then sign_01_3=1;
   else if log2_fc_01gy_v_Mock_3h <= -1 then sign_01_3=-1;
   else sign_01_3=0;

   if log2_fc_01gy_v_Mock_24h = . then sign_01_24=.;
   else if log2_fc_01gy_v_Mock_24h >= 1 then sign_01_24=1;
   else if log2_fc_01gy_v_Mock_24h <= -1 then sign_01_24=-1;
   else sign_01_24=0;

   if log2_fc_01gy_v_Mock_72h = . then sign_01_72=.;
   else if log2_fc_01gy_v_Mock_72h >= 1 then sign_01_72=1;
   else if log2_fc_01gy_v_Mock_72h <= -1 then sign_01_72=-1;
   else sign_01_72=0;

   if log2_fc_1gy_v_Mock_lh = . then sign_1_1=.;
   else if log2_fc_1gy_v_Mock_lh >= 1 then sign_1_1=1;
   else if log2_fc_1gy_v_Mock_lh <= -1 then sign_1_1=-1;
   else sign_1_1=0;

   if log2_fc_1gy_v_Mock_3h = . then sign_1_3=.;
   else if log2_fc_1gy_v_Mock_3h >= 1 then sign_1_3=1;
   else if log2_fc_1gy_v_Mock_3h <= -1 then sign_1_3=-1;
   else sign_1_3=0;

   if log2_fc_1gy_v_Mock_24h = . then sign_1_24=.;
   else if log2_fc_1gy_v_Mock_24h >= 1 then sign_1_24=1;
   else if log2_fc_1gy_v_Mock_24h <= -1 then sign_1_24=-1;
   else sign_1_24=0;

   if log2_fc_1gy_v_Mock_72h = . then sign_1_72=.;
   else if log2_fc_1gy_v_Mock_72h >= 1 then sign_1_72=1;
   else if log2_fc_1gy_v_Mock_72h <= -1 then sign_1_72=-1;
   else sign_1_72=0;

run;

proc freq data=de_genes_counts;
  tables flag_01gy_v_Mock_1h_fdr05*sign_01_1 
         flag_01gy_v_Mock_3h_fdr05*sign_01_3 
         flag_01gy_v_Mock_24h_fdr05*sign_01_24 
         flag_01gy_v_Mock_72h_fdr05*sign_01_72 
         flag_1gy_v_Mock_1h_fdr05*sign_01_1 
         flag_1gy_v_Mock_3h_fdr05*sign_01_3 
         flag_1gy_v_Mock_24h_fdr05*sign_01_24 
         flag_1gy_v_Mock_72h_fdr05*sign_01_72 ;
run;


/*
                FC1         ANY
Dose    Time    Down    Up  Down    Up  
10      1       109     387 2058    2508
10      3       11      22  58      52
10      24      6       23  110     394
10      72      14      50  787     758
100     1       66      263 1269    1905
100     3       14      65  719     1375
100     24      0       0   31      16
100     72      4       12  28      107

*/




data de_genes_counts;
   set arabrna.arab_results_by_gene;
   if log2_fc_01gy_v_Mock_lh = . then sign_01_1=.;
   else if log2_fc_01gy_v_Mock_lh > 0 then sign_01_1=1;
   else if log2_fc_01gy_v_Mock_lh < 0 then sign_01_1=-1;
   else sign_01_1=0;

   if log2_fc_01gy_v_Mock_3h = . then sign_01_3=.;
   else if log2_fc_01gy_v_Mock_3h > 0 then sign_01_3=1;
   else if log2_fc_01gy_v_Mock_3h < 0 then sign_01_3=-1;
   else sign_01_3=0;

   if log2_fc_01gy_v_Mock_24h = . then sign_01_24=.;
   else if log2_fc_01gy_v_Mock_24h > 0 then sign_01_24=1;
   else if log2_fc_01gy_v_Mock_24h < 0 then sign_01_24=-1;
   else sign_01_24=0;

   if log2_fc_01gy_v_Mock_72h = . then sign_01_72=.;
   else if log2_fc_01gy_v_Mock_72h > 0 then sign_01_72=1;
   else if log2_fc_01gy_v_Mock_72h < 0 then sign_01_72=-1;
   else sign_01_72=0;

   if log2_fc_1gy_v_Mock_lh = . then sign_1_1=.;
   else if log2_fc_1gy_v_Mock_lh > 0 then sign_1_1=1;
   else if log2_fc_1gy_v_Mock_lh < 0 then sign_1_1=-1;
   else sign_1_1=0;

   if log2_fc_1gy_v_Mock_3h = . then sign_1_3=.;
   else if log2_fc_1gy_v_Mock_3h > 0 then sign_1_3=1;
   else if log2_fc_1gy_v_Mock_3h < 0 then sign_1_3=-1;
   else sign_1_3=0;

   if log2_fc_1gy_v_Mock_24h = . then sign_1_24=.;
   else if log2_fc_1gy_v_Mock_24h > 0 then sign_1_24=1;
   else if log2_fc_1gy_v_Mock_24h < 0 then sign_1_24=-1;
   else sign_1_24=0;

   if log2_fc_1gy_v_Mock_72h = . then sign_1_72=.;
   else if log2_fc_1gy_v_Mock_72h > 0 then sign_1_72=1;
   else if log2_fc_1gy_v_Mock_72h <=0 then sign_1_72=-1;
   else sign_1_72=0;

run;

proc freq data=de_genes_counts;
  tables flag_01gy_v_Mock_1h_fdr05*sign_01_1 
         flag_01gy_v_Mock_3h_fdr05*sign_01_3 
         flag_01gy_v_Mock_24h_fdr05*sign_01_24 
         flag_01gy_v_Mock_72h_fdr05*sign_01_72 
         flag_1gy_v_Mock_1h_fdr05*sign_01_1 
         flag_1gy_v_Mock_3h_fdr05*sign_01_3 
         flag_1gy_v_Mock_24h_fdr05*sign_01_24 
         flag_1gy_v_Mock_72h_fdr05*sign_01_72 ;
run;


data early_01 early_1 late_01 late_1;
   set arabrna.arab_results_by_gene;
   if (flag_01gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >=1) then output early_01;
   if (flag_1gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >=1) then output early_1;
   if (flag_01gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >=1) then output late_01;
   if (flag_1gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >=1) then output late_1;
  keep gene_id;
run;

proc sort data=early_01;
 by gene_id;
proc sort data=early_1;
 by gene_id;
proc sort data=late_01;
 by gene_id;
proc sort data=late_1;
 by gene_id;
run;

data compare_early_vs_late;
  merge early_01 (in=in1) early_1 (in=in2) late_01 (in=in3) late_1 (in=in4);
  by gene_id;
  if in1 then DE_01gy_1h=1; else DE_01gy_1h=0;
  if in2 then DE_1gy_1h=1; else DE_1gy_1h=0;
  if in3 then DE_01gy_72h=1; else DE_01gy_72h=0;
  if in4 then DE_1gy_72h=1; else DE_1gy_72h=0;
run;

proc freq data=compare_early_vs_late noprint;
  tables DE_01gy_1h*DE_1gy_1h*DE_01gy_72h*DE_1gy_72h / out=compare;
proc print data=compare;
run;

/*
DE_01gy_    DE_1gy_    DE_01gy_    DE_1gy_
   1h          1h         72h        72h      COUNT

    0          0           0          1          7
    0          0           1          0         51
    0          0           1          1          8
    0          1           0          0         13
    1          0           0          0        176
    1          0           1          0          3
    1          0           1          1          1
    1          1           0          0        315
    1          1           1          0          1

*/


data early_01 early_1 late_01 late_1;
   set arabrna.arab_results_by_gene;
   if (flag_01gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >0) then output early_01;
   if (flag_1gy_v_Mock_1h_fdr05=1 and abs(log2_fc_01gy_v_Mock_lh) >0) then output early_1;
   if (flag_01gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >0) then output late_01;
   if (flag_1gy_v_Mock_72h_fdr05=1 and abs(log2_fc_01gy_v_Mock_72h) >0) then output late_1;
  keep gene_id;
run;

proc sort data=early_01;
 by gene_id;
proc sort data=early_1;
 by gene_id;
proc sort data=late_01;
 by gene_id;
proc sort data=late_1;
 by gene_id;
run;

data compare_early_vs_late;
  merge early_01 (in=in1) early_1 (in=in2) late_01 (in=in3) late_1 (in=in4);
  by gene_id;
  if in1 then DE_01gy_1h=1; else DE_01gy_1h=0;
  if in2 then DE_1gy_1h=1; else DE_1gy_1h=0;
  if in3 then DE_01gy_72h=1; else DE_01gy_72h=0;
  if in4 then DE_1gy_72h=1; else DE_1gy_72h=0;
run;

proc freq data=compare_early_vs_late noprint;
  tables DE_01gy_1h*DE_1gy_1h*DE_01gy_72h*DE_1gy_72h / out=compare;
proc print data=compare;
run;

/*
DE_01gy_    DE_1gy_    DE_01gy_    DE_1gy_
   1h          1h         72h        72h      COUNT

    0          0           0          1          49
    0          0           1          0        1069
    0          0           1          1          35
    0          1           0          0         519
    0          1           1          0          32
    0          1           1          1           3
    1          0           0          0        1738
    1          0           0          1           7
    1          0           1          0         189
    1          0           1          1          12
    1          1           0          0        2405
    1          1           0          1          10
    1          1           1          0         186
    1          1           1          1          19


*/







