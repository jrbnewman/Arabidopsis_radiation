/* Counts for arabidopsis paper */

libname brassMAP '/TB14/TB14/sandbox/dtra_sandbox/brass_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

/* (1) Import BLASTN results for Br to At transcripts
   (2) Import Br to At orthologs list
   (3) Remove from BLASTN results any pairs not in ortholog list
   (4) Calculate distribution of mismatches and gaps 

   (5) Import At TE to Br BLASTN results
   (6) Use parameters derived from (4) to filter hits
   (7) Filter for >90% matches
   (8) Make a Brassica rapa TE database

*/

proc import datafile="/TB14/TB14/sandbox/brassica_sandbox/blast_output/Br_transcripts_At_blastn.tsv"
   out=xs_br2at dbms=tab replace;
   guessingrows=max;
   getnames=no;
run;


proc import datafile="!HOME/concannon/useful_brassica_data/Brapa_R500_v1_6_PEdgar/Brapa_R500_v1.6/R500_Ath_Orthologues.csv"
   out=ortho dbms=csv replace;
   guessingrows=max;
   getnames=yes;
run;

data xs_br2at2;
  set xs_br2at;
  length br_geneID $20.;
  length at_geneID $20.;
  br_geneID = compress(scan(VAR1,1,":"));
  at_geneID = compress(scan(VAR2,1,"."));
  perc_match = (1- VAR5 / VAR4) * 100 ; 
  rename VAR1=br_query
         VAR2=at_ref
         VAR3=perc_identity
         VAR4=length
         VAR5=mismatch
         VAR6=gap_open
         VAR7=q_start
         VAR8=q_end
         VAR9=s_start
         VAR10=s_end
         VAR11=evalue
         VAR12=bitscore;
run;

data ortho_2;
  set ortho;
  if Ath_ortho = "NS" then delete;
  if Ath_ortho = "NA" then delete;
  if Ath_ortho = " " then delete;
  if Ath_ortho = "" then delete;
  keep BRA Ath_ortho;
  rename BRA=br_geneID ath_ortho=at_geneID;
run;

proc sort data=ortho_2 nodup;
  by br_geneID at_geneID;
proc sort data=xs_br2at2;
  by br_geneID at_geneID;
run;

data xs_br2at3;
  merge xs_br2at2 (in=in1) ortho_2 (in=in2);
  by br_geneID at_geneID;
  if in1 and in2;
run;

proc means data=xs_br2at3 noprint;
  var mismatch;
  output out=mm_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;


proc means data=xs_br2at3 noprint;
  var perc_match;
  output out=pmm_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;

proc means data=xs_br2at3 noprint;
  var gap_open;
  output out=gap_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;

proc means data=xs_br2at3 noprint;
  var perc_identity;
  output out=ID_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;

proc print data=mm_stats;
proc print data=pmm_stats;
proc print data=gap_stats;
proc print data=ID_stats;
run;

/*

Mismatch

_FREQ_            mean              sd             min              q1           media              q3             max

 61241    118.55306086    107.71976896               0              35              92             171            1427



Perc match ((1 - MM/length) * 100)

 Obs    _TYPE_    _FREQ_      mean        sd      min       q1       media        q3      max

  1        0       61241    88.1291    2.85053     80    86.2963    88.2143    89.9865    100


Gaps


_FREQ_            mean              sd             min              q1           media              q3             max

 61241    5.9365457782     7.058484128               0               1               3               9              80


Perc ID

  _FREQ_            mean              sd             min              q1           media              q3             max

   61241    85.626166131    3.4700098176              80          82.857          85.335          87.967             100


Parameters:
ignore number of mismatches (for now)
Minimum % Identity is 90%
Max gaps is 3 (median of transcript matches)


*/




data xs_br2at3;
  merge xs_br2at2 (in=in1) ortho_2 (in=in2);
  by br_geneID at_geneID;
  if in1 and in2;
run;

proc sort data=xs_br2at3;
  by br_geneID at_geneID descending perc_identity;
run;

data xs_br2at4;
  set xs_br2at3;
  by br_geneID at_geneID;
  if first.at_geneID then output;
run;


proc means data=xs_br2at4 noprint;
  var mismatch;
  output out=mm_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;


proc means data=xs_br2at4 noprint;
  var perc_match;
  output out=pmm_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;

proc means data=xs_br2at4 noprint;
  var gap_open;
  output out=gap_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;

proc means data=xs_br2at4 noprint;
  var perc_identity;
  output out=ID_stats mean=mean stddev=sd min=min q1=q1 median=media q3=q3 max=max;
run;

proc print data=mm_stats;
proc print data=pmm_stats;
proc print data=gap_stats;
proc print data=ID_stats;
run;



/*

Mismatch

  _FREQ_            mean              sd             min              q1           media              q3             max

   25235    121.96437488    103.63527197               0              40              99             176            1427

Perc match ((1 - MM/length) * 100)


   _FREQ_      mean        sd      min       q1       media     q3    max

    25235    88.1377    2.90183     80    86.2396    88.1720    90    100


Gaps

 _FREQ_            mean              sd             min              q1           media              q3             max

  25235    5.7444422429    6.8021173431               0               1               3               9              57

Perc ID


 _FREQ_            mean              sd             min              q1           media              q3             max

  25235    85.885979077    3.5529629016              80            83.1          85.578          88.235             100

Parameters:
ignore number of mismatches (for now)
Minimum % Identity is 90%
Max gaps is 3 (median of transcript matches)

*/


proc import datafile="/TB14/TB14/sandbox/brassica_sandbox/blast_output/At_TE_fragments_Br_blastn.tsv"
   out=te_at2br dbms=tab replace;
   guessingrows=max;
   getnames=no;
run;




data te_at2br2;
  set te_at2br;
  length at_TE_ID $20.;
  q_length = abs(VAR8 - VAR7) + 1;
  s_length = abs(VAR10 - VAR9) + 1;
  TE_length = scan(scan(VAR1,3,":"),2,"-") - scan(scan(VAR1,3,":"),1,"-"); 
  len_diff = TE_length - VAR4;
  rename VAR1=at_query
         VAR2=br_ref
         VAR3=perc_identity
         VAR4=length
         VAR5=mismatch
         VAR6=gap_open
         VAR7=q_start
         VAR8=q_end
         VAR9=s_start
         VAR10=s_end
         VAR11=evalue
         VAR12=bitscore;
  at_TE_ID = compress(scan(VAR1,1,":"));
run;


data te_at2br2a;
  set te_at2br2;
  if abs(len_diff) / TE_length > 0.25 then delete;
run;


/* Check gap sizes */

data te_at2br3;
   set te_at2br2;
   diff_len_qlen = abs(length - q_length);
   diff_len_slen = abs(length - s_length);
run;




data te_at2br2_filtered_85 te_at2br2_filtered_90 te_at2br2_filtered_95 te_at2br2_filtered_100;
   set te_at2br2;
   where  gap_open <= 3   ;
   if perc_identity >= 85 then output te_at2br2_filtered_85;
   if perc_identity >= 90 then output te_at2br2_filtered_90;
   if perc_identity >= 95 then output te_at2br2_filtered_95;
   if perc_identity >= 100 then output te_at2br2_filtered_100;
run;


data te_at2br2_filtered_85_A te_at2br2_filtered_90_A te_at2br2_filtered_95_A te_at2br2_filtered_100_A;
   set te_at2br2;
   where br_ref ? "A" and gap_open <= 3;
   if perc_identity >= 85 then output te_at2br2_filtered_85_A;
   if perc_identity >= 90 then output te_at2br2_filtered_90_A;
   if perc_identity >= 95 then output te_at2br2_filtered_95_A;
   if perc_identity >= 100 then output te_at2br2_filtered_100_A;
run;



data te_at2br2_filtered_85_Agap te_at2br2_filtered_90_Agap te_at2br2_filtered_95_Agap te_at2br2_filtered_100_Agap; 
   set te_at2br2;
   where br_ref ? "A" ;
   if perc_identity >= 85 then output te_at2br2_filtered_85_Agap;
   if perc_identity >= 90 then output te_at2br2_filtered_90_Agap;
   if perc_identity >= 95 then output te_at2br2_filtered_95_Agap;
   if perc_identity >= 100 then output te_at2br2_filtered_100_Agap;
run;



data te_at2br2_filtered_85_gap te_at2br2_filtered_90_gap te_at2br2_filtered_95_gap te_at2br2_filtered_100_gap;
   set te_at2br2;

   if perc_identity >= 85 then output te_at2br2_filtered_85_gap;
   if perc_identity >= 90 then output te_at2br2_filtered_90_gap;
   if perc_identity >= 95 then output te_at2br2_filtered_95_gap;
   if perc_identity >= 100 then output te_at2br2_filtered_100_gap;
run;

/*

There were 612757 observations read from the data set WORK.TE_AT2BR2.
WHERE gap_open<=3;
The data set WORK.TE_AT2BR2_FILTERED_85 has 589678 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_90 has 550421 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_95 has 272640 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_100 has 41820 observations and 15 variables.

There were 64583 observations read from the data set WORK.TE_AT2BR2.
 WHERE br_ref contains 'A' and (gap_open<=3);
 The data set WORK.TE_AT2BR2_FILTERED_85_A has 48420 observations and 15 variables.
 The data set WORK.TE_AT2BR2_FILTERED_90_A has 37873 observations and 15 variables.
 The data set WORK.TE_AT2BR2_FILTERED_95_A has 16676 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_100_A has 2896 observations and 15 variables.

There were 69611 observations read from the data set WORK.TE_AT2BR2.
WHERE br_ref contains 'A';
The data set WORK.TE_AT2BR2_FILTERED_85_AGAP has 49147 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_90_AGAP has 37878 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_95_AGAP has 16676 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_100_A has 2896 observations and 15 variables.


 There were 619674 observations read from the data set WORK.TE_AT2BR2.
 The data set WORK.TE_AT2BR2_FILTERED_85_GAP has 590675 observations and 15 variables.
 The data set WORK.TE_AT2BR2_FILTERED_90_GAP has 550427 observations and 15 variables.
 The data set WORK.TE_AT2BR2_FILTERED_95_GAP has 272640 observations and 15 variables.
The data set WORK.TE_AT2BR2_FILTERED_100 has 41820 observations and 15 variables.
*/

proc freq data=te_at2br2 noprint;
   where br_ref ? "A" and gap_open <= 3;
 tables at_query / out=cnt_hits_per_TE_80;
run;

proc means data=cnt_hits_per_TE_80 noprint;
 var count;
 output out=distrib_hits_80 mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;

proc freq data=te_at2br2_filtered_85_A noprint;
 tables at_query / out=cnt_hits_per_TE_85;
run;

proc means data=cnt_hits_per_TE_85 noprint;
 var count;
 output out=distrib_hits_85 mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;


proc freq data=te_at2br2_filtered_90_A noprint;
 tables at_query / out=cnt_hits_per_TE_90;
run;

proc means data=cnt_hits_per_TE_90 noprint;
 var count;
 output out=distrib_hits_90 mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;


proc freq data=te_at2br2_filtered_95_A noprint;
 tables at_query / out=cnt_hits_per_TE_95;
run;

proc means data=cnt_hits_per_TE_95 noprint;
 var count;
 output out=distrib_hits_95 mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;


proc freq data=te_at2br2_filtered_100_A noprint;
 tables at_query / out=cnt_hits_per_TE_100;
run;

proc means data=cnt_hits_per_TE_100 noprint;
 var count;
 output out=distrib_hits_100 mean=mean stddev=sd min=min q1=q1 median=median q3=q3 max=max;
run;

proc print data=distrib_hits_80;
proc print data=distrib_hits_85;
proc print data=distrib_hits_90;
proc print data=distrib_hits_95;
proc print data=distrib_hits_100;
run;

/*

%ID  _FREQ_      mean        sd      min    q1    median    q3    max

80  4771     13.5366    26.3501     1      1       3      10    195
85  2884     16.7892    30.9351     1      1       3      16    195
90  1434     26.4107    37.6170     1      1       4      34    195
95   524     31.8244    38.1783     1      1      7.5     67    156
100  90      32.1778    33.0347     1      2      28      67    100

Question: how many hits do arabidopsis TEs have to their own genome at 100 percent identity?

*/




