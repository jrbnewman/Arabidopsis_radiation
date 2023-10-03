/* Check concordance of results between gene-level analyses */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

/* To check:
(1) Genes on at each condition
(2) Genes DE (P<0.05, FDR P<0.05) each comparison
(3) Genes significantly up/down each comparison
(4) Final U/P/off/etc. calls each comparison
*/


/* (1) check genes on per condition. These should be close to identical */

data old_dtct_flags;
  set rs.arab_flag_gene_on_gt0;
run;

data new_dtct_flags;
  set rs.arab_flag_gene_on_cpm_gt0;
run;

proc sort data=old_dtct_flags;
  by gene_id;
proc sort data=new_dtct_flags;
  by gene_id;
run;

data dtct_flag_compare oops1 oops2;
  merge old_dtct_flags (in=in1) new_dtct_flags (in=in2);
  by gene_id;
  if in1 and in2 then output dtct_flag_compare;
  else if in1 then output oops1;
  else output oops2; *these are genes that were "off" at the fusion level, so I will look at just the intersection;
run;

data dtct_flag_compare2;
  set dtct_flag_compare;
  array change _numeric_;
  do over change;
    if change=. then change=0;
    end;
run;

proc freq data=dtct_flag_compare2;
   tables flag_gene_on_mock_apn0*flag_gene_on_mock_cpm0 flag_gene_on_01gy_apn0*flag_gene_on_01gy_cpm0
          flag_gene_on_1gy_apn0*flag_gene_on_1gy_cpm0 flag_gene_on_1hr_apn0*flag_gene_on_1hr_cpm0
          flag_gene_on_3hr_apn0*flag_gene_on_3hr_cpm0 flag_gene_on_34hr_apn0*flag_gene_on_24hr_cpm0
          flag_gene_on_72hr_apn0*flag_gene_on_72hr_cpm0 flag_gene_on_mock_1hr_apn0*flag_gene_on_mock_1hr_cpm0
          flag_gene_on_mock_3hr_apn0*flag_gene_on_mock_3hr_cpm0 flag_gene_on_mock_24hr_apn0*flag_gene_on_mock_24hr_cpm0
          flag_gene_on_mock_72hr_apn0*flag_gene_on_mock_72hr_cpm0 flag_gene_on_01gy_1hr_apn0*flag_gene_on_01gy_1hr_cpm0
          flag_gene_on_01gy_3hr_apn0*flag_gene_on_01gy_3hr_cpm0 flag_gene_on_01gy_24hr_apn0*flag_gene_on_01gy_24hr_cpm0
          flag_gene_on_01gy_72hr_apn0*flag_gene_on_01gy_72hr_cpm0 flag_gene_on_1gy_1hr_apn0*flag_gene_on_1gy_1hr_cpm0
          flag_gene_on_1gy_3hr_apn0*flag_gene_on_1gy_3hr_cpm0 flag_gene_on_1gy_24hr_apn0*flag_gene_on_1gy_24hr_cpm0
          flag_gene_on_1gy_72hr_apn0*flag_gene_on_1gy_72hr_cpm0;
run;

data oops2_2;
  set oops2;
    array change _numeric_;
  do over change;
    if change=. then change=0;
    end;
run;

proc freq data=oops2_2;
   tables flag_gene_on_mock_cpm0 flag_gene_on_01gy_cpm0
          flag_gene_on_1gy_cpm0 flag_gene_on_1hr_cpm0
          flag_gene_on_3hr_cpm0 flag_gene_on_24hr_cpm0
          flag_gene_on_72hr_cpm0 flag_gene_on_mock_1hr_cpm0
          flag_gene_on_mock_3hr_cpm0 flag_gene_on_mock_24hr_cpm0
          flag_gene_on_mock_72hr_cpm0 flag_gene_on_01gy_1hr_cpm0
          flag_gene_on_01gy_3hr_cpm0 flag_gene_on_01gy_24hr_cpm0
          flag_gene_on_01gy_72hr_cpm0 flag_gene_on_1gy_1hr_cpm0
          flag_gene_on_1gy_3hr_cpm0 flag_gene_on_1gy_24hr_cpm0
          flag_gene_on_1gy_72hr_cpm0;
run;

/*
CONDITION	OFF		APN0	CP160M0	BOTH	CPM0-only (off for fusions)
Mock		7653	892		164		20355	1868
0.1gy		7613	928		161		20362	1855
1gy			7490	941		159		20474	1873
1hr			7798	899		160		20207	1818
3hr			7712	896		160		20296	1806
24hr		7707	920		168		20269	1839
72hr		7688	901		159		20316	1850

Mock 1hr	7749	883		162		20270	1854
Mock 3hr	7711	873		159		20321	1829
Mock 24hr	7671	984		172		20237	1856
Mock 72hr	7739	898		177		20250	1843
0.1gy 1hr	7733	953		173		20205	1822
0.1gy 3hr	7663	880		170		20351	1834
0.1gy 24hr	7595	964		175		20330	1825
0.1gy 72hr	7669	916		181		20298	1860
1gy 1hr		7674	916		167		20307	1819
1gy 3hr		7590	936		163		20375	1821
1gy 24hr	7593	912		177		20382	1895
1gy 72hr	7371	969		177		20547	1883

Close enough. of the ones in common event-based picks up more, although in total it's actually less. Might not report this piece
*/

data new_dtct_flags2;
  set new_dtct_flags;
  array change _numeric_;
    do over change;
    if change=. then change=0;
    end;
run;



proc freq data=new_dtct_flags2 noprint;
  tables  flag_gene_on_mock_cpm0*flag_gene_on_01gy_cpm0*flag_gene_on_1gy_cpm0 / out=freqs;
proc print data=freqs;
run;

proc freq data=new_dtct_flags2 noprint;
  tables  flag_gene_on_1hr_cpm0*flag_gene_on_3hr_cpm0*flag_gene_on_24hr_cpm0*flag_gene_on_72hr_cpm0 / out=freqs;
proc print data=freqs;
run;


proc freq data=new_dtct_flags2 noprint;
 tables flag_gene_on_mock_1hr_cpm0*flag_gene_on_mock_3hr_cpm0*flag_gene_on_mock_24hr_cpm0*
          flag_gene_on_mock_72hr_cpm0*flag_gene_on_01gy_1hr_cpm0*
          flag_gene_on_01gy_3hr_cpm0*flag_gene_on_01gy_24hr_cpm0*
          flag_gene_on_01gy_72hr_cpm0*flag_gene_on_1gy_1hr_cpm0*
          flag_gene_on_1gy_3hr_cpm0*flag_gene_on_1gy_24hr_cpm0*
          flag_gene_on_1gy_72hr_cpm0 / out=freqs;
proc print data=freqs;
run;


/*
flag_gene_    flag_gene_    flag_gene_
 on_mock_      on_01gy_       on_1gy_
   cpm0          cpm0          cpm0       COUNT    PERCENT

     0             0             0        11088    32.6338
     0             0             1          197     0.5798
     0             1             0          141     0.4150
     0             1             1          164     0.4827
     1             0             0          146     0.4297
     1             0             1          168     0.4945
     1             1             0           96     0.2825
     1             1             1        21977    64.6820

 flag_gene_    flag_gene_    flag_gene_    flag_gene_
   on_1hr_       on_3hr_      on_24hr_      on_72hr_
    cpm0          cpm0          cpm0          cpm0       COUNT    PERCENT

      0             0             0             0        10975    32.3013
      0             0             0             1          154     0.4532
      0             0             1             0          133     0.3914
      0             0             1             1           76     0.2237
      0             1             0             0          138     0.4062
      0             1             0             1           95     0.2796
      0             1             1             0           68     0.2001
      0             1             1             1          153     0.4503
      1             0             0             0          119     0.3502
      1             0             0             1           62     0.1825
      1             0             1             0           78     0.2296
      1             0             1             1          118     0.3473
      1             1             0             0           47     0.1383
      1             1             0             1          111     0.3267
      1             1             1             0           94     0.2767
      1             1             1             1        21556    63.4429



1372 rows
20787 genes on in all conditions
9732 genes off in all conditions
*/

proc export data=freqs outfile="!PATCON/arabidopsis/reports/genes_on_by_trt_time_cpm.csv"
   dbms=csv replace;
run;






/* (2) Genes DE (P<0.05, FDR P<0.05) each comparison */

data old_flags;
  set rs.fdr_by_gene_v2;
   if p_treatment < 0.05  and p_treatment > 0 then flag_treatment_p05_old=1; else flag_treatment_p05_old=0;
   if p_time < 0.05 and p_time > 0  then flag_time_p05_old=1; else flag_time_p05_old=0;
   if p_trt_by_time < 0.05 and p_trt_by_time > 0  then flag_trt_by_time_p05_old=1; else flag_trt_by_time_p05_old=0;

   if p_01gy_v_mock_1h < 0.05 and p_01gy_v_mock_1h > 0 then flag_01gy_v_mock_1h_p05_old=1; else flag_01gy_v_mock_1h_p05_old=0;
   if p_01gy_v_mock_3h < 0.05 and p_01gy_v_mock_3h > 0 then flag_01gy_v_mock_3h_p05_old=1; else flag_01gy_v_mock_3h_p05_old=0;
   if p_01gy_v_mock_24h < 0.05 and p_01gy_v_mock_24h > 0 then flag_01gy_v_mock_24h_p05_old=1; else flag_01gy_v_mock_24h_p05_old=0;
   if p_01gy_v_mock_72h < 0.05 and p_01gy_v_mock_72h > 0 then flag_01gy_v_mock_72h_p05_old=1; else flag_01gy_v_mock_72h_p05_old=0;

   if p_1gy_v_mock_1h < 0.05 and p_1gy_v_mock_1h > 0 then flag_1gy_v_mock_1h_p05_old=1; else flag_1gy_v_mock_1h_p05_old=0;
   if p_1gy_v_mock_3h < 0.05 and p_1gy_v_mock_3h > 0 then flag_1gy_v_mock_3h_p05_old=1; else flag_1gy_v_mock_3h_p05_old=0;
   if p_1gy_v_mock_24h < 0.05 and p_1gy_v_mock_24h > 0 then flag_1gy_v_mock_24h_p05_old=1; else flag_1gy_v_mock_24h_p05_old=0;
   if p_1gy_v_mock_72h < 0.05 and p_1gy_v_mock_72h > 0 then flag_1gy_v_mock_72h_p05_old=1; else flag_1gy_v_mock_72h_p05_old=0;

   drop p_: fdr_: ;
   rename flag_treatment_fdr05=flag_treatment_fdr05_old
          flag_time_fdr05=flag_time_fdr05_old
          flag_trt_by_time_fdr05=flag_trt_by_time_fdr05_old

          flag_01gy_v_mock_1h_fdr05=flag_p_01gy_v_mock_1h_p05_old
          flag_01gy_v_mock_3h_fdr05=flag_p_01gy_v_mock_3h_p05_old
          flag_01gy_v_mock_24h_fdr05=flag_p_01gy_v_mock_24h_p05_old
          flag_01gy_v_mock_72h_fdr05=flag_p_01gy_v_mock_72h_p05_old

          flag_1gy_v_mock_1h_fdr05=flag_p_1gy_v_mock_1h_p05_old
          flag_1gy_v_mock_3h_fdr05=flag_p_1gy_v_mock_3h_p05_old
          flag_1gy_v_mock_24h_fdr05=flag_p_1gy_v_mock_24h_p05_old
          flag_1gy_v_mock_72h_fdr05=flag_p_1gy_v_mock_72h_p05_old ;
run;


data new_flags;
   set rs.fdr_by_gene_log_cpm;
   if p_treatment < 0.05  and p_treatment > 0 then flag_treatment_p05_new=1; else flag_treatment_p05_new=0;
   if p_time < 0.05 and p_time > 0  then flag_time_p05_new=1; else flag_time_p05_new=0;
   if p_trt_by_time < 0.05 and p_trt_by_time > 0  then flag_trt_by_time_p05_new=1; else flag_trt_by_time_p05_new=0;

   if p_01gy_v_mock_1h < 0.05 and p_01gy_v_mock_1h > 0 then flag_01gy_v_mock_1h_p05_new=1; else flag_01gy_v_mock_1h_p05_new=0;
   if p_01gy_v_mock_3h < 0.05 and p_01gy_v_mock_3h > 0 then flag_01gy_v_mock_3h_p05_new=1; else flag_01gy_v_mock_3h_p05_new=0;
   if p_01gy_v_mock_24h < 0.05 and p_01gy_v_mock_24h > 0 then flag_01gy_v_mock_24h_p05_new=1; else flag_01gy_v_mock_24h_p05_new=0;
   if p_01gy_v_mock_72h < 0.05 and p_01gy_v_mock_72h > 0 then flag_01gy_v_mock_72h_p05_new=1; else flag_01gy_v_mock_72h_p05_new=0;

   if p_1gy_v_mock_1h < 0.05 and p_1gy_v_mock_1h > 0 then flag_1gy_v_mock_1h_p05_new=1; else flag_1gy_v_mock_1h_p05_new=0;
   if p_1gy_v_mock_3h < 0.05 and p_1gy_v_mock_3h > 0 then flag_1gy_v_mock_3h_p05_new=1; else flag_1gy_v_mock_3h_p05_new=0;
   if p_1gy_v_mock_24h < 0.05 and p_1gy_v_mock_24h > 0 then flag_1gy_v_mock_24h_p05_new=1; else flag_1gy_v_mock_24h_p05_new=0;
   if p_1gy_v_mock_72h < 0.05 and p_1gy_v_mock_72h > 0 then flag_1gy_v_mock_72h_p05_new=1; else flag_1gy_v_mock_72h_p05_new=0;

   drop p_: fdr_: ;
   rename flag_treatment_fdr05=flag_treatment_fdr05_new
          flag_time_fdr05=flag_time_fdr05_new
          flag_trt_by_time_fdr05=flag_trt_by_time_fdr05_new

          flag_01gy_v_mock_1h_fdr05=flag_p_01gy_v_mock_1h_p05_new
          flag_01gy_v_mock_3h_fdr05=flag_p_01gy_v_mock_3h_p05_new
          flag_01gy_v_mock_24h_fdr05=flag_p_01gy_v_mock_24h_p05_new
          flag_01gy_v_mock_72h_fdr05=flag_p_01gy_v_mock_72h_p05_new

          flag_1gy_v_mock_1h_fdr05=flag_p_1gy_v_mock_1h_p05_new
          flag_1gy_v_mock_3h_fdr05=flag_p_1gy_v_mock_3h_p05_new
          flag_1gy_v_mock_24h_fdr05=flag_p_1gy_v_mock_24h_p05_new
          flag_1gy_v_mock_72h_fdr05=flag_p_1gy_v_mock_72h_p05_new ;
run;

proc sort data=old_flags;
   by gene_id;
proc sort data=new_flags;
   by gene_id;
run;

data old_v_new_flags_compare new_only;
   merge old_flags (in=in1) new_flags (in=in2);
   by gene_id;
   if in1 and in2 then output old_v_new_flags_compare;
   else if in2 then output new_only;
run;

data old_v_new_flags_compare2;
   set old_v_new_flags_compare;
   array change _numeric_ ;
   do over change;
     if change=. then change=0;
     end;
run;


data new_only2;
   set new_only;
   array change _numeric_ ;
   do over change;
     if change=. then change=0;
     end;
run;

proc freq data=old_v_new_flags_compare2;
   tables 
flag_treatment_p05_new*flag_treatment_p05_old
flag_time_p05_new*flag_time_p05_old
flag_trt_by_time_p05_new*flag_trt_by_time_p05_old
flag_01gy_v_mock_1h_p05_new*flag_01gy_v_mock_1h_p05_old
flag_01gy_v_mock_3h_p05_new*flag_01gy_v_mock_3h_p05_old
flag_01gy_v_mock_24h_p05_new*flag_01gy_v_mock_24h_p05_old
flag_01gy_v_mock_72h_p05_new*flag_01gy_v_mock_72h_p05_old
flag_1gy_v_mock_1h_p05_new*flag_1gy_v_mock_1h_p05_old
flag_1gy_v_mock_3h_p05_new*flag_1gy_v_mock_3h_p05_old
flag_1gy_v_mock_24h_p05_new*flag_1gy_v_mock_24h_p05_old
flag_1gy_v_mock_72h_p05_new*flag_1gy_v_mock_72h_p05_old

flag_treatment_fdr05_new*flag_treatment_fdr05_old
flag_time_fdr05_new*flag_time_fdr05_old
flag_trt_by_time_fdr05_new*flag_trt_by_time_fdr05_old
flag_p_01gy_v_mock_1h_p05_new*flag_p_01gy_v_mock_1h_p05_old
flag_p_01gy_v_mock_3h_p05_new*flag_p_01gy_v_mock_3h_p05_old
flag_p_01gy_v_mock_24h_p05_new*flag_p_01gy_v_mock_24h_p05_old
flag_p_01gy_v_mock_72h_p05_new*flag_p_01gy_v_mock_72h_p05_old
flag_p_1gy_v_mock_1h_p05_new*flag_p_1gy_v_mock_1h_p05_old
flag_p_1gy_v_mock_3h_p05_new*flag_p_1gy_v_mock_3h_p05_old
flag_p_1gy_v_mock_24h_p05_new*flag_p_1gy_v_mock_24h_p05_old
flag_p_1gy_v_mock_72h_p05_new*flag_p_1gy_v_mock_72h_p05_old ;
run;


P<0.05
CONTRAST	NONE	APN0	CPM0	BOTH
0.1vM 1h	13239	1019	956		5514
0.1vM 3h	17040	760		778		2150
0.1vM 24h	16121	1025	1670	1912
0.1vM 72h	15468	899		1058	3303
1vM 1h		14646	835		914		4333
1vM 3h		14089	1357	2654	2628
1vM 24h		17293	816		872		1747
1vM 72h		16596	845		987		2300

FDR<0.05
CONTRAST	NONE	APN0	CPM0	BOTH
0.1vM 1h	15750	635		757		3586
0.1vM 3h	20576	49		25		78
0.1vM 24h	20132	117		288		191
0.1vM 72h	18893	353		411		1071
1vM 1h		17335	392		567		2434
1vM 3h		18473	248		1482	525
1vM 24h		20664	19		26		19
1vM 72h		20585	13		90		40

*/



proc freq data= new_only2;
   tables 
flag_01gy_v_mock_1h_p05_new
flag_01gy_v_mock_3h_p05_new
flag_01gy_v_mock_24h_p05_new
flag_01gy_v_mock_72h_p05_new
flag_1gy_v_mock_1h_p05_new
flag_1gy_v_mock_3h_p05_new
flag_1gy_v_mock_24h_p05_new
flag_1gy_v_mock_72h_p05_new

flag_p_01gy_v_mock_1h_p05_new
flag_p_01gy_v_mock_3h_p05_new
flag_p_01gy_v_mock_24h_p05_new
flag_p_01gy_v_mock_72h_p05_new
flag_p_1gy_v_mock_1h_p05_new
flag_p_1gy_v_mock_3h_p05_new
flag_p_1gy_v_mock_24h_p05_new
flag_p_1gy_v_mock_72h_p05_new;
run;




/*

      flag_01gy_v_                             Cumulative    Cumulative
   mock_1h_p05_new    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        1788       82.74          1788        82.74
                 1         373       17.26          2161       100.00


      flag_01gy_v_                             Cumulative    Cumulative
   mock_3h_p05_new    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        1952       90.33          1952        90.33
                 1         209        9.67          2161       100.00


      flag_01gy_v_                             Cumulative    Cumulative
  mock_24h_p05_new    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        1951       90.28          1951        90.28
                 1         210        9.72          2161       100.00


     flag_01gy_v_                             Cumulative    Cumulative
 mock_72h_p05_new    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        1945       90.00          1945        90.00
                1         216       10.00          2161       100.00


 flag_1gy_v_mock_                             Cumulative    Cumulative
       1h_p05_new    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        1830       84.68          1830        84.68
                1         331       15.32          2161       100.00


 flag_1gy_v_mock_                             Cumulative    Cumulative
       3h_p05_new    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        1855       85.84          1855        85.84
                1         306       14.16          2161       100.00


  flag_1gy_v_mock_                             Cumulative    Cumulative
       24h_p05_new    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        1997       92.41          1997        92.41
                 1         164        7.59          2161       100.00


  flag_1gy_v_mock_                             Cumulative    Cumulative
       72h_p05_new    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        1976       91.44          1976        91.44
                 1         185        8.56          2161       100.00






    flag_p_01gy_v_                             Cumulative    Cumulative
   mock_1h_p05_new    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        1938       89.68          1938        89.68
                 1         223       10.32          2161       100.00


     flag_p_01gy_v_                             Cumulative    Cumulative
    mock_3h_p05_new    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        2154       99.68          2154        99.68
                  1           7        0.32          2161       100.00


     flag_p_01gy_v_                             Cumulative    Cumulative
   mock_24h_p05_new    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        2136       98.84          2136        98.84
                  1          25        1.16          2161       100.00


     flag_p_01gy_v_                             Cumulative    Cumulative
   mock_72h_p05_new    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        2098       97.08          2098        97.08
                  1          63        2.92          2161       100.00


       flag_p_1gy_v_                             Cumulative    Cumulative
     mock_1h_p05_new    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0        1986       91.90          1986        91.90
                   1         175        8.10          2161       100.00


       flag_p_1gy_v_                             Cumulative    Cumulative
     mock_3h_p05_new    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0        2073       95.93          2073        95.93
                   1          88        4.07          2161       100.00


       flag_p_1gy_v_                             Cumulative    Cumulative
    mock_24h_p05_new    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0        2159       99.91          2159        99.91
                   1           2        0.09          2161       100.00

    flag_p_1gy_v_                             Cumulative    Cumulative
 mock_72h_p05_new    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        2156       99.77          2156        99.77
                1           5        0.23          2161       100.00


*/





/* Look at up down indicators. Might just do the "final" U/D/etc business */

data old_sign;
  set rs.arab_sign_by_contrast_gene;
  keep gene_id sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
               sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h;
  rename sign_01gy_v_Mock_1h =sign_01gy_v_Mock_1h_old
               sign_01gy_v_Mock_3h =sign_01gy_v_Mock_3h_old
               sign_01gy_v_Mock_24h=sign_01gy_v_Mock_24h_old 
               sign_01gy_v_Mock_72h=sign_01gy_v_Mock_72h_old
               sign_1gy_v_Mock_1h=sign_1gy_v_Mock_1h_old
               sign_1gy_v_Mock_3h=sign_1gy_v_Mock_3h_old
               sign_1gy_v_Mock_24h=sign_1gy_v_Mock_24h_old
               sign_1gy_v_Mock_72h=sign_1gy_v_Mock_72h_old;
run;

data new_sign;
  set rs.arab_sign_by_contrast_gene_cpm;
  keep gene_id sign_01gy_v_Mock_1h sign_01gy_v_Mock_3h sign_01gy_v_Mock_24h sign_01gy_v_Mock_72h
               sign_1gy_v_Mock_1h sign_1gy_v_Mock_3h sign_1gy_v_Mock_24h sign_1gy_v_Mock_72h;
  rename sign_01gy_v_Mock_1h =sign_01gy_v_Mock_1h_new
               sign_01gy_v_Mock_3h =sign_01gy_v_Mock_3h_new
               sign_01gy_v_Mock_24h=sign_01gy_v_Mock_24h_new 
               sign_01gy_v_Mock_72h=sign_01gy_v_Mock_72h_new
               sign_1gy_v_Mock_1h=sign_1gy_v_Mock_1h_new
               sign_1gy_v_Mock_3h=sign_1gy_v_Mock_3h_new
               sign_1gy_v_Mock_24h=sign_1gy_v_Mock_24h_new
               sign_1gy_v_Mock_72h=sign_1gy_v_Mock_72h_new;
run;

proc sort data=old_sign;
  by gene_id;
proc sort data=new_sign;
  by gene_id;
run;

data old_v_new_sign;
  merge old_sign (in=in1) new_sign (in=in2);
  by gene_id;
  if in1 and in2;
run;

data old_v_new_sign2;
 set old_v_new_sign;
 array change _character_;
 do over change;
   if change ="" then change="0";
 end;
run;

proc freq data=old_v_new_sign2 noprint;
  tables sign_01gy_v_Mock_1h_old*sign_01gy_v_Mock_1h_new /out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables sign_01gy_v_Mock_3h_old*sign_01gy_v_Mock_3h_new/out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables 
         sign_01gy_v_Mock_24h_old*sign_01gy_v_Mock_24h_new/out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables 
         sign_01gy_v_Mock_72h_old*sign_01gy_v_Mock_72h_new/out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables 
         sign_1gy_v_Mock_1h_old*sign_1gy_v_Mock_1h_new/out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables 
         sign_1gy_v_Mock_3h_old*sign_1gy_v_Mock_3h_new/out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables 
         sign_1gy_v_Mock_24h_old*sign_1gy_v_Mock_24h_new/out=freqs;
proc print data=freqs;
run;
proc freq data=old_v_new_sign2 noprint;
  tables 
         sign_1gy_v_Mock_72h_old*sign_1gy_v_Mock_72h_new/out=freqs;
proc print data=freqs;
run;




/*
sign_01gy_    sign_01gy_
 v_Mock_       v_Mock_
  1h_old        1h_new      COUNT

    D             D          2800
    D             N           798
    N             D           276
    N             N         13152
    N             U           703
    U             N           243
    U             U          2756

sign_01gy_    sign_01gy_
 v_Mock_       v_Mock_
  3h_old        3h_new      COUNT

    D             D          1175
    D             N           516
    N             D           508
    N             N         16976
    N             U           289
    U             N           257
    U             U          1007

sign_01gy_    sign_01gy_
 v_Mock_       v_Mock_
 24h_old       24h_new      COUNT

    D             D           624
    D             N           824
    D             U             1
    N             D           110
    N             N         16024
    N             U          1605
    U             N           210
    U             U          1330

 sign_01gy_    sign_01gy_
  v_Mock_       v_Mock_
  72h_old       72h_new      COUNT

     D             D          1720
     D             N           665
     N             D           258
     N             N         15403
     N             U           824
     U             N           249
     U             U          1609

 sign_1gy_    sign_1gy_
  v_Mock_      v_Mock_
  1h_old       1h_new      COUNT

     D            D         1963
     D            N          575
     N            D          258
     N            N        14590
     N            U          672
     U            D            2
     U            N          274
     U            U         2394

   

sign_1gy_    sign_1gy_
 v_Mock_      v_Mock_
 3h_old       3h_new      COUNT

    D            D          740
    D            N         1226
    D            U            1
    N            D           83
    N            N        14018
    N            U         2599
    U            N          144
    U            U         1917

sign_1gy_    sign_1gy_
 v_Mock_      v_Mock_
 24h_old      24h_new     COUNT

    D            D          970
    D            N          349
    N            D          623
    N            N        17237
    N            U          276
    U            N          475
    U            U          798

sign_1gy_    sign_1gy_
 v_Mock_      v_Mock_
 72h_old      72h_new     COUNT

    D            D          664
    D            N          630
    N            D          146
    N            N        16479
    N            U          867
    U            N          230
    U            U         1712
*/



data old_sign;
  set rs.arab_sign_by_contrast_gene_fdr2;
  keep gene_id sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_72h_fdr05
               sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_72h_fdr05;
  rename sign_01gy_v_Mock_1h_fdr05 =sign_01gy_v_Mock_1h_old
               sign_01gy_v_Mock_3h_fdr05 =sign_01gy_v_Mock_3h_old
               sign_01gy_v_Mock_24h_fdr05=sign_01gy_v_Mock_24h_old 
               sign_01gy_v_Mock_72h_fdr05=sign_01gy_v_Mock_72h_old
               sign_1gy_v_Mock_1h_fdr05=sign_1gy_v_Mock_1h_old
               sign_1gy_v_Mock_3h_fdr05=sign_1gy_v_Mock_3h_old
               sign_1gy_v_Mock_24h_fdr05=sign_1gy_v_Mock_24h_old
               sign_1gy_v_Mock_72h_fdr05=sign_1gy_v_Mock_72h_old;
run;

data new_sign;
  set rs.arab_sign_by_cntrst_gene_cpm_fdr;
  keep gene_id sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_72h_fdr05
               sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_72h_fdr05;
  rename sign_01gy_v_Mock_1h_fdr05 =sign_01gy_v_Mock_1h_new
               sign_01gy_v_Mock_3h_fdr05 =sign_01gy_v_Mock_3h_new
               sign_01gy_v_Mock_24h_fdr05=sign_01gy_v_Mock_24h_new 
               sign_01gy_v_Mock_72h_fdr05=sign_01gy_v_Mock_72h_new
               sign_1gy_v_Mock_1h_fdr05=sign_1gy_v_Mock_1h_new
               sign_1gy_v_Mock_3h_fdr05=sign_1gy_v_Mock_3h_new
               sign_1gy_v_Mock_24h_fdr05=sign_1gy_v_Mock_24h_new
               sign_1gy_v_Mock_72h_fdr05=sign_1gy_v_Mock_72h_new;
run;

proc sort data=old_sign;
  by gene_id;
proc sort data=new_sign;
  by gene_id;
run;

data old_v_new_sign;
  merge old_sign new_sign;
  by gene_id;
run;

data old_v_new_sign2;
 set old_v_new_sign;
 array change _character_;
 do over change;
   if change ="" then change="0";
 end;
run;

proc freq data=old_v_new_sign2;
  tables sign_01gy_v_Mock_1h_old*sign_01gy_v_Mock_1h_new
         sign_01gy_v_Mock_3h_old*sign_01gy_v_Mock_3h_new
         sign_01gy_v_Mock_24h_old*sign_01gy_v_Mock_24h_new
         sign_01gy_v_Mock_72h_old*sign_01gy_v_Mock_72h_new
         sign_1gy_v_Mock_1h_old*sign_1gy_v_Mock_1h_new
         sign_1gy_v_Mock_3h_old*sign_1gy_v_Mock_3h_new
         sign_1gy_v_Mock_24h_old*sign_1gy_v_Mock_24h_new
         sign_1gy_v_Mock_72h_old*sign_1gy_v_Mock_72h_new;
run;


/*

sign_01gy_v_Mock_1h_old     sign_01gy_v_Mock_1h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |0.1gy   |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    463 |    130 |     87 |    152 |   1559 |    136 |   2527
         |   1.94 |   0.54 |   0.36 |   0.64 |   6.52 |   0.57 |  10.57
         |  18.32 |   5.14 |   3.44 |   6.02 |  61.69 |   5.38 |
         |  32.58 |  34.21 |   4.23 |  31.21 |   9.15 |   5.42 |
---------+--------+--------+--------+--------+--------+--------+
0.1gy    |    193 |    141 |      0 |      4 |     14 |      1 |    353
         |   0.81 |   0.59 |   0.00 |   0.02 |   0.06 |   0.00 |   1.48
         |  54.67 |  39.94 |   0.00 |   1.13 |   3.97 |   0.28 |
         |  13.58 |  37.11 |   0.00 |   0.82 |   0.08 |   0.04 |
---------+--------+--------+--------+--------+--------+--------+
D        |      9 |      1 |   1764 |      2 |    499 |      0 |   2275
         |   0.04 |   0.00 |   7.38 |   0.01 |   2.09 |   0.00 |   9.52
         |   0.40 |   0.04 |  77.54 |   0.09 |  21.93 |   0.00 |
         |   0.63 |   0.26 |  85.71 |   0.41 |   2.93 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    195 |      0 |      0 |    200 |     25 |      0 |    420
         |   0.82 |   0.00 |   0.00 |   0.84 |   0.10 |   0.00 |   1.76
         |  46.43 |   0.00 |   0.00 |  47.62 |   5.95 |   0.00 |
         |  13.72 |   0.00 |   0.00 |  41.07 |   0.15 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    553 |     97 |    207 |    129 |  14829 |    549 |  16364
         |   2.31 |   0.41 |   0.87 |   0.54 |  62.05 |   2.30 |  68.47
         |   3.38 |   0.59 |   1.26 |   0.79 |  90.62 |   3.35 |
         |  38.92 |  25.53 |  10.06 |  26.49 |  87.00 |  21.89 |
---------+--------+--------+--------+--------+--------+--------+
U        |      8 |     11 |      0 |      0 |    119 |   1822 |   1960
         |   0.03 |   0.05 |   0.00 |   0.00 |   0.50 |   7.62 |   8.20
         |   0.41 |   0.56 |   0.00 |   0.00 |   6.07 |  92.96 |
         |   0.56 |   2.89 |   0.00 |   0.00 |   0.70 |  72.65 |
---------+--------+--------+--------+--------+--------+--------+
Total        1421      380     2058      487    17045     2508    23899
             5.95     1.59     8.61     2.04    71.32    10.49   100.00


sign_01gy_v_Mock_3h_old     sign_01gy_v_Mock_3h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |0.1gy   |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    463 |    121 |      6 |    106 |   1804 |      1 |   2501
         |   1.94 |   0.51 |   0.03 |   0.44 |   7.55 |   0.00 |  10.46
         |  18.51 |   4.84 |   0.24 |   4.24 |  72.13 |   0.04 |
         |  33.43 |  30.40 |  10.34 |  28.88 |   8.34 |   1.92 |
---------+--------+--------+--------+--------+--------+--------+
0.1gy    |    160 |    195 |      0 |      0 |     11 |      0 |    366
         |   0.67 |   0.82 |   0.00 |   0.00 |   0.05 |   0.00 |   1.53
         |  43.72 |  53.28 |   0.00 |   0.00 |   3.01 |   0.00 |
         |  11.55 |  48.99 |   0.00 |   0.00 |   0.05 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      1 |      0 |     40 |      0 |     38 |      0 |     79
         |   0.00 |   0.00 |   0.17 |   0.00 |   0.16 |   0.00 |   0.33
         |   1.27 |   0.00 |  50.63 |   0.00 |  48.10 |   0.00 |
         |   0.07 |   0.00 |  68.97 |   0.00 |   0.18 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    181 |      2 |      0 |    146 |     19 |      0 |    348
         |   0.76 |   0.01 |   0.00 |   0.61 |   0.08 |   0.00 |   1.46
         |  52.01 |   0.57 |   0.00 |  41.95 |   5.46 |   0.00 |
         |  13.07 |   0.50 |   0.00 |  39.78 |   0.09 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    578 |     80 |     12 |    115 |  19757 |     13 |  20555
         |   2.42 |   0.33 |   0.05 |   0.48 |  82.67 |   0.05 |  86.01
         |   2.81 |   0.39 |   0.06 |   0.56 |  96.12 |   0.06 |
         |  41.73 |  20.10 |  20.69 |  31.34 |  91.30 |  25.00 |
---------+--------+--------+--------+--------+--------+--------+
U        |      2 |      0 |      0 |      0 |     10 |     38 |     50
         |   0.01 |   0.00 |   0.00 |   0.00 |   0.04 |   0.16 |   0.21
         |   4.00 |   0.00 |   0.00 |   0.00 |  20.00 |  76.00 |
         |   0.14 |   0.00 |   0.00 |   0.00 |   0.05 |  73.08 |
---------+--------+--------+--------+--------+--------+--------+
Total        1385      398       58      367    21639       52    23899
             5.80     1.67     0.24     1.54    90.54     0.22   100.00


sign_01gy_v_Mock_24h_old     sign_01gy_v_Mock_24h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |0.1gy   |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    446 |    105 |      8 |    135 |   1785 |     17 |   2496
         |   1.87 |   0.44 |   0.03 |   0.56 |   7.47 |   0.07 |  10.44
         |  17.87 |   4.21 |   0.32 |   5.41 |  71.51 |   0.68 |
         |  32.13 |  24.36 |   7.27 |  35.43 |   8.42 |   4.31 |
---------+--------+--------+--------+--------+--------+--------+
0.1gy    |    175 |    184 |      0 |      0 |     18 |      0 |    377
         |   0.73 |   0.77 |   0.00 |   0.00 |   0.08 |   0.00 |   1.58
         |  46.42 |  48.81 |   0.00 |   0.00 |   4.77 |   0.00 |
         |  12.61 |  42.69 |   0.00 |   0.00 |   0.08 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      1 |      0 |     75 |      0 |     94 |      0 |    170
         |   0.00 |   0.00 |   0.31 |   0.00 |   0.39 |   0.00 |   0.71
         |   0.59 |   0.00 |  44.12 |   0.00 |  55.29 |   0.00 |
         |   0.07 |   0.00 |  68.18 |   0.00 |   0.44 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    173 |      3 |      0 |    125 |     16 |      0 |    317
         |   0.72 |   0.01 |   0.00 |   0.52 |   0.07 |   0.00 |   1.33
         |  54.57 |   0.95 |   0.00 |  39.43 |   5.05 |   0.00 |
         |  12.46 |   0.70 |   0.00 |  32.81 |   0.08 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    592 |    138 |     27 |    121 |  19260 |    261 |  20399
         |   2.48 |   0.58 |   0.11 |   0.51 |  80.59 |   1.09 |  85.36
         |   2.90 |   0.68 |   0.13 |   0.59 |  94.42 |   1.28 |
         |  42.65 |  32.02 |  24.55 |  31.76 |  90.87 |  66.24 |
---------+--------+--------+--------+--------+--------+--------+
U        |      1 |      1 |      0 |      0 |     22 |    116 |    140
         |   0.00 |   0.00 |   0.00 |   0.00 |   0.09 |   0.49 |   0.59
         |   0.71 |   0.71 |   0.00 |   0.00 |  15.71 |  82.86 |
         |   0.07 |   0.23 |   0.00 |   0.00 |   0.10 |  29.44 |
---------+--------+--------+--------+--------+--------+--------+
Total        1388      431      110      381    21195      394    23899
             5.81     1.80     0.46     1.59    88.69     1.65   100.00


sign_01gy_v_Mock_72h_old     sign_01gy_v_Mock_72h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |0.1gy   |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    799 |     52 |     28 |    107 |   1784 |     35 |   2805
         |   3.34 |   0.22 |   0.12 |   0.45 |   7.46 |   0.15 |  11.74
         |  28.48 |   1.85 |   1.00 |   3.81 |  63.60 |   1.25 |
         |  46.67 |  50.00 |   3.56 |  29.16 |   8.84 |   4.62 |
---------+--------+--------+--------+--------+--------+--------+
0.1gy    |     57 |     38 |      0 |      0 |      2 |      0 |     97
         |   0.24 |   0.16 |   0.00 |   0.00 |   0.01 |   0.00 |   0.41
         |  58.76 |  39.18 |   0.00 |   0.00 |   2.06 |   0.00 |
         |   3.33 |  36.54 |   0.00 |   0.00 |   0.01 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      2 |      0 |    614 |      0 |    273 |      0 |    889
         |   0.01 |   0.00 |   2.57 |   0.00 |   1.14 |   0.00 |   3.72
         |   0.22 |   0.00 |  69.07 |   0.00 |  30.71 |   0.00 |
         |   0.12 |   0.00 |  78.02 |   0.00 |   1.35 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    149 |      0 |      0 |    140 |     20 |      0 |    309
         |   0.62 |   0.00 |   0.00 |   0.59 |   0.08 |   0.00 |   1.29
         |  48.22 |   0.00 |   0.00 |  45.31 |   6.47 |   0.00 |
         |   8.70 |   0.00 |   0.00 |  38.15 |   0.10 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    702 |     13 |    145 |    120 |  18015 |    266 |  19261
         |   2.94 |   0.05 |   0.61 |   0.50 |  75.38 |   1.11 |  80.59
         |   3.64 |   0.07 |   0.75 |   0.62 |  93.53 |   1.38 |
         |  41.00 |  12.50 |  18.42 |  32.70 |  89.31 |  35.09 |
---------+--------+--------+--------+--------+--------+--------+
U        |      3 |      1 |      0 |      0 |     77 |    457 |    538
         |   0.01 |   0.00 |   0.00 |   0.00 |   0.32 |   1.91 |   2.25
         |   0.56 |   0.19 |   0.00 |   0.00 |  14.31 |  84.94 |
         |   0.18 |   0.96 |   0.00 |   0.00 |   0.38 |  60.29 |
---------+--------+--------+--------+--------+--------+--------+
Total        1712      104      787      367    20171      758    23899
             7.16     0.44     3.29     1.54    84.40     3.17   100.00


sign_1gy_v_Mock_1h_old     sign_1gy_v_Mock_1h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |1gy     |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    463 |    115 |     84 |    128 |   1628 |     91 |   2509
         |   1.94 |   0.48 |   0.35 |   0.54 |   6.81 |   0.38 |  10.50
         |  18.45 |   4.58 |   3.35 |   5.10 |  64.89 |   3.63 |
         |  32.91 |  29.19 |   6.54 |  33.07 |   8.78 |   4.81 |
---------+--------+--------+--------+--------+--------+--------+
1gy      |    181 |    168 |      0 |      3 |     19 |      0 |    371
         |   0.76 |   0.70 |   0.00 |   0.01 |   0.08 |   0.00 |   1.55
         |  48.79 |  45.28 |   0.00 |   0.81 |   5.12 |   0.00 |
         |  12.86 |  42.64 |   0.00 |   0.78 |   0.10 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      4 |      0 |   1044 |      2 |    280 |      0 |   1330
         |   0.02 |   0.00 |   4.37 |   0.01 |   1.17 |   0.00 |   5.57
         |   0.30 |   0.00 |  78.50 |   0.15 |  21.05 |   0.00 |
         |   0.28 |   0.00 |  81.25 |   0.52 |   1.51 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    185 |      1 |      0 |    138 |     19 |      0 |    343
         |   0.77 |   0.00 |   0.00 |   0.58 |   0.08 |   0.00 |   1.44
         |  53.94 |   0.29 |   0.00 |  40.23 |   5.54 |   0.00 |
         |  13.15 |   0.25 |   0.00 |  35.66 |   0.10 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    572 |    108 |    157 |    115 |  16483 |    410 |  17845
         |   2.39 |   0.45 |   0.66 |   0.48 |  68.97 |   1.72 |  74.67
         |   3.21 |   0.61 |   0.88 |   0.64 |  92.37 |   2.30 |
         |  40.65 |  27.41 |  12.22 |  29.72 |  88.93 |  21.68 |
---------+--------+--------+--------+--------+--------+--------+
U        |      2 |      2 |      0 |      1 |    106 |   1390 |   1501
         |   0.01 |   0.01 |   0.00 |   0.00 |   0.44 |   5.82 |   6.28
         |   0.13 |   0.13 |   0.00 |   0.07 |   7.06 |  92.60 |
         |   0.14 |   0.51 |   0.00 |   0.26 |   0.57 |  73.51 |
---------+--------+--------+--------+--------+--------+--------+
Total        1407      394     1285      387    18535     1891    23899
             5.89     1.65     5.38     1.62    77.56     7.91   100.00



sign_1gy_v_Mock_3h_old     sign_1gy_v_Mock_3h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |1gy     |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    434 |    127 |     32 |    127 |   1700 |     56 |   2476
         |   1.82 |   0.53 |   0.13 |   0.53 |   7.11 |   0.23 |  10.36
         |  17.53 |   5.13 |   1.29 |   5.13 |  68.66 |   2.26 |
         |  31.59 |  31.05 |  10.39 |  35.08 |   8.65 |   3.13 |
---------+--------+--------+--------+--------+--------+--------+
1gy      |    187 |    191 |      0 |      2 |     11 |      0 |    391
         |   0.78 |   0.80 |   0.00 |   0.01 |   0.05 |   0.00 |   1.64
         |  47.83 |  48.85 |   0.00 |   0.51 |   2.81 |   0.00 |
         |  13.61 |  46.70 |   0.00 |   0.55 |   0.06 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      1 |      0 |    219 |      0 |    233 |      0 |    453
         |   0.00 |   0.00 |   0.92 |   0.00 |   0.97 |   0.00 |   1.90
         |   0.22 |   0.00 |  48.34 |   0.00 |  51.43 |   0.00 |
         |   0.07 |   0.00 |  71.10 |   0.00 |   1.19 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    156 |      0 |      0 |    117 |     12 |      0 |    285
         |   0.65 |   0.00 |   0.00 |   0.49 |   0.05 |   0.00 |   1.19
         |  54.74 |   0.00 |   0.00 |  41.05 |   4.21 |   0.00 |
         |  11.35 |   0.00 |   0.00 |  32.32 |   0.06 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    594 |     90 |     57 |    116 |  17689 |   1425 |  19971
         |   2.49 |   0.38 |   0.24 |   0.49 |  74.02 |   5.96 |  83.56
         |   2.97 |   0.45 |   0.29 |   0.58 |  88.57 |   7.14 |
         |  43.23 |  22.00 |  18.51 |  32.04 |  89.98 |  79.74 |
---------+--------+--------+--------+--------+--------+--------+
U        |      2 |      1 |      0 |      0 |     14 |    306 |    323
         |   0.01 |   0.00 |   0.00 |   0.00 |   0.06 |   1.28 |   1.35
         |   0.62 |   0.31 |   0.00 |   0.00 |   4.33 |  94.74 |
         |   0.15 |   0.24 |   0.00 |   0.00 |   0.07 |  17.12 |
---------+--------+--------+--------+--------+--------+--------+
Total        1374      409      308      362    19659     1787    23899
             5.75     1.71     1.29     1.51    82.26     7.48   100.00


sign_1gy_v_Mock_24h_old     sign_1gy_v_Mock_24h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |1gy     |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    390 |    131 |      1 |    109 |   1832 |      1 |   2464
         |   1.63 |   0.55 |   0.00 |   0.46 |   7.67 |   0.00 |  10.31
         |  15.83 |   5.32 |   0.04 |   4.42 |  74.35 |   0.04 |
         |  28.74 |  28.35 |   3.23 |  36.09 |   8.43 |   6.25 |
---------+--------+--------+--------+--------+--------+--------+
1gy      |    204 |    185 |      0 |      1 |     19 |      0 |    409
         |   0.85 |   0.77 |   0.00 |   0.00 |   0.08 |   0.00 |   1.71
         |  49.88 |  45.23 |   0.00 |   0.24 |   4.65 |   0.00 |
         |  15.03 |  40.04 |   0.00 |   0.33 |   0.09 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      1 |      0 |      7 |      0 |      6 |      0 |     14
         |   0.00 |   0.00 |   0.03 |   0.00 |   0.03 |   0.00 |   0.06
         |   7.14 |   0.00 |  50.00 |   0.00 |  42.86 |   0.00 |
         |   0.07 |   0.00 |  22.58 |   0.00 |   0.03 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    160 |      2 |      0 |    112 |     19 |      0 |    293
         |   0.67 |   0.01 |   0.00 |   0.47 |   0.08 |   0.00 |   1.23
         |  54.61 |   0.68 |   0.00 |  38.23 |   6.48 |   0.00 |
         |  11.79 |   0.43 |   0.00 |  37.09 |   0.09 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    602 |    144 |     23 |     80 |  19842 |      3 |  20694
         |   2.52 |   0.60 |   0.10 |   0.33 |  83.02 |   0.01 |  86.59
         |   2.91 |   0.70 |   0.11 |   0.39 |  95.88 |   0.01 |
         |  44.36 |  31.17 |  74.19 |  26.49 |  91.31 |  18.75 |
---------+--------+--------+--------+--------+--------+--------+
U        |      0 |      0 |      0 |      0 |     13 |     12 |     25
         |   0.00 |   0.00 |   0.00 |   0.00 |   0.05 |   0.05 |   0.10
         |   0.00 |   0.00 |   0.00 |   0.00 |  52.00 |  48.00 |
         |   0.00 |   0.00 |   0.00 |   0.00 |   0.06 |  75.00 |
---------+--------+--------+--------+--------+--------+--------+
Total        1357      462       31      302    21731       16    23899
             5.68     1.93     0.13     1.26    90.93     0.07   100.00


sign_1gy_v_Mock_72h_old     sign_1gy_v_Mock_72h_new

Frequency|
Percent  |
Row Pct  |
Col Pct  |0       |1gy     |D       |Mock    |N       |U       |  Total
---------+--------+--------+--------+--------+--------+--------+
0        |    358 |    138 |      3 |     96 |   1833 |      2 |   2430
         |   1.50 |   0.58 |   0.01 |   0.40 |   7.67 |   0.01 |  10.17
         |  14.73 |   5.68 |   0.12 |   3.95 |  75.43 |   0.08 |
         |  27.58 |  26.64 |  11.11 |  34.66 |   8.46 |   1.85 |
---------+--------+--------+--------+--------+--------+--------+
1gy      |    201 |    249 |      0 |      2 |     20 |      0 |    472
         |   0.84 |   1.04 |   0.00 |   0.01 |   0.08 |   0.00 |   1.97
         |  42.58 |  52.75 |   0.00 |   0.42 |   4.24 |   0.00 |
         |  15.49 |  48.07 |   0.00 |   0.72 |   0.09 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
D        |      0 |      0 |      9 |      0 |      6 |      0 |     15
         |   0.00 |   0.00 |   0.04 |   0.00 |   0.03 |   0.00 |   0.06
         |   0.00 |   0.00 |  60.00 |   0.00 |  40.00 |   0.00 |
         |   0.00 |   0.00 |  33.33 |   0.00 |   0.03 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
Mock     |    132 |      1 |      0 |     94 |     14 |      0 |    241
         |   0.55 |   0.00 |   0.00 |   0.39 |   0.06 |   0.00 |   1.01
         |  54.77 |   0.41 |   0.00 |  39.00 |   5.81 |   0.00 |
         |  10.17 |   0.19 |   0.00 |  33.94 |   0.06 |   0.00 |
---------+--------+--------+--------+--------+--------+--------+
N        |    607 |    130 |     15 |     85 |  19791 |     75 |  20703
         |   2.54 |   0.54 |   0.06 |   0.36 |  82.81 |   0.31 |  86.63
         |   2.93 |   0.63 |   0.07 |   0.41 |  95.59 |   0.36 |
         |  46.76 |  25.10 |  55.56 |  30.69 |  91.32 |  69.44 |
---------+--------+--------+--------+--------+--------+--------+
U        |      0 |      0 |      0 |      0 |      7 |     31 |     38
         |   0.00 |   0.00 |   0.00 |   0.00 |   0.03 |   0.13 |   0.16
         |   0.00 |   0.00 |   0.00 |   0.00 |  18.42 |  81.58 |
         |   0.00 |   0.00 |   0.00 |   0.00 |   0.03 |  28.70 |
---------+--------+--------+--------+--------+--------+--------+
Total        1298      518       27      277    21671      108    23899
             5.43     2.17     0.11     1.16    90.68     0.45   100.00

*/



data new_sign;
  set rs.arab_sign_by_cntrst_gene_cpm_fdr;
  keep gene_id sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_72h_fdr05
               sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_3h_fdr05 sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_72h_fdr05;
  rename sign_01gy_v_Mock_1h_fdr05 =sign_01gy_v_Mock_1h_new
               sign_01gy_v_Mock_3h_fdr05 =sign_01gy_v_Mock_3h_new
               sign_01gy_v_Mock_24h_fdr05=sign_01gy_v_Mock_24h_new 
               sign_01gy_v_Mock_72h_fdr05=sign_01gy_v_Mock_72h_new
               sign_1gy_v_Mock_1h_fdr05=sign_1gy_v_Mock_1h_new
               sign_1gy_v_Mock_3h_fdr05=sign_1gy_v_Mock_3h_new
               sign_1gy_v_Mock_24h_fdr05=sign_1gy_v_Mock_24h_new
               sign_1gy_v_Mock_72h_fdr05=sign_1gy_v_Mock_72h_new;
run;


proc freq data=new_sign;
  tables sign_01gy_v_Mock_1h_new  sign_01gy_v_Mock_3h_new sign_01gy_v_Mock_24h_new sign_01gy_v_Mock_72h_new
         sign_1gy_v_Mock_1h_new sign_1gy_v_Mock_3h_new sign_1gy_v_Mock_24h_new sign_1gy_v_Mock_72h_new;
run;

/*
sign_01gy_
v_Mock_                                Cumulative    Cumulative
1h_new        Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              380        1.69           380         1.69
D                 2058        9.16          2438        10.85
Mock               487        2.17          2925        13.01
N                17045       75.83         19970        88.84
U                 2508       11.16         22478       100.00

                    Frequency Missing = 411


sign_01gy_
v_Mock_                                Cumulative    Cumulative
3h_new        Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              398        1.77           398         1.77
D                   58        0.26           456         2.03
Mock               367        1.63           823         3.66
N                21639       96.11         22462        99.77
U                   52        0.23         22514       100.00

                    Frequency Missing = 375




*/sign_01gy_
v_Mock_                                Cumulative    Cumulative
24h_new       Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              431        1.91           431         1.91
D                  110        0.49           541         2.40
Mock               381        1.69           922         4.10
N                21195       94.15         22117        98.25
U                  394        1.75         22511       100.00

                    Frequency Missing = 378


sign_01gy_
v_Mock_                                Cumulative    Cumulative
72h_new       Frequency     Percent     Frequency      Percent
---------------------------------------------------------------
0.1gy              104        0.47           104         0.47
D                  787        3.55           891         4.02
Mock               367        1.65          1258         5.67
N                20171       90.91         21429        96.58
U                  758        3.42         22187       100.00

                    Frequency Missing = 702

sign_1gy_
v_Mock_                               Cumulative    Cumulative
1h_new       Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               394        1.75           394         1.75
D                1285        5.71          1679         7.46
Mock              387        1.72          2066         9.19
N               18535       82.41         20601        91.59
U                1891        8.41         22492       100.00

                   Frequency Missing = 397


sign_1gy_
v_Mock_                               Cumulative    Cumulative
3h_new       Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               409        1.82           409         1.82
D                 308        1.37           717         3.18
Mock              362        1.61          1079         4.79
N               19659       87.28         20738        92.07
U                1787        7.93         22525       100.00

                   Frequency Missing = 364

sign_1gy_
v_Mock_                               Cumulative    Cumulative
24h_new      Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               462        2.05           462         2.05
D                  31        0.14           493         2.19
Mock              302        1.34           795         3.53
N               21731       96.40         22526        99.93
U                  16        0.07         22542       100.00

                   Frequency Missing = 347


sign_1gy_
v_Mock_                               Cumulative    Cumulative
72h_new      Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
1gy               518        2.29           518         2.29
D                  27        0.12           545         2.41
Mock              277        1.23           822         3.64
N               21671       95.89         22493        99.52
U                 108        0.48         22601       100.00

                   Frequency Missing = 288

*/
