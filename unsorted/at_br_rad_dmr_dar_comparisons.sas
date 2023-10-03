ods listing; ods html close;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_sig_regions.csv"
   out=at_regions dbms=csv replace;
   guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/br_rad_sig_regions.csv"
   out=br_regions dbms=csv replace;
   guessingrows=max;
run;


proc sort data=at_regions nodup;
   by _all_;
proc sort data=at_regions;
   by chr region_start region_stop;
run;

data at_regions_group;
    retain region_num superregion_start superregion_stop;
    set at_regions;
    by chr;
    prev_start=lag1(region_start);
    prev_stop=lag1(region_stop);
    if first.chr then do;
         superregion_start=region_start;
         superregion_stop=region_stop;
         region_num=1;
         end;
    else do;
         if region_start < prev_stop then do;
              region_num=region_num;
              if region_stop > superregion_stop then superregion_stop=region_stop;
              else superregion_stop=superregion_stop;
              superregion_start=superregion_start;
              end;
         else do;
              region_num=region_num + 1;
              superregion_start=region_start;
              superregion_stop=region_stop;
               end;
         end;
     if site_type="GC" then do;
      flag_dmr=0; flag_dar=1;
      end;
     else do;
      flag_dmr=1; flag_dar=0;
      end;
run;


/* Mingqi's desired table output:
   chr, region_start, region_stop, N DAR 0.1Gy 72h, N DMR 0.1Gy 72h, N DAR 1Gy 72h, N DMR 1Gy 72h 
   flag DAR 0.1Gy 72h, flag DMR 0.1Gy 72h, flag DAR 1Gy 72h, flag DMR 1Gy 72h*/


data at_regions_group2;
   set at_regions_group;
   if flag_dmr=1 and treatment="0.1Gy" then flag_DMR_01_72=1; else flag_DMR_01_72=0;
   if flag_dar=1 and treatment="0.1Gy" then flag_DAR_01_72=1; else flag_DAR_01_72=0;
   if flag_dmr=1 and treatment="1Gy" then flag_DMR_1_72=1; else flag_DMR_1_72=0;
   if flag_dar=1 and treatment="1Gy" then flag_DAR_1_72=1; else flag_DAR_1_72=0;
run;

proc sort data=at_regions_group2;
     by chr region_num;
proc means data=at_regions_group2 noprint;
     by chr region_num;
     var superregion_start superregion_stop flag_DMR_01_72 flag_DAR_01_72 flag_DMR_1_72 flag_DAR_1_72;
     output out=at_regions_collapsed(drop=_TYPE_ _FREQ_) min(superregion_start)=superregion_start
                                     max(superregion_stop)=superregion_stop
                                     sum(flag_DMR_01_72)=N_DMR_01_72 max(flag_DMR_01_72)=flag_DMR_01_72
                                     sum(flag_DAR_01_72)=N_DAR_01_72 max(flag_DAR_01_72)=flag_DAR_01_72
                                     sum(flag_DMR_1_72)=N_DMR_1_72 max(flag_DMR_1_72)=flag_DMR_1_72
                                     sum(flag_DAR_1_72)=N_DAR_1_72 max(flag_DAR_1_72)=flag_DAR_1_72;
run;


data at_Regions_collapsed2;
   retain chr superregion_start superregion_stop sum_: flag_: ;
   set at_regions_collapsed;
   rename superregion_start=region_start superregion_stop=region_stop;
run;



proc freq data=at_regions_collapsed2 noprint;
  tables flag_DAR_01_72*flag_DMR_01_72*flag_DAR_1_72*flag_DMR_1_72 / out=at_Region_cnt;
run;

proc print data=at_region_cnt;
run;

/*
   flag_      flag_
  DAR_01_    DMR_01_      flag_       flag_
     72         72      DAR_1_72    DMR_1_72    COUNT

     0          0           0           1         459
     0          0           1           0        8041
     0          0           1           1          15
     0          1           0           0         215
     0          1           0           1           9
     0          1           1           0           6
     1          0           0           0       28163
     1          0           0           1          21
     1          0           1           0       11809
     1          0           1           1           5
     1          1           0           0          14
     1          1           0           1           1
     1          1           1           0           4

*/

   

proc sort data=br_regions nodup;
   by _all_;
proc sort data=br_regions;
   by chr region_start region_stop;
run;

data br_regions_group;
    retain region_num superregion_start superregion_stop;
    set br_regions;
    by chr;
    prev_start=lag1(region_start);
    prev_stop=lag1(region_stop);
    if first.chr then do;
         superregion_start=region_start;
         superregion_stop=region_stop;
         region_num=1;
         end;
    else do;
         if region_start < prev_stop then do;
              region_num=region_num;
              if region_stop > superregion_stop then superregion_stop=region_stop;
              else superregion_stop=superregion_stop;
              superregion_start=superregion_start;
              end;
         else do;
              region_num=region_num + 1;
              superregion_start=region_start;
              superregion_stop=region_stop;
               end;
         end;
     if site_type="GC" then do;
      flag_dmr=0; flag_dar=1;
      end;
     else do;
      flag_dmr=1; flag_dar=0;
      end;
run;


/* Mingqi's desired table output:
   chr, region_start, region_stop, N DAR 0.1Gy 72h, N DMR 0.1Gy 72h, N DAR 1Gy 72h, N DMR 1Gy 72h 
   flag DAR 0.1Gy 72h, flag DMR 0.1Gy 72h, flag DAR 1Gy 72h, flag DMR 1Gy 72h*/

1.4cGy_1h
1.4cGy_72h
10cGy_1h
10cGy_72h*/

data br_regions_group2;
   set br_regions_group;
   if flag_dmr=1 and treatment="1.4cGy_1h" then flag_DMR_1p4_1=1; else flag_DMR_1p4_1=0;
   if flag_dar=1 and treatment="1.4cGy_1h" then flag_DAR_1p4_1=1; else flag_DAR_1p4_1=0;
   if flag_dmr=1 and treatment="1.4cGy_72h" then flag_DMR_1p4_72=1; else flag_DMR_1p4_72=0;
   if flag_dar=1 and treatment="1.4cGy_72h" then flag_DAR_1p4_72=1; else flag_DAR_1p4_72=0;
   if flag_dmr=1 and treatment="10cGy_1h" then flag_DMR_10_1=1; else flag_DMR_10_1=0;
   if flag_dar=1 and treatment="10cGy_1h" then flag_DAR_10_1=1; else flag_DAR_10_1=0;
   if flag_dmr=1 and treatment="10cGy_72h" then flag_DMR_10_72=1; else flag_DMR_10_72=0;
   if flag_dar=1 and treatment="10cGy_72h" then flag_DAR_10_72=1; else flag_DAR_10_72=0;
run;

proc sort data=br_regions_group2;
     by chr region_num;
proc means data=br_regions_group2 noprint;
     by chr region_num;
     var superregion_start superregion_stop 
flag_DMR_1p4_1 flag_DAR_1p4_1 flag_DMR_1p4_72 flag_DAR_1p4_72
flag_DMR_10_1 flag_DAR_10_1 flag_DAR_10_1 flag_DMR_10_72 flag_DAR_10_72
;
     output out=br_regions_collapsed(drop=_TYPE_ _FREQ_) min(superregion_start)=superregion_start
                                     max(superregion_stop)=superregion_stop
                                     sum(flag_DMR_1p4_1)=N_DMR_1p4_1 max(flag_DMR_1p4_1)=flag_DMR_1p4_1
                                     sum(flag_DAR_1p4_1)=N_DAR_1p4_1 max(flag_DAR_1p4_1)=flag_DAR_1p4_1
                                     sum(flag_DMR_1p4_72)=N_DMR_1p4_72 max(flag_DMR_1p4_72)=flag_DMR_1p4_72
                                     sum(flag_DAR_1p4_72)=N_DAR_1p4_72 max(flag_DAR_1p4_72)=flag_DAR_1p4_72
                                     sum(flag_DMR_10_1)=N_DMR_10_1 max(flag_DMR_10_1)=flag_DMR_10_1
                                     sum(flag_DAR_10_1)=N_DAR_10_1 max(flag_DAR_10_1)=flag_DAR_10_1
                                     sum(flag_DMR_10_72)=N_DMR_10_72 max(flag_DMR_10_72)=flag_DMR_10_72
                                     sum(flag_DAR_10_72)=N_DAR_10_72 max(flag_DAR_10_72)=flag_DAR_10_72;
run;


data br_Regions_collapsed2;
   retain chr superregion_start superregion_stop sum_: flag_: ;
   set br_regions_collapsed;
   rename superregion_start=region_start superregion_stop=region_stop;
run;

proc freq data=br_Regions_collapsed2 noprint;
  tables flag_DAR_1p4_1*flag_DMR_1p4_1*
         flag_DAR_1p4_72*flag_DMR_1p4_72*
         flag_DAR_10_1*flag_DMR_10_1*
         flag_DAR_10_72*flag_DMR_10_72 / out=br_Region_cnt;
run;

proc print data=br_region_cnt;
run;

/*

  flag_      flag_      flag_      flag_                           flag_     flag_
DAR_1p4_   DMR_1p4_   DAR_1p4_   DMR_1p4_     flag_      flag_    DAR_10_   DMR_10_
    1          1         72         72      DAR_10_1   DMR_10_1      72        72     COUNT

    0          0          0          0          0          0         0         1        872
    0          0          0          0          0          0         1         0      19659
    0          0          0          0          0          0         1         1         24
    0          0          0          0          0          1         0         0        106
    0          0          0          0          0          1         1         0          2
    0          0          0          0          1          0         0         0      33408
    0          0          0          0          1          0         0         1         37
    0          0          0          0          1          0         1         0       1704
    0          0          0          0          1          0         1         1          6
    0          0          0          0          1          1         0         0         11
    0          0          0          0          1          1         1         0          3
    0          0          0          1          0          0         0         0        169
    0          0          0          1          0          0         0         1         13
    0          0          0          1          0          0         1         0          2
    0          0          0          1          1          0         0         0         19
    0          0          0          1          1          0         0         1          2
    0          0          0          1          1          0         1         0          1
    0          0          1          0          0          0         0         0      31358
    0          0          1          0          0          0         0         1         28
    0          0          1          0          0          0         1         0       6226
    0          0          1          0          0          0         1         1          9
    0          0          1          0          1          0         0         0       1805
    0          0          1          0          1          0         0         1          6
    0          0          1          0          1          0         1         0        739
    0          0          1          0          1          0         1         1          4
    0          0          1          0          1          1         0         0          1
    0          0          1          0          1          1         1         0          1
    0          0          1          1          0          0         0         0          7
    0          0          1          1          0          0         1         0          2
    0          0          1          1          1          0         0         0          1
    0          1          0          0          0          0         0         0        113
    0          1          0          0          0          0         0         1          2
    0          1          0          0          0          0         1         0          2
    0          1          0          0          0          0         1         1          1
    0          1          0          0          0          1         0         0          9
    0          1          0          0          1          0         0         0          6
    0          1          0          0          1          1         0         0          1
    0          1          0          0          1          1         1         0          1
    0          1          0          1          0          0         0         0          1
    0          1          1          0          1          0         1         0          1
    0          1          1          0          1          1         0         0          1
    0          1          1          0          1          1         1         0          1
    1          0          0          0          0          0         0         0      35070
    1          0          0          0          0          0         0         1         31
    1          0          0          0          0          0         1         0       1197
    1          0          0          0          0          0         1         1          2
    1          0          0          0          0          1         0         0          2
    1          0          0          0          1          0         0         0       6045
    1          0          0          0          1          0         0         1          9
    1          0          0          0          1          0         1         0        520
    1          0          0          0          1          0         1         1          3
    1          0          0          0          1          1         0         0          3
    1          0          0          0          1          1         1         1          1
    1          0          0          1          0          0         0         0          4
    1          0          0          1          0          0         0         1          2
    1          0          0          1          0          0         1         0          1
    1          0          0          1          1          0         0         0          4
    1          0          0          1          1          0         1         0          1
    1          0          1          0          0          0         0         0       1643
    1          0          1          0          0          0         1         0        617
    1          0          1          0          0          0         1         1          2
    1          0          1          0          1          0         0         0        579
    1          0          1          0          1          0         1         0        294
    1          0          1          0          1          0         1         1          1
    1          0          1          0          1          1         1         0          1
    1          0          1          1          0          0         0         0          1
    1          0          1          1          1          0         1         0          1
    1          1          0          0          0          0         0         0          7
    1          1          0          0          1          0         0         0          4
    1          1          1          0          1          0         1         0          2

*/


proc export data=at_Regions_collapsed2
    outfile="!PATCON/DTRA/at_rad_merged_regions.csv"
    dbms=csv replace;
run;

proc export data=br_Regions_collapsed2
    outfile="!PATCON/DTRA/br_rad_merged_regions.csv"
    dbms=csv replace;
run;


/* For AT: significant in 3 tests
   For BR: significant in 4 tests */

data at_regions_keep;
   set at_regions_collapsed2;
   if sum (of flag_:) >= 3 then output;
   keep chr region_start region_stop;
run;


data br_regions_keep;
   set br_regions_collapsed2;
   if sum (of flag_:) >= 4 then output;
   keep chr region_start region_stop;
run;

proc export data=at_regions_keep
      outfile="!PATCON/DTRA/arabidopsis_merged_regions_3sig.bed"
      dbms=tab replace;
      putnames=no;
run;


proc export data=br_regions_keep
      outfile="!PATCON/DTRA/brassica_merged_regions_4sig.bed"
      dbms=tab replace;
      putnames=no;
run;



/* Output BED files for MEME analysis */


data at_Regions_collapsed3;
   retain chr superregion_start superregion_stop  ;
   set at_regions_collapsed;
   num_sig=sum(of flag_: );
   keep chr  superregion_start  superregion_stop region_num num_sig;
run;

data at_region2num;
  set at_regions_group2;
  keep chr region_num site_type region_start region_stop;
run;

proc sort data=at_regions_collapsed3 nodup;
   by chr region_num;
proc sort data=at_region2num;
   by chr region_num;
run;

data at_region2num2;
  merge at_regions_collapsed3 (in=in1)  at_region2num (in=in2);
  by chr region_num;
  if in1 and in2;
run;

data at_region2num_id;
   set at_region2num2;
   length merged_ID $32.;
   length region_ID $32.;
   merged_ID = catx("_",chr,region_num,superregion_start,superregion_stop);
   region_ID = catx("_",chr,region_num,region_start,region_stop);
   superregion_start2=superregion_start-6;
   superregion_stop2=superregion_stop+6;
   region_start2=region_start-6;
   region_stop2=region_stop+6;
   if num_sig >=3;
run;

data at_super_meme;
  retain chr superregion_start2 superregion_stop2 merged_ID;
  set at_region2num_id;
  keep chr superregion_start2 superregion_stop2 merged_ID;
run;

data at_region_meme;
  retain chr region_start2 region_stop2 region_ID;
  set at_region2num_id;
  keep  chr region_start2 region_stop2 region_ID;
run;

proc sort data=at_super_meme nodup;
   by chr superregion_start2 superregion_stop2;
run;

proc sort data=at_region_meme nodup;
   by chr region_start2 region_stop2;
run;







data br_Regions_collapsed3;
   retain chr superregion_start superregion_stop  ;
   set br_regions_collapsed;
   num_sig=sum(of flag_: );
   keep chr  superregion_start  superregion_stop region_num num_sig;
run;

data br_region2num;
  set br_regions_group2;
  keep chr region_num site_type region_start region_stop;
run;

proc sort data=br_regions_collapsed3 nodup;
   by chr region_num;
proc sort data=br_region2num;
   by chr region_num;
run;

data br_region2num2;
  merge br_regions_collapsed3 (in=in1)  br_region2num (in=in2);
  by chr region_num;
  if in1 and in2;
run;

data br_region2num_id;
   set br_region2num2;
   length merged_ID $32.;
   length region_ID $32.;
   merged_ID = catx("_",chr,region_num,superregion_start,superregion_stop);
   region_ID = catx("_",chr,region_num,region_start,region_stop);
   superregion_start2=superregion_start-6;
   superregion_stop2=superregion_stop+6;
   region_start2=region_start-6;
   region_stop2=region_stop+6;
   if num_sig >=4;
run;

data br_super_meme;
  retain chr superregion_start2 superregion_stop2 merged_ID;
  set br_region2num_id;
  keep chr superregion_start2 superregion_stop2 merged_ID;
run;

data br_region_meme;
  retain chr region_start2 region_stop2 region_ID;
  set br_region2num_id;
  keep  chr region_start2 region_stop2 region_ID;
run;

proc sort data=br_super_meme nodup;
   by chr superregion_start2 superregion_stop2;
run;

proc sort data=br_region_meme nodup;
   by chr region_start2 region_stop2;
run;

proc export data=at_super_meme
    outfile="!PATCON/DTRA/At_sig3_merged_regions_for_meme.bed"
    dbms=tab replace; putnames=no;
run;

proc export data=at_region_meme
    outfile="!PATCON/DTRA/At_sig3_regions_for_meme.bed"
    dbms=tab replace; putnames=no;
run;
proc export data=br_super_meme
    outfile="!PATCON/DTRA/Br_sig4_merged_regions_for_meme.bed"
    dbms=tab replace; putnames=no;
run;

proc export data=br_region_meme
    outfile="!PATCON/DTRA/Br_sig4_regions_for_meme.bed"
    dbms=tab replace; putnames=no;
run;




   
