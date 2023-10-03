ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

/* Import design file for Arabidopsis radiation sensitivity experiment */

proc import datafile="!PATCON/arabidopsis/design_files/arabidopsis_rs_design_file.csv"
    out=design dbms=csv replace; guessingrows=300;
run;

/* Check: are there the same number of treat*time*rep ?*/

proc freq data=design noprint;
  tables treatment*time*replicate / out=sample_count;
run;

proc print data=sample_count;
run;

/*
  Obs    treatment            time       replicate    COUNT    PERCENT

    1      0.1gy                 1               1      5      2.77778
    2      0.1gy                 1               2      5      2.77778
    3      0.1gy                 1               3      5      2.77778
    4      0.1gy                 3               1      5      2.77778
    5      0.1gy                 3               2      5      2.77778
    6      0.1gy                 3               3      5      2.77778
    7      0.1gy                24               1      5      2.77778
    8      0.1gy                24               2      5      2.77778
    9      0.1gy                24               3      5      2.77778
   10      0.1gy                72               1      5      2.77778
   11      0.1gy                72               2      5      2.77778
   12      0.1gy                72               3      5      2.77778
   13      1gy                   1               1      5      2.77778
   14      1gy                   1               2      5      2.77778
   15      1gy                   1               3      5      2.77778
   16      1gy                   3               1      5      2.77778
   17      1gy                   3               2      5      2.77778
   18      1gy                   3               3      5      2.77778
   19      1gy                  24               1      5      2.77778
   20      1gy                  24               2      5      2.77778
   21      1gy                  24               3      5      2.77778
   22      1gy                  72               1      5      2.77778
   23      1gy                  72               2      5      2.77778
   24      1gy                  72               3      5      2.77778
   25      Mock                  1               1      5      2.77778
   26      Mock                  1               2      5      2.77778
   27      Mock                  1               3      5      2.77778
   28      Mock                  3               1      5      2.77778
   29      Mock                  3               2      5      2.77778
   30      Mock                  3               3      5      2.77778
   31      Mock                 24               1      5      2.77778
   32      Mock                 24               2      5      2.77778
   33      Mock                 24               3      5      2.77778
   34      Mock                 72               1      5      2.77778
   35      Mock                 72               2      5      2.77778
   36      Mock                 72               3      5      2.77778

Looks okay!
*/

/* Make permenant */

data rs.arab_design_file;
  set design;
run;

