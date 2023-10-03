libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Create design file */

data design_file;
   format treatment $10.;
   format rep best12.;
   format units $10.;
   input treatment $ rep units $;
   datalines;
   0Gy 1 0U
   0Gy 2 0U 
   01Gy 1 0U
   01Gy 2 0U
   1Gy 1 0U
   1Gy 2 0U
   0Gy 1 100U
   0Gy 2 100U 
   01Gy 1 100U
   01Gy 2 100U
   1Gy 1 100U
   1Gy 2 100U 
   ;
run;


data wgbsA.design_file;
  set design_file;
run;




