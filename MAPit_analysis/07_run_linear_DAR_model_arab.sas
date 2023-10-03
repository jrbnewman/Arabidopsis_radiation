
/* Run binomial model on DMR or DAR */

*libname wgbslocA "/blue/concannon/share/jnewman/brassica_wgbs/sas_data2";
*libname wgbslocB "/blue/concannon/share/jnewman/brassica_wgbs/sas_data";
libname wgbslocA "/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs";


ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


data cyto2meth_01gy;
   set wgbslocA.cytosine_to_acc_region_ge2_v4 ;
   where comparison="01Gy_0Gy" ;
   keep region_num site_type chr start_pos stop_pos comparison;
run;

data counts_01gy;
   set wgbslocA.meth_data_gc_: ;
   length group $3.;
   where treatment="01Gy" or treatment="0Gy";
   if treatment="0Gy" then group="CTL";
   else group="TRT";
run;

proc sort data=cyto2meth_01gy;
   by site_type chr start_pos stop_pos;
proc sort data=counts_01gy;
   by site_type chr start_pos stop_pos;
run;

data data_01gy;
  merge cyto2meth_01gy (in=in1) counts_01gy (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;




data cyto2meth_1gy;
   set wgbslocA.cytosine_to_acc_region_ge2_v4 ;
   where comparison="1Gy_0Gy" ;
   keep region_num site_type chr start_pos stop_pos comparison;
run;

data counts_1gy;
   set wgbslocA.meth_data_gc_: ;
   length group $3.;
   where treatment="1Gy" or treatment="0Gy";
   if treatment="0Gy" then group="CTL";
   else group="TRT";
run;

proc sort data=cyto2meth_1gy;
   by site_type chr start_pos stop_pos;
proc sort data=counts_1gy;
   by site_type chr start_pos stop_pos;
run;

data data_1gy;
  merge cyto2meth_1gy (in=in1) counts_1gy (in=in2);
  by site_type chr start_pos stop_pos;
  if in1 and in2;
run;


data data_for_models;
   set data_01gy data_1gy;
run;

proc sort data=data_for_models;
  by comparison site_type chr region_num  start_pos stop_pos  group units rep;
run;

proc means data=data_For_models noprint;
   by comparison site_type chr region_num  start_pos stop_pos group units;
   var perc_methyl;
   output out=data_for_models2 mean=perc_methyl;
run;

proc transpose data=data_for_models2 out=data_For_models_sbys;
   by comparison site_type chr region_num  start_pos stop_pos group;
    id units;
     var perc_methyl;
run;


data data_for_models3;
  set data_for_models_sbys;
  perc_access = _100U - _0U ;
run;



proc sort data=data_for_models2;
  by comparison site_type chr region_num group units start_pos stop_pos;
proc sort data=data_for_models3;
  by comparison site_type chr region_num group  start_pos stop_pos;
run;

ods listing close;
proc mixed data=data_for_models2 ;
  by comparison site_type chr region_num;
  class group units;
  model perc_methyl = group|units / htype=3;
 *output out=gmxout pred=pred Resid=Resid student=stu;
 ods output tests3=anova1 FitStatistics=fits;
 run;
 quit;


proc mixed data=data_for_models3 ;
  by comparison site_type chr region_num;
  class group;
  model perc_access = group / htype=3;
 *output out=gmxout pred=pred Resid=Resid student=stu;
 ods output tests3=anova2 FitStatistics=fits2;
 run;
 quit;


data wgbslocA.mixedmeth_dar_anova_all;   set anova; run;
data wgbslocA.mixedmeth_dar_fitstats_all;   set fits; run;



data wgbslocA.mixedacc_dar_anova_all;   set anova2; run;
data wgbslocA.mixedacc_dar_fitstats_all;   set fits2; run;




