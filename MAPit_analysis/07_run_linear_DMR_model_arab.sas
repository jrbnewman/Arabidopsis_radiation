

/* Run binomial model on DMR or DAR */

libname wgbslocA "/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs";
*libname wgbslocB "/blue/concannon/share/jnewman/brassica_wgbs/sas_data";

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

%let chrom=%getsys(chrom);

data cyto2meth_01gy;
   set wgbslocA.cytosine_to_meth_region_ge2_v4 ;
   *set wgbslocA.cytosine2acc_&chrom. ;
   where comparison="0Gy_vs_01G" and chr="&chrom.";
   keep region_num site_type chr start_pos stop_pos comparison;
run;

data counts_01gy;
   set wgbslocA.meth_data_cg_chg_chh_&chrom. ;
   *set wgbslocA.meth_data_gc_&chrom. ;
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
   set wgbslocA.cytosine_to_meth_region_ge2_v4 ;
   *set wgbslocA.cytosine2acc_&chrom. ;
   where comparison="0Gy_vs_1Gy" and chr="&chrom.";
   keep region_num site_type chr start_pos stop_pos comparison;
run;

data counts_1gy;
   set wgbslocA.meth_data_cg_chg_chh_&chrom. ;
   *set wgbslocA.meth_data_gc_&chrom. ;
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
  by comparison site_type chr region_num group units start_pos stop_pos rep;
run;

proc means data=data_For_models noprint;
  by comparison site_type chr region_num group units start_pos stop_pos;
  var perc_methyl;
  output out=data_for_models2 mean=perc_methyl;
run;
ods listing close;
proc mixed data=data_for_models2 ;
  by comparison site_type chr region_num;
  class group;
  model perc_methyl = group / htype=3;
 *output out=gmxout pred=pred Resid=Resid student=stu;
 ods output tests3=anova FitStatistics=fits;
 run;
 quit;


data wgbslocA.mixed_dmr_anova_&chrom.;   set anova; run;
data wgbslocA.mixed_dmr_fitstats_&chrom.;   set fits; run;

