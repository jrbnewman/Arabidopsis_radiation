

/* Run binomial model on DMR or DAR */

*libname wgbslocA "/blue/concannon/share/jnewman/brassica_wgbs/sas_data2";
*libname wgbslocB "/blue/concannon/share/jnewman/brassica_wgbs/sas_data";
libname wgbslocA "/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs";

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

*%let chrom=%getsys(chrom);

*%let chrom=5;

data cyto2meth_01gy;
   set wgbslocA.cytosine2meth_: ;
   *set wgbslocA.cytosine2acc_&chrom. ;
   where comparison="0Gy_vs_01G";
   keep region_num site_type chr start_pos stop_pos comparison;
run;

data counts_01gy;
   set wgbslocA.meth_data_cg_chg_chh_: ;
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
   set wgbslocA.cytosine2meth_: ;
   *set wgbslocA.cytosine2acc_&chrom. ;
   where comparison="0Gy_vs_1Gy";
   keep region_num site_type chr start_pos stop_pos comparison;
run;

data counts_1gy;
   set wgbslocA.meth_data_cg_chg_chh_: ;
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

ods listing close;

proc fmm data=data_for_models;
   by comparison site_type chr region_num;
   class group start_pos;
   model methyl_c / total_c = group start_pos / dist=bb ;
   *random group*start_pos / solution ;

   /* order                     CTL TRT */
   *contrast "TRT_v_CTL" group   -1  1    ;
   *estimate "TRT_v_CTL" group   -1  1    ;

 output out=gmxout pred=gpredy pred=pred Resid=Resid student=stu;
 ods output tests3=anova  FitStatistics=fits;
 run;
 quit;


/* Make permanent */

data wgbslocA.betabin_dmr_pred;   set gmxout; run;
data wgbslocA.betabin_dmr_anova;   set anova; run;
data wgbslocA.betabin_dmr_fitstats;   set fits; run;




