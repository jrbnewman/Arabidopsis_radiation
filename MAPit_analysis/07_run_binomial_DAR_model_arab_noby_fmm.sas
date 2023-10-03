

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


data cyto2meth_01gy;
   set wgbslocA.cytosine2acc_: ;
   where comparison="01Gy_0Gy";
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
   set wgbslocA.cytosine2acc_: ;
   where comparison="1Gy_0Gy";
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
  by comparison site_type chr region_num group units start_pos stop_pos rep;
run;


proc freq data=data_for_models noprint;
 tables group*units / out=check;
run;


ods listing close;

proc glimmix data=data_for_models;
   by comparison site_type chr region_num;
   class group units start_pos;
   model methyl_c / total_c = group|units start_pos / dist=binomial link=logit ddfm=kr chisq   oddsratio(diff=all);
   random group*units*start_pos / solution ;

   /* order                     CTL_0U CTL_100U TRT_0U TRT_100U */
   contrast "TRT_v_CTL" group*units        1       -1    -1        1  ;
   contrast "TRT_100U_0U" group*units      0        0    -1        1  ;
   contrast "CTL_100U_0U" group*units     -1        1     0        0  ;

   estimate "TRT_v_CTL" group*units        1       -1    -1        1 ;
   estimate "TRT_100U_0U" group*units      0        0    -1        1  ;
   estimate "CTL_100U_0U" group*units     -1        1     0        0  ;

 output out=gmxout pred(ilink)=gpredy lcl(ilink)=lower ucl(ilink)=upper
 pred=pred Resid=Resid student=stu;
 ods output tests3=anova
 contrasts=con estimates=est FitStatistics=fits  OddsRatios=odds;
 run;
 quit;


/* Make permanent */

data wgbslocA.bin_dar_pred2_all;   set gmxout; run;
data wgbslocA.bin_dar_anova2_all;   set anova; run;
data wgbslocA.bin_dar_contrasts2_all;   set con; run;
data wgbslocA.bin_dar_estimates2_all;   set est; run;
data wgbslocA.bin_dar_fitstats2_all;   set fits; run;
data wgbslocA.bin_dar_odds2_all;   set odds; run;



