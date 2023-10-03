

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

ods listing close;

proc glimmix data=data_for_models;
   by comparison site_type chr region_num;
   class group start_pos;
   model methyl_c / total_c = group  / dist=binomial link=logit ddfm=bw;
   random intercept / subject=start_pos ;

   /* order                     CTL TRT */
   contrast "TRT_v_CTL" group   -1  1    ;
   estimate "TRT_v_CTL" group   -1  1    ;

 output out=gmxout pred(ilink)=gpredy lcl(ilink)=lower ucl(ilink)=upper
 pred=pred Resid=Resid student=stu;
 ods output tests3=anova
 contrasts=con estimates=est FitStatistics=fits;
 run;
 quit;


/* Make permanent */

data wgbslocA.bin2_dmr_pred_&chrom.;   set gmxout; run;
data wgbslocA.bin2_dmr_anova_&chrom.;   set anova; run;
data wgbslocA.bin2_dmr_contrasts_&chrom.;   set con; run;
data wgbslocA.bin2_dmr_estimates_&chrom.;   set est; run;
data wgbslocA.bin2_dmr_fitstats_&chrom.;   set fits; run;



