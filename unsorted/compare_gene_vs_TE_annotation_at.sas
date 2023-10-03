/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;




data cg_up_dmr_01_72 cg_up_dmr_1_72 cg_dn_dmr_01_72 cg_dn_dmr_1_72
     chg_up_dmr_01_72 chg_up_dmr_1_72 chg_dn_dmr_01_72 chg_dn_dmr_1_72
     chh_up_dmr_01_72 chh_up_dmr_1_72 chh_dn_dmr_01_72 chh_dn_dmr_1_72;
     set arabMAP.results_by_dmr_annot;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     /* FET */
     if site_type="CG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72;
     end;


     /* FET */
     if site_type="CHG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72;
     end;

     /* FET */
     if site_type="CHH" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72;
     end;


     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
run;

proc sort data=cg_up_dmr_01_72 nodup; by _all_; run;
proc sort data=cg_dn_dmr_01_72 nodup; by _all_; run;
proc sort data=cg_up_dmr_1_72 nodup; by _all_; run;
proc sort data=cg_dn_dmr_1_72 nodup; by _all_; run;

proc sort data=chg_up_dmr_01_72 nodup; by _all_; run;
proc sort data=chg_dn_dmr_01_72 nodup; by _all_; run;
proc sort data=chg_up_dmr_1_72 nodup; by _all_; run;
proc sort data=chg_dn_dmr_1_72 nodup; by _all_; run;

proc sort data=chh_up_dmr_01_72 nodup; by _all_; run;
proc sort data=chh_dn_dmr_01_72 nodup; by _all_; run;
proc sort data=chh_up_dmr_1_72 nodup; by _all_; run;
proc sort data=chh_dn_dmr_1_72 nodup; by _all_; run;




data up_dar_01_72 up_dar_1_72 dn_dar_01_72 dn_dar_1_72;
     set arabMAP.results_by_dar_annot;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72;
     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
run;



proc sort data=up_dar_01_72 nodup; by _all_; run;
proc sort data=dn_dar_01_72 nodup; by _all_; run;
proc sort data=up_dar_1_72 nodup; by _all_; run;
proc sort data=dn_dar_1_72 nodup; by _all_; run;




data cg_up_dmr_01_72_te cg_up_dmr_1_72_te cg_dn_dmr_01_72_te cg_dn_dmr_1_72_te
     chg_up_dmr_01_72_te chg_up_dmr_1_72_te chg_dn_dmr_01_72_te chg_dn_dmr_1_72_te
     chh_up_dmr_01_72_te chh_up_dmr_1_72_te chh_dn_dmr_01_72_te chh_dn_dmr_1_72_te;
     set arabMAP.results_by_dmr_annot_TE;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     /* FET */
     if site_type="CG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72_te;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72_te;
     end;


     /* FET */
     if site_type="CHG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72_te;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72_te;
     end;

     /* FET */
     if site_type="CHH" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_te;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_te;
     end;


     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
     rename feature=feature_TE distance_to_tss=distance_to_tss_TE geneID=geneID_TE;
run;

proc sort data=cg_up_dmr_01_72_te nodup; by _all_; run;
proc sort data=cg_dn_dmr_01_72_te nodup; by _all_; run;
proc sort data=cg_up_dmr_1_72_te nodup; by _all_; run;
proc sort data=cg_dn_dmr_1_72_te nodup; by _all_; run;

proc sort data=chg_up_dmr_01_72_te nodup; by _all_; run;
proc sort data=chg_dn_dmr_01_72_te nodup; by _all_; run;
proc sort data=chg_up_dmr_1_72_te nodup; by _all_; run;
proc sort data=chg_dn_dmr_1_72_te nodup; by _all_; run;

proc sort data=chh_up_dmr_01_72_te nodup; by _all_; run;
proc sort data=chh_dn_dmr_01_72_te nodup; by _all_; run;
proc sort data=chh_up_dmr_1_72_te nodup; by _all_; run;
proc sort data=chh_dn_dmr_1_72_te nodup; by _all_; run;




data up_dar_01_72_te up_dar_1_72_te dn_dar_01_72_te dn_dar_1_72_te;
     set arabMAP.results_by_dar_annot_te;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_72_te;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_72_te;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72_te;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72_te;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72_te;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72_te;
     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
     rename feature=feature_TE distance_to_tss=distance_to_tss_TE geneID=geneID_TE;
run;




proc sort data=up_dar_01_72_te nodup; by _all_; run;
proc sort data=dn_dar_01_72_te nodup; by _all_; run;
proc sort data=up_dar_1_72_te nodup; by _all_; run;
proc sort data=dn_dar_1_72_te nodup; by _all_; run;


%macro mergeData(inData);
proc sort data=&inData.;
  by comparison site_type chr region_start region_stop;
proc sort data=&inData._TE;
  by comparison site_type chr region_start region_stop;
run;


data &inData._gn_vs_te;
  merge &inData. (in=in1) &inData._TE (in=in2);
  by comparison site_type chr region_start region_stop;
  if in1 and in2;
run;

data flag_&inData._gn_vs_te;
  set &inData._gn_vs_te;
  if feature ne "Intergenic" then flag_has_gene=1; else flag_has_gene=0;
  if feature="exon" or feature="intron" then flag_gene_body=1; else flag_gene_body=0;

  if feature_TE ne "Intergenic" then flag_has_TE=1; else flag_has_TE=0;
  if feature_TE="exon" or feature="intron" then flag_TE_body=1; else flag_TE_body=0;
run;

proc freq data=flag_&inData._gn_vs_te noprint;
  table flag_has_gene*flag_has_TE / out=ctab1;
  table flag_has_gene*flag_gene_body*flag_has_TE*flag_TE_body / out=ctab2;
  table flag_has_gene*flag_has_TE*feature*feature_TE / out=ctab3;
run;

proc print data=ctab1;
proc print data=ctab2;
proc print data=ctab3;
run;

%mend;
%mergeData(cg_up_dmr_01_72);
%mergeData(cg_up_dmr_1_72);
%mergeData(cg_dn_dmr_01_72);
%mergeData(cg_dn_dmr_1_72);

%mergeData(chg_up_dmr_01_72);
%mergeData(chg_up_dmr_1_72);
%mergeData(chg_dn_dmr_01_72);
%mergeData(chg_dn_dmr_1_72);

%mergeData(chh_up_dmr_01_72);
%mergeData(chh_up_dmr_1_72);
%mergeData(chh_dn_dmr_01_72);
%mergeData(chh_dn_dmr_1_72);

%mergeData(up_dar_01_72);
%mergeData(up_dar_1_72);
%mergeData(dn_dar_01_72);
%mergeData(dn_dar_1_72);

/* Get regions and methylation levels for DMR/DAR plots:

CHH hyper 10cGy, genic only
CHH hypo 10cGy, genic only
CHH hyper 10cGy, TE only
CHH hypo 10cGy, TE only
CHH hyper 100cGy, genic only
CHH hypo 100cGy, genic only
CHH hyper 100cGy, TE only
CHH hypo 100cGy, TE only

DAR hyper 10cGy, genic only
DAR hypo 10cGy, genic only
DAR hyper 10cGy, TE only
DAR hypo 10cGy, TE only
DAR hyper 100cGy, genic only
DAR hypo 100cGy, genic only
DAR hyper 100cGy, TE only
DAR hypo 100cGy, TE only */

* GC methylation data;
data meth_data;
  set arabMAP.methylation_data_cg_chg_chh;
  where site_type="CHH";
  keep chr stop_pos treatment rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by chr stop_pos treatment  rep ;
proc means data=meth_data noprint;
  by chr stop_pos treatment   ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;

data meth_data3;
  set meth_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=meth_data3 out=meth_sbys10;
  where total_C >= 10;
  by chr stop_pos;
  id treatment ;
  var perc_methyl2;
run;

data meth_sbys10_2;
  set meth_sbys10;
  if _01Gy ne . and _0Gy ne . then _01Gy_common=_01Gy; else _01Gy_common=.;
  if _1Gy ne . and _0Gy ne . then _1Gy_common=_1Gy; else _1Gy_common=.;
  if (_1Gy ne . or _01Gy ne .) and _0Gy ne . then _0Gy_common=_0Gy; else _0Gy_common=.;

  if _01Gy_common ne . and _1Gy_common ne .  then do;
    _01Gy_common_all = _01Gy_common;
    _1Gy_common_all = _1Gy_common;
    _0Gy_common_all = _0Gy_common;
    end;
  else do;
   _01Gy_common_all = .;
   _1Gy_common_all = .;
   _0Gy_common_all = . ;
    end;
  rename stop_pos=pos;
run;


%macro getMethyl(inRegion,methVar);

data genic te;
   set flag_&inRegion._gn_vs_te;
   dar_center=int((region_start + region_stop) / 2);
   plot_start=dar_center - 999;
   plot_stop=dar_center + 999;
   if feature="Intergenic" and feature_te ne "Intergenic" then output te;
   else if feature ne "Intergenic" and feature_te = "Intergenic" then output genic;
   keep chr dar_center plot_start plot_stop;
run;


data genic_2;
  set genic;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

data te_2;
  set te;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=genic_2;
   by chr pos;
proc sort data=te_2;
   by chr pos;
proc sort data=meth_sbys10_2;
   by chr pos;
run;

data genic_w_meth;
  merge genic_2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;
data te_w_meth;
  merge te_2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data genic_w_meth2;
  set genic_w_meth;
  distance_to_center=dar_Center-pos;
run;


data genic_w_meth3;
  set genic_w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=genic_w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=genic_w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_&inRegion._genic
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_&inRegion._genic ;
  by grouped_pos2;
run;


proc means data=mean_&inRegion._genic noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_&inRegion._genic_2
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_&inRegion._genic_1;
  set mean_&inRegion._genic;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_&inRegion._genic_3;
  set mean_&inRegion._genic_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;




data te_w_meth2;
  set te_w_meth;
  distance_to_center=dar_Center-pos;
run;


data te_w_meth3;
  set te_w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=te_w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=te_w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_&inRegion._te
  mean(_01Gy_common)=mean_01Gy_common_te
  mean(_1Gy_common)=mean_1Gy_common_te
  mean( _0Gy_common)=mean_0Gy_common_te
  mean(_01Gy_common_all)=mean_01Gy_common_all_te
  mean(_1Gy_common_all)=mean_1Gy_common_all_te
  mean( _0Gy_common_all)=mean_0Gy_common_all_te;
run;

proc sort data=mean_&inRegion._te ;
  by grouped_pos2;
run;


proc means data=mean_&inRegion._te noprint;
  by grouped_pos2  ;
  var  mean_01Gy_common_te mean_1Gy_common_te  mean_0Gy_common_te
       mean_01Gy_common_all_te mean_1Gy_common_all_te  mean_0Gy_common_all_te ;
  output out=mean_&inRegion._te_2
  mean(mean_01Gy_common_te)=mean_01Gy_common_te
  mean(mean_1Gy_common_te)=mean_1Gy_common_te
  mean(mean_0Gy_common_te)=mean_0Gy_common_te
  mean(mean_01Gy_common_all_te)=mean_01Gy_common_all_te
  mean(mean_1Gy_common_all_te)=mean_1Gy_common_all_te
  mean(mean_0Gy_common_all_te)=mean_0Gy_common_all_te;
run;


data mean_&inRegion._te_1;
  set mean_&inRegion._te;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_&inRegion._te_3;
  set mean_&inRegion._te_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;



proc sort data=mean_&inRegion._genic_1;
   by pos;
proc sort data=mean_&inRegion._te_1;
   by pos;
proc sort data=mean_&inRegion._genic_3;
   by pos;
proc sort data=mean_&inRegion._te_3;
   by pos;
run;


data dmr_&inRegion._100;
  merge mean_&inRegion._genic_3 mean_&inRegion._te_3;
   by pos;
  keep pos mean_0Gy_common &methVar. mean_0Gy_common_te &methVar._te ;
run;

data dmr_&inRegion._100_2;
retain   pos mean_0Gy_common &methVar. mean_0Gy_common_te &methVar._te ;
  set dmr_&inRegion._100;
run;



proc export data=dmr_&inRegion._100_2
     outfile="!HOME/concannon/DTRA/dmr_plot_&inRegion._100_bins.csv"
     dbms=csv replace;
run;


%mend;


%getMethyl(chh_up_dmr_01_72,mean_01Gy_common);
%getMethyl(chh_up_dmr_1_72,mean_1Gy_common);
%getMethyl(chh_dn_dmr_01_72,mean_01Gy_common);
%getMethyl(chh_dn_dmr_1_72,mean_1Gy_common);










* GC methylation data;
data gc_data;
  set arabMAP.methylation_data_gc;
  where flag_normalized=1;
  keep chr stop_pos treatment units rep total_C_norm methyl_C_norm perc_methyl_norm;
run;

proc sort data=gc_data;
  by chr stop_pos treatment units rep ;
proc means data=gc_data noprint;
  by chr stop_pos treatment units  ;
  var total_C_norm methyl_C_norm perc_methyl_norm;
  output out=gc_data2 sum(total_C_norm)=total_C sum(methyl_C_norm)=methyl_C mean(perc_methyl_norm)=perc_methyl;
run;

data gc_data3;
  set gc_data2;
  perc_methyl2=(methyl_C / total_C) * 100 ;
run;


proc transpose data=gc_data3 out=gc_sbys10;
  where total_C >= 10;
  by chr stop_pos;
  id treatment units;
  var perc_methyl2;
run;

data gc_sbys10_2;
  set gc_sbys10;
  if _01Gy0U ne  . and _01Gy100U ne . then _01Gy_100U_0U=_01Gy100U - _01Gy0U;
  else _01Gy_100U_0U=.;

  if _1Gy0U ne . and _1Gy100U ne . then _1Gy_100U_0U=_1Gy100U - _1Gy0U;
  else _1Gy_100U_0U=.;

  if _0Gy0U ne . and _0Gy100U ne . then _0Gy_100U_0U=_0Gy100U - _0Gy0U;
  else _0Gy_100U_0U=.;

  if _01Gy_100U_0U ne . and _0Gy_100U_0U ne . then _01Gy_100U_0U_common=_01Gy_100U_0U; else _01Gy_100U_0U_common=.;
  if _1Gy_100U_0U ne . and _0Gy_100U_0U ne . then _1Gy_100U_0U_common=_1Gy_100U_0U; else _1Gy_100U_0U_common=.;
  if (_1Gy_100U_0U ne . or _01Gy_100U_0U ne .) and _0Gy_100U_0U ne . then _0Gy_100U_0U_common=_0Gy_100U_0U; else _0Gy_100U_0U_common=.;

  if _01Gy_100U_0U_common ne . and _1Gy_100U_0U_common ne .  then do;
    _01Gy_100U_0U_common_all = _01Gy_100U_0U_common;
    _1Gy_100U_0U_common_all = _1Gy_100U_0U_common;
    _0Gy_100U_0U_common_all = _0Gy_100U_0U_common;
    end;
  else do;
   _01Gy_100U_0U_common_all = .;
   _1Gy_100U_0U_common_all = .;
   _0Gy_100U_0U_common_all = . ;
    end;
  rename stop_pos=pos;
run;


%macro getMethyl(inRegion,methVar);

data genic te;
   set flag_&inRegion._gn_vs_te;
   dar_center=int((region_start + region_stop) / 2);
   plot_start=dar_center - 999;
   plot_stop=dar_center + 999;
   if feature="Intergenic" and feature_te ne "Intergenic" then output te;
   else if feature ne "Intergenic" and feature_te = "Intergenic" then output genic;
   keep chr dar_center plot_start plot_stop;
run;


data genic_2;
  set genic;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

data te_2;
  set te;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=genic_2;
   by chr pos;
proc sort data=te_2;
   by chr pos;
proc sort data=gc_sbys10_2;
   by chr pos;
run;

data genic_w_meth;
  merge genic_2 (in=in1) gc_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;
data te_w_meth;
  merge te_2 (in=in1) gc_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data genic_w_meth2;
  set genic_w_meth;
  distance_to_center=dar_Center-pos;
run;


data genic_w_meth3;
  set genic_w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=genic_w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=genic_w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_100U_0U_common _1Gy_100U_0U_common  _0Gy_100U_0U_common
       _01Gy_100U_0U_common_all _1Gy_100U_0U_common_all  _0Gy_100U_0U_common_all ;
  output out=mean_&inRegion._genic
  mean(_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean( _0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean( _0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;

proc sort data=mean_&inRegion._genic ;
  by grouped_pos2;
run;




proc means data=mean_&inRegion._genic noprint;
  by grouped_pos2  ;
  var  mean_01Gy_100U_0U_common mean_1Gy_100U_0U_common  mean_0Gy_100U_0U_common
       mean_01Gy_100U_0U_common_all mean_1Gy_100U_0U_common_all  mean_0Gy_100U_0U_common_all ;
  output out=mean_&inRegion._genic_2
  mean(mean_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common
  mean(mean_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common
  mean(mean_0Gy_100U_0U_common)=mean_0Gy_100U_0U_common
  mean(mean_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all
  mean(mean_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all
  mean(mean_0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all;
run;



data mean_&inRegion._genic_1;
  set mean_&inRegion._genic;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_&inRegion._genic_3;
  set mean_&inRegion._genic_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;




data te_w_meth2;
  set te_w_meth;
  distance_to_center=dar_Center-pos;
run;


data te_w_meth3;
  set te_w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=te_w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=te_w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy_100U_0U_common _1Gy_100U_0U_common  _0Gy_100U_0U_common
       _01Gy_100U_0U_common_all _1Gy_100U_0U_common_all  _0Gy_100U_0U_common_all ;
  output out=mean_&inRegion._te
  mean(_01Gy_100U_0U_common)=mean_01Gy_100U_0U_common_te
  mean(_1Gy_100U_0U_common)=mean_1Gy_100U_0U_common_te
  mean( _0Gy_100U_0U_common)=mean_0Gy_100U_0U_common_te
  mean(_01Gy_100U_0U_common_all)=mean_01Gy_100U_0U_common_all_te
  mean(_1Gy_100U_0U_common_all)=mean_1Gy_100U_0U_common_all_te
  mean( _0Gy_100U_0U_common_all)=mean_0Gy_100U_0U_common_all_te;
run;

proc sort data=mean_&inRegion._te ;
  by grouped_pos2;
run;


proc means data=mean_&inRegion._te noprint;
  by grouped_pos2  ;
  var  mean_01Gy_100U_0U_common_te mean_1Gy_100U_0U_common_te  mean_0Gy_100U_0U_common_te
       mean_01Gy_100U_0U_common_all_te mean_1Gy_100U_0U_common_all_te  mean_0Gy_100U_0U_common_all_te ;
  output out=mean_&inRegion._te_2
  mean(mean_01Gy_100U_0U_common_te)=mean_01Gy_100U_0U_common_te
  mean(mean_1Gy_100U_0U_common_te)=mean_1Gy_100U_0U_common_te
  mean(mean_0Gy_100U_0U_common_te)=mean_0Gy_100U_0U_common_te
  mean(mean_01Gy_100U_0U_common_all_te)=mean_01Gy_100U_0U_common_all_te
  mean(mean_1Gy_100U_0U_common_all_te)=mean_1Gy_100U_0U_common_all_te
  mean(mean_0Gy_100U_0U_common_all_te)=mean_0Gy_100U_0U_common_all_te;
run;


data mean_&inRegion._te_1;
  set mean_&inRegion._te;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_&inRegion._te_3;
  set mean_&inRegion._te_2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;



proc sort data=mean_&inRegion._genic_1;
   by pos;
proc sort data=mean_&inRegion._te_1;
   by pos;
proc sort data=mean_&inRegion._genic_3;
   by pos;
proc sort data=mean_&inRegion._te_3;
   by pos;
run;


data dmr_&inRegion._100;
  merge mean_&inRegion._genic_3 mean_&inRegion._te_3;
   by pos;
  keep pos mean_0Gy_100U_0U_common &methVar. mean_0Gy_100U_0U_common_te &methVar._te ;
run;

data dmr_&inRegion._100_2;
retain   pos mean_0Gy_100U_0U_common &methVar. mean_0Gy_100U_0U_common_te &methVar._te ;
  set dmr_&inRegion._100;
run;



proc export data=dmr_&inRegion._100_2
     outfile="!HOME/concannon/DTRA/dar_plot_&inRegion._100_bins.csv"
     dbms=csv replace;
run;


%mend;


%getMethyl(up_dar_01_72,mean_01Gy_100U_0U_common);
%getMethyl(up_dar_1_72,mean_1Gy_100U_0U_common);
%getMethyl(dn_dar_01_72,mean_01Gy_100U_0U_common);
%getMethyl(dn_dar_1_72,mean_1Gy_100U_0U_common);







