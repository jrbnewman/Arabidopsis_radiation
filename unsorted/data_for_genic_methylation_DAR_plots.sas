/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/Arabidopsis_thaliana.TAIR10.37.gtf"
  out=gtf dbms=tab replace;
  guessingrows=max;
  getnames=no;
run;


data gene;
  set gtf;
  where VAR3 = "gene";
  length gene_id $15.;
  gene_id=compress(tranwrd(tranwrd(scan(VAR9,2," "), ";", ""), '"', ''));
  keep gene_id VAR1 VAR4 VAR5 VAR7; 
  rename VAR7=strand VAR1=chr VAR4=gene_start VAR5=gene_stop;
run;


/* prep methylation and accessibility counts */
%macro filterData(siteType);


data meth_data;
  set arabMAP.methylation_data_CG_CHG_CHH;
  where site_Type = "&siteType.";
  keep site_type chr stop_pos treatment units rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by site_type chr stop_pos treatment  rep ;
proc means data=meth_data noprint;
  by site_type chr stop_pos treatment   ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;


proc transpose data=meth_data2 out=meth_sbys10;
  where total_C >= 10;
  by site_type chr stop_pos;
  id treatment  ;
  var perc_methyl;
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

proc sort data=gene;
  by chr gene_id;
run;

data gene_2;
  set gene;
  by chr gene_id;
  gene_length=gene_stop-gene_start;
  do pos = gene_start to gene_stop;
  output;
  end;
run;

proc sort data=gene_2;
  by chr pos;
proc sort data=meth_sbys10_2;
  by chr pos;
run;

data meth_sbys10_gene;
  merge meth_sbys10_2 (in=in1) gene_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data meth_sbys10_gene2;
  set meth_sbys10_gene;
  if strand="+" then pos_bin=int((pos-gene_start)/gene_length * 1000) + 1000 ;
  else if strand="-" then pos_bin=int((gene_stop-pos)/gene_length * 1000) + 1000 ;
run;

data promoter;
  set gene;
  if strand="+" then do;
    promoter_start=gene_start-1000;
    promoter_stop=gene_start-1;
    end;
  else do;
    promoter_start=gene_stop+1000;
    promoter_stop=gene_stop+1;
    end;
run;

data downstream;
  set gene;
  if strand="-" then do;
    downstream_start=gene_start-1;
    downstream_stop=gene_start-1000;
    end;
  else do;
    downstream_start=gene_stop+1;
    downstream_stop=gene_stop+1000;
    end;
run;

proc sort data=promoter;
  by chr gene_id ;
proc sort data=downstream;
  by chr gene_id ;
run;


data promoter_2;
  set promoter;
  by chr gene_id;
  if strand="+" then do pos = promoter_start to promoter_stop;
  output;
  end;
  else if strand="-" then do pos = promoter_stop to promoter_start;
  output;
  end;
run;


data downstream_2;
  set downstream;
  by chr gene_id;
  if strand="+" then do pos = downstream_start to downstream_stop;
  output;
  end;
  else if strand="-" then do pos = downstream_stop to downstream_start;
  output;
  end;
run;



proc sort data=promoter_2;
  by chr pos;
proc sort data=downstream_2;
  by chr pos;
run;

data meth_sbys10_promoter;
  merge meth_sbys10_2 (in=in1) promoter_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_downstream;
  merge meth_sbys10_2 (in=in1) downstream_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data meth_sbys10_promoter2;
  set meth_sbys10_promoter;
  if strand="+" then pos_bin=(pos-promoter_start);
  else if strand="-" then pos_bin=(promoter_start-pos) ;
run;

data meth_sbys10_downstream2;
  set meth_sbys10_downstream;
  if strand="+" then pos_bin=(pos-downstream_start) + 2000;
  else if strand="-" then pos_bin=(downstream_start-pos) + 2000;
run;



data meth_sbys10_genestack;
   set meth_sbys10_gene2 (in=in1) meth_sbys10_promoter2 (in=in2) meth_sbys10_downstream2 (in=in3);
run;

proc sort data=meth_sbys10_genestack;
  by gene_id chr pos_bin;
run;

/* two plots to make:
(1) Line plot of average accessibility across genic region, +/- 1kb (rescale gene length to 1000)
(2) For each genic region, calculate accessibiilty and then plot the distribution 

*/
* distrib of accessibilty per gene;

proc sort data= meth_sbys10_genestack;
   by gene_id;
run;

proc means data=meth_sbys10_genestack noprint;
  by gene_id;
  var _01Gy _1Gy _0Gy;
  output out=mean_acc_per_gene(drop=_TYPE_ _FREQ_) mean=;
run;

proc transpose data=mean_acc_per_gene out=mean_acc_tall(rename=(_NAME_=dose COL1=methylation));
  by gene_id;
  var  _01Gy _1Gy _0Gy;
run;


data mean_acc_tall2;
   set mean_acc_tall;
   length dose2 $10.;
   if dose = "_01Gy" then dose2="2_10cGy";
   if dose = "_0Gy" then dose2="1_Mock";
   if dose = "_1Gy" then dose2="3_100cGy";
   drop dose;
   rename dose2=dose;
run;

proc sort data=mean_acc_tall2;
  by dose gene_id;
run;


proc means data=mean_acc_tall2 noprint;
  by dose;
  var methylation;
  output out=distrib_acc mean=mean stddev=stddev min=min q1=q1 median=median q3=q3 max=max;
run;


proc export data=mean_acc_tall2 outfile="!HOME/concannon/DTRA/at_mean_&siteType._methylation_per_gene.csv"
dbms=csv replace;
run;




/* plot 2: data for lineplot across genic space */


data meth_data2;
  set meth_sbys10_genestack;
  grouped_pos=int(pos_bin/10)*10;
run;


proc sort data=meth_data2;
  by pos_bin grouped_pos;
run;


proc means data=meth_data2 noprint;
  by  pos_bin grouped_pos  ;
  var _01Gy _1Gy _0Gy   ;
  output out=meth_data3
  mean=;
run;


proc sort data=meth_data3 ;
  by grouped_pos;
run;


proc means data=meth_data3 noprint;
  by grouped_pos  ;
  var  _01Gy _1Gy _0Gy    
       ;
  output out=meth_data4
  mean=;
run;


data meth_data_export;
  set meth_data4;
  drop _TYPE_ _FREQ_ ;
  rename grouped_pos=pos;
run;


proc sort data=meth_data_export;
  by pos;
proc transpose data=meth_data_export out=mean_data_tall(rename=(_NAME_=dose COL1=methylation));
  by pos;
  var  _01Gy _1Gy _0Gy;
run;

data mean_data_tall2;
   set mean_data_tall;
   length dose2 $10.;
   if dose = "_01Gy" then dose2="2_10cGy";
   if dose = "_0Gy" then dose2="1_Mock";
   if dose = "_1Gy" then dose2="3_100cGy";
   drop dose;
   rename dose2=dose;
run;


proc sort data=mean_data_tall2;
  by dose pos;
run;


proc means data=mean_data_tall2 noprint;
  var methylation;
  output out=check min=min max=max;
run;

proc export data=mean_data_tall2 outfile="!HOME/concannon/DTRA/at_mean_&siteType._methylation_across_gene.csv"
dbms=csv replace;
run;



  
           
 /* plot 3: average methylation in intergenic vs genic DARs */


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER_annotation.txt"
   out=dar_annot dbms=tab replace;
   guessingrows=max;
run;

data dar_annot2;
  set dar_annot;
  length comparison $12.;
  length site_type $4.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'|'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'|'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'|'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'|')) + 0;
  feature = scan(annotation, 1, " ");
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;


data results_by_dar;
   set wgbsA.results_by_dar_5sites;
run;

proc freq data=results_by_dar;
  tables comparison*flag_fdr05;
run;


proc sort data=dar_annot2;
  by comparison site_type chr  region_num;
proc sort data=results_by_dar;
  by comparison site_type chr  region_num;
run;


data dar_w_annot;
  merge dar_annot2 (in=in1) results_by_dar (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


data dar_w_annot2;
  set dar_w_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;



/* methylation in DARs */


data meth_data;
  set arabMAP.methylation_data_CG_CHG_CHH;
  where site_Type="&siteType.";
  keep chr stop_pos treatment units rep total_C methyl_C perc_methyl;
run;

proc sort data=meth_data;
  by chr stop_pos treatment   rep ;
proc means data=meth_data noprint;
  by chr stop_pos treatment    ;
  var total_C methyl_C perc_methyl;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl)=perc_methyl;
run;

data meth_data3;
  set meth_data2;
  perc_methyl2=perc_methyl * 100 ;
run;



proc transpose data=meth_data3 out=meth_sbys10;
  where total_C >= 10;
  by  chr stop_pos;
  id treatment  ;
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



/* Get DARs */

data up_dar_01_1kb_ig up_dar_1_1kb_ig dn_dar_01_1kb_ig dn_dar_1_1kb_ig
     up_dar_01_1kb_gn up_dar_1_1kb_gn dn_dar_01_1kb_gn dn_dar_1_1kb_gn;
     set dar_w_annot2;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";

     dar_center=int((region_start + region_stop) / 2);
     plot_start=dar_center - 999;
     plot_stop=dar_center + 999;
     if feature="Intergenic" then do;
     /* all DMRs by site type and comparison */
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output up_dar_01_1kb_ig;
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output dn_dar_01_1kb_ig;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output up_dar_1_1kb_ig;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output dn_dar_1_1kb_ig;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output up_dar_01_1kb_ig;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output dn_dar_01_1kb_ig;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output up_dar_1_1kb_ig;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output dn_dar_1_1kb_ig;
     end;
     else do;
     /* all DMRs by site type and comparison */
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output up_dar_01_1kb_gn;
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output dn_dar_01_1kb_gn;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output up_dar_1_1kb_gn;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output dn_dar_1_1kb_gn;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output up_dar_01_1kb_gn;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output dn_dar_01_1kb_gn;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output up_dar_1_1kb_gn;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output dn_dar_1_1kb_gn;
     end;



     keep comparison chr dar_Center plot_start plot_stop ;
run;

proc sort data=up_dar_01_1kb_ig nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_01_1kb_ig nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_1_1kb_ig nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_1_1kb_ig nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_01_1kb_gn nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_01_1kb_gn nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=up_dar_1_1kb_gn nodup;  by chr dar_center plot_start plot_stop; run;
proc sort data=dn_dar_1_1kb_gn nodup;  by chr dar_center plot_start plot_stop; run;






%macro mergeMETH(inName);

data &inName._2;
  set &inName.;
  by chr dar_center;
  do pos = plot_start to plot_stop;
  output;
  end;
run;

proc sort data=&inName._2;
   by chr pos;
proc sort data=meth_sbys10_2;
   by chr pos;
run;

data &inName._w_meth;
  merge &inName._2 (in=in1) meth_sbys10_2 (in=in2);
  by chr pos;
  if in1 and in2;
run;


data &inName._w_meth2;
  set &inName._w_meth;
  distance_to_center=dar_Center-pos;
run;


data &inName._w_meth3;
  set &inName._w_meth2;
  grouped_pos=int(distance_to_center/10) * 10;
  if distance_to_center < 0 then grouped_pos2=grouped_pos -10;
  else grouped_pos2=grouped_pos;
run;

/* For each, calculate the mean and SD */

proc sort data=&inName._w_meth3;
  by distance_to_center grouped_pos2 ;
run;


proc means data=&inName._w_meth3 noprint;
  by distance_to_center grouped_pos2  ;
  var  _01Gy _0Gy _1Gy _01Gy_common _1Gy_common  _0Gy_common
       _01Gy_common_all _1Gy_common_all  _0Gy_common_all ;
  output out=mean_diff_&inName.
  mean(_01Gy)=mean_01Gy
  mean(_1Gy)=mean_1Gy
  mean( _0Gy)=mean_0Gy
  mean(_01Gy_common)=mean_01Gy_common
  mean(_1Gy_common)=mean_1Gy_common
  mean( _0Gy_common)=mean_0Gy_common
  mean(_01Gy_common_all)=mean_01Gy_common_all
  mean(_1Gy_common_all)=mean_1Gy_common_all
  mean( _0Gy_common_all)=mean_0Gy_common_all;
run;

proc sort data=mean_diff_&inName. ;
  by grouped_pos2;
run;


proc means data=mean_diff_&inName. noprint;
  by grouped_pos2  ;
  var  mean_01Gy mean_1Gy mean_0Gy mean_01Gy_common mean_1Gy_common  mean_0Gy_common
       mean_01Gy_common_all mean_1Gy_common_all  mean_0Gy_common_all ;
  output out=mean_diff_&inName._2
  mean(mean_01Gy)=mean_01Gy
  mean(mean_1Gy)=mean_1Gy
  mean( mean_0Gy)=mean_0Gy
  mean(mean_01Gy_common)=mean_01Gy_common
  mean(mean_1Gy_common)=mean_1Gy_common
  mean(mean_0Gy_common)=mean_0Gy_common
  mean(mean_01Gy_common_all)=mean_01Gy_common_all
  mean(mean_1Gy_common_all)=mean_1Gy_common_all
  mean(mean_0Gy_common_all)=mean_0Gy_common_all;
run;


data mean_diff_&inName._1;
  set mean_diff_&inName.;
  drop _TYPE_ _FREQ_;
  keep distance_To_center mean_: ;
  rename distance_to_center=pos;
run;

data mean_diff_&inName._3;
  set mean_diff_&inName._2;
  drop _TYPE_ _FREQ_;
  rename grouped_pos2=pos;
run;
%mend;



%mergeMETH(up_dar_01_1kb_ig);
%mergeMETH(dn_dar_01_1kb_ig);
%mergeMETH(up_dar_1_1kb_ig);
%mergeMETH(dn_dar_1_1kb_ig);

%mergeMETH(up_dar_01_1kb_gn);
%mergeMETH(dn_dar_01_1kb_gn);
%mergeMETH(up_dar_1_1kb_gn);
%mergeMETH(dn_dar_1_1kb_gn);



/* Export for plots --- Make 0.1 and 1Gy comparisons separately for now */

%macro exportLine(inData, var1, var2, outName);

data export;
  retain pos &var1. &var2.;
  set &inData.;
  keep pos &var1. &var2.;
run;

proc export data=export
     outfile="!HOME/concannon/DTRA/at_rad_DAR_line_plots/input_data/DAR_&siteType._methylation_&outName._new.csv"
     dbms=csv replace;
run;

%mend;

%exportLine( mean_diff_up_dar_01_1kb_ig_3, mean_0Gy, mean_01Gy, hyper_DAR_1kb_01Gy_0Gy_notCommon_binned_intergenic);
%exportLine( mean_diff_up_dar_01_1kb_ig_3, mean_0Gy_common, mean_01Gy_common, hyper_DAR_1kb_01Gy_0Gy_binned_intergenic);
%exportLine( mean_diff_up_dar_01_1kb_ig_3, mean_0Gy_common, mean_01Gy_common_all, hyper_DAR_1kb_01Gy_0Gy_common_binned_intergenic);

%exportLine( mean_diff_up_dar_1_1kb_ig_3, mean_0Gy, mean_1Gy, hyper_DAR_1kb_1Gy_0Gy_notCommon_binned_intergenic);
%exportLine( mean_diff_up_dar_1_1kb_ig_3, mean_0Gy_common, mean_1Gy_common, hyper_DAR_1kb_1Gy_0Gy_binned_intergenic);
%exportLine( mean_diff_up_dar_1_1kb_ig_3, mean_0Gy_common, mean_1Gy_common_all, hyper_DAR_1kb_1Gy_0Gy_common_binned_intergenic);

%exportLine( mean_diff_dn_dar_01_1kb_ig_3, mean_0Gy, mean_01Gy, hypo_DAR_1kb_01Gy_0Gy_notCommon_binned_intergenic);
%exportLine( mean_diff_dn_dar_01_1kb_ig_3, mean_0Gy_common, mean_01Gy_common, hypo_DAR_1kb_01Gy_0Gy_binned_intergenic);
%exportLine( mean_diff_dn_dar_01_1kb_ig_3, mean_0Gy_common, mean_01Gy_common_all, hypo_DAR_1kb_01Gy_0Gy_common_binned_intergenic);

%exportLine( mean_diff_dn_dar_1_1kb_ig_3, mean_0Gy, mean_1Gy, hypo_DAR_1kb_1Gy_0Gy_notCommon_binned_intergenic);
%exportLine( mean_diff_dn_dar_1_1kb_ig_3, mean_0Gy_common, mean_1Gy_common, hypo_DAR_1kb_1Gy_0Gy_binned_intergenic);
%exportLine( mean_diff_dn_dar_1_1kb_ig_3, mean_0Gy_common, mean_1Gy_common_all, hypo_DAR_1kb_1Gy_0Gy_common_binned_intergenic);

%exportLine( mean_diff_up_dar_01_1kb_gn_3, mean_0Gy, mean_01Gy, hyper_DAR_1kb_01Gy_0Gy_notCommon_binned_genic);
%exportLine( mean_diff_up_dar_01_1kb_gn_3, mean_0Gy_common, mean_01Gy_common, hyper_DAR_1kb_01Gy_0Gy_binned_genic);
%exportLine( mean_diff_up_dar_01_1kb_gn_3, mean_0Gy_common, mean_01Gy_common_all, hyper_DAR_1kb_01Gy_0Gy_common_binned_genic);

%exportLine( mean_diff_up_dar_1_1kb_gn_3, mean_0Gy, mean_1Gy, hyper_DAR_1kb_1Gy_0Gy_notCommon_binned_genic);
%exportLine( mean_diff_up_dar_1_1kb_gn_3, mean_0Gy_common, mean_1Gy_common, hyper_DAR_1kb_1Gy_0Gy_binned_genic);
%exportLine( mean_diff_up_dar_1_1kb_gn_3, mean_0Gy_common, mean_1Gy_common_all, hyper_DAR_1kb_1Gy_0Gy_common_binned_genic);

%exportLine( mean_diff_dn_dar_01_1kb_gn_3, mean_0Gy, mean_01Gy, hypo_DAR_1kb_01Gy_0Gy_notCommon_binned_genic);
%exportLine( mean_diff_dn_dar_01_1kb_gn_3, mean_0Gy_common, mean_01Gy_common, hypo_DAR_1kb_01Gy_0Gy_binned_genic);
%exportLine( mean_diff_dn_dar_01_1kb_gn_3, mean_0Gy_common, mean_01Gy_common_all, hypo_DAR_1kb_01Gy_0Gy_common_binned_genic);

%exportLine( mean_diff_dn_dar_1_1kb_gn_3, mean_0Gy, mean_1Gy, hypo_DAR_1kb_1Gy_0Gy_notCommon_binned_genic);
%exportLine( mean_diff_dn_dar_1_1kb_gn_3, mean_0Gy_common, mean_1Gy_common, hypo_DAR_1kb_1Gy_0Gy_binned_genic);
%exportLine( mean_diff_dn_dar_1_1kb_gn_3, mean_0Gy_common, mean_1Gy_common_all, hypo_DAR_1kb_1Gy_0Gy_common_binned_genic);

%mend;


%filterData(CG);
%filterData(CHG);
%filterData(CHH);
























