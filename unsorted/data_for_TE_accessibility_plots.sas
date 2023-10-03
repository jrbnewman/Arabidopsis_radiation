/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


proc import datafile="!HOME/concannon/DTRA/TAIR10_transposable_elements.noChr.gtf"
  out=gtf dbms=tab replace;
  guessingrows=max;
  getnames=no;
run;


data gene;
  set gtf;
  where VAR3 = "exon";
  length chr $2.;
  length gene_id $15.;
  gene_id=compress(tranwrd(tranwrd(scan(VAR9,2," "), ";", ""), '"', ''));
  chr=compress(VAR1);
  keep gene_id chr  VAR4 VAR5 VAR7; 
  rename VAR7=strand  VAR4=gene_start VAR5=gene_stop;
run;


/* prep methylation and accessibility counts */



data meth_data;
  set arabMAP.methylation_data_gc;
  keep site_type chr stop_pos treatment units rep total_C methyl_C perc_methyl_norm;
run;

proc sort data=meth_data;
  by site_type chr stop_pos treatment units  rep ;
proc means data=meth_data noprint;
  by site_type chr stop_pos treatment  units  ;
  var total_C methyl_C perc_methyl_norm;
  output out=meth_data2 sum(total_C)=total_C sum(methyl_C)=methyl_C mean(perc_methyl_norm)=perc_methyl;
run;


proc transpose data=meth_data2 out=meth_sbys10;
  where total_C >= 10;
  by site_type chr stop_pos;
  id treatment units ;
  var perc_methyl;
run;

data meth_sbys10_2;
  set meth_sbys10;
  if _01Gy100U = . or _01Gy0U = . then _01Gy=.; else _01Gy = _01Gy100U - _01Gy0U;
  if _1Gy100U = . or _1Gy0U = . then _1Gy=.; else _1Gy = _1Gy100U - _1Gy0U;
  if _0Gy100U = . or _0Gy0U = . then _0Gy=.; else _0Gy = _0Gy100U - _0Gy0U;


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

proc transpose data=mean_acc_per_gene out=mean_acc_tall(rename=(_NAME_=dose COL1=accessibility));
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
  var accessibility;
  output out=distrib_acc mean=mean stddev=stddev min=min q1=q1 median=median q3=q3 max=max;
run;


proc export data=mean_acc_tall2 outfile="!HOME/concannon/DTRA/at_mean_accessiblity_per_gene.csv"
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
proc transpose data=meth_data_export out=mean_data_tall(rename=(_NAME_=dose COL1=accessibility));
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
  var accessibility;
  output out=check min=min max=max;
run;

proc export data=mean_data_tall2 outfile="!HOME/concannon/DTRA/at_mean_accessiblity_across_TEs.csv"
dbms=csv replace;
run;


/* Methylation */



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

proc export data=mean_data_tall2 outfile="!HOME/concannon/DTRA/at_mean_&siteType._methylation_across_TEs.csv"
dbms=csv replace;
run;


%mend;


%filterData(CG);
%filterData(CHG);
%filterData(CHH);

























