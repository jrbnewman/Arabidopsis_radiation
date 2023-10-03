/* Import BED files for each site and sample */

ods listing; ods html close;
libname wgbs '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname wgbsloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";

/* Merge methylation data with regions */

data pos2region;
  set wgbs.pos2region_index;
run;

proc sort data=pos2region;
  by chr pos gene_id region_type strand region_start region_stop;
run;

data region_cat;
   array gene[15] $15.;
   array type[15] $10.;
   array orientation[15] $1.;
   array start[15] $9.;
   array stop[15] $9.;
   retain gene1-gene15 type1-type15 orientation1-orientation15 start1-start15 stop1-stop15;
   set pos2region ;
   by chr pos;
   if first.pos then do;
       call missing(of gene1-gene15 );
       call missing(of type1-type15 );
       call missing(of orientation1-orientation15 );
       call missing(of start1-start15  );
       call missing(of stop1-stop15 );
       records=0;
       end;
   records + 1;
   gene[records]=gene_id;
   type[records]=region_type;
   orientation[records]=strand;
   start[records]=region_start;
   stop[records]=region_stop;
   if last.pos then output;
run;

data region_cat2;
  set region_cat;
  length gene_id_cat $240.;
  length region_type_cat $160.;
  length strand_cat $60.;
  length region_start_cat $160.;
  length region_stop_cat $160.;
  gene_id_cat = catx("|", OF gene1-gene15);
  region_type_cat = catx("|", OF type1-type15);
  strand_cat = catx("|", OF orientation1-orientation15);
  region_start_cat = catx("|", OF start1-start15);
  region_stop_cat = catx("|", OF stop1-stop15);
  keep chr pos gene_id_cat region_type_cat strand_cat
       region_start_cat region_stop_cat;
run;

data meth_data;
  set wgbsloc.methylation_data;
  drop pos;
  rename pos_end=pos;
run;

proc sort data=region_cat2;
   by chr pos;
proc sort data=meth_data;
  by chr pos;
run;

data meth_data_regions;
  merge region_cat2 (in=in1) meth_data (in=in2);
  by chr pos;
  if in1 and in2;
run;

data wgbsloc.meth_data2region;
  set meth_data_regions;
run;


