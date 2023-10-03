/* Create 3 regions per gene:  
    (1) Promoter: TSS-1500 to TSS+500
    (2) Gene body: TSS to TES
    (3) 3' downstream: TES-500 to TES+1500

    Then create an index of all possible sites within each region for merging methylation data */

ods listing; ods html close;
libname wgbs '!PATCON/arabidopsis_wgbs_cold/sas_data';

/* (1) Import gene coordinates */

proc import datafile="!PATCON/useful_arabidopsis_data/TAIR10/output/tair20_genes.bed"
   out=gene_coord dbms=tab replace;
   guessingrows=max; getnames=no;
run;

data gene_coord2;
  set gene_coord;
  rename VAR1=chr VAR2=gene_start VAR3=gene_stop VAR4=gene_id VAR6=strand;
  drop VAR5;
run;

/* (2) Create index of regions */

data promoter_regions;
   set gene_coord2;
   length region_type $10.;
   region_type="promoter";
   region_start=gene_start - 1500;
   region_stop = gene_start + 500;
   drop gene_start gene_stop;
run;

data genebody_regions;
   set gene_coord2;
   length region_type $10.;
   region_type="gene_body";
   region_start=gene_start;
   region_stop = gene_stop;
   drop gene_start gene_stop;
run;


data all_regions;
  set promoter_regions genebody_regions ;
run;

/* (3) Create index of positions per region */

proc sort data=all_regions;
   by gene_id region_type strand chr region_start region_stop;
run;

data position_index;
  set all_regions;
  by gene_id region_type strand chr region_start region_stop;
  do pos = region_start to region_stop ;
  output;
  end;
 run;


/* Make permenant */

data wgbs.pos2region_index;
  set position_index;
run;

data wgbs.promoter_regions;
  set promoter_regions;
run;

data wgbs.genebody_regions;
  set genebody_regions;
run;


data wgbs.gene_coordinates;
  set gene_coord2;
run;




/* Check: how frequently is each position assigned to a region? */

proc freq data=position_index noprint;
  tables chr*pos / out=pos_count;
proc sort data=pos_count;
  by descending count;
run; *15 total;

/* create a unique pos to region index */

/*
proc sort data=position_index;
  by chr pos gene_id region_type strand region_start region_stop;
run;

data region_cat;
   array gene[15] $15.;
   array type[15] $10.;
   array orientation[15] $1.;
   array start[15] $9.;
   array stop[15] $9.;
   retain gene1-gene15 type1-type15 orientation1-orientation15 start1-start15 stop1-stop15;
   set position_index ;
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
*/



