/* Import BED files for each site and sample */

ods listing; ods html close;
libname wgbs '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname wgbsloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";

/* Merge methylation data with regions -- cold MapIT data only */

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
  if temperature ne "FANS";
  drop pos;
  rename pos_end=pos;
run;

/* Any site with good coverage per rep */
data good_coverage;
   set wgbs.flag_coverage_ge10;
  if site_type="GC" then do;
     if flag_22C_0U_coverage_ge10=1 and flag_22C_100U_coverage_ge10=1 and
        flag_4C_0U_coverage_ge10=1 and flag_4C_100U_coverage_ge10=1 then output;
     end;
   else do;
        if flag_22C_0U_coverage_ge10=1 and flag_4C_0U_coverage_ge10=1 then output;
   end;
 keep chr pos site_type;
run;

/* Any site with evidence of methylation */
data has_meth;
  set wgbs.flag_methylation_gt0;
  if flag_22C_0U_methylation_gt0=1
  or flag_22C_100U_methylation_gt0=1
  or flag_4C_0U_methylation_gt0=1
  or flag_4C_100U_methylation_gt0=1;
 keep chr pos site_type;
run;

proc sort data=meth_data;
  by chr pos site_type;
proc sort data=good_coverage;
  by chr pos site_type;
proc sort data=has_meth;
  by chr pos site_type;
run;

data meth_data2;
  merge good_coverage (in=in1) has_meth (in=in2) meth_data (in=in3);
  by chr pos site_type;
  if in1 and in2 and in3;
run;

proc sort data=region_cat2;
   by chr pos;
proc sort data=meth_data2;
   by chr pos;
run;

data meth_data_regions;
  merge region_cat2 (in=in1) meth_data2 (in=in2);
  by chr pos;
  if in1 and in2 then output meth_data_regions;
run;

/* Count and flag regions with at least 3 testable positions per region */

data flag_testable;
  set meth_data_regions;
  if perc_methyl > 0 then flag_testable=1;
  else flag_testable=0;
  keep chr gene_id_cat region_type_cat strand_cat site_type region_start_cat
       region_stop_cat
       dose units rep pos flag_testable;
run;

/* Per position and site type, count the number of testable observations
   Should be max 6 for non-GC sites, and max 12 for GC sites */
  
proc sort data=flag_testable;
   by chr gene_id_cat region_type_cat strand_cat site_type region_start_cat
       region_stop_cat pos;
proc means data=flag_testable noprint;
   by chr gene_id_cat region_type_cat strand_cat site_type region_start_cat
       region_stop_cat pos;
   var flag_testable;
   output out=flag_testable_by_pos
          sum=;
run;

data flag_testable_by_pos2;
  set flag_testable_by_pos;
  if flag_testable > 0 then flag_testable2=1;
  else flag_testable2=0;
run;

proc freq data=flag_testable_by_pos2;
  tables flag_testable flag_testable2;
run;
 

/* Expand out regions */

data flag_testable_by_pos3;
  set flag_testable_by_pos2;
  length gene_id $15.;
  length region_type $10.;
  length strand $1.;
  format region_start best12.;
  format region_stop best12.;
  do i=1 by 1 while(scan(gene_id_cat,i,"|") ^= "");
      gene_id=scan(gene_id_cat,i,"|");
      region_type=scan(region_type_cat,i,"|");
      strand=scan(strand_cat,i,"|");
      region_start=scan(region_start_cat,i,"|");
      region_stop=scan(region_stop_cat,i,"|");
      output;
      end;
run;

/* Count number of testable positions per region, and flag if less than 3 */

proc sort data=flag_testable_by_pos3;
  by chr gene_id region_type strand region_start region_stop site_type;
proc means data=flag_testable_by_pos3 noprint;
  by chr gene_id region_type strand region_start region_stop site_type;
  var flag_testable2;
  output out=num_testable_pos_by_region sum=num_testable;
run;




data num_testable_pos_by_region2;
  set num_testable_pos_by_region;
  if num_testable >=3 then flag_num_testable_ge3=1;
  else flag_num_testable_ge3=0;
run;



proc freq data=num_testable_pos_by_region2;
  tables site_type*flag_num_testable_ge3;
run;

/*

*/
 

/* Make perm */

data wgbs.num_testable_pos_by_region;
  set num_testable_pos_by_region2;
run;

data wgbsloc.meth_data2region;
  set meth_data_regions;
run;


/* FANS data */


  


data meth_data;
  set wgbsloc.methylation_data;
  if temperature = "FANS" then output;
  if temperature = "22C" and units="100U" then output;
  drop pos;
  rename pos_end=pos;
run;

/* Any site with good coverage per rep */
data good_coverage;
   set wgbs.flag_coverage_ge10;
   if flag_22C_100U_coverage_ge10=1 and (flag_FANS_0p5U_coverage_ge10=1
    or flag_FANS_1p5U_coverage_ge10=1 or flag_FANS_5U_coverage_ge10=1 
    or flag_FANS_25U_coverage_ge10 =1) then output;
 keep chr pos site_type;
run;

/* Any site with evidence of methylation */
data has_meth;
  set wgbs.flag_methylation_fans_gt0;
  if flag_22C_100U_methylation_gt0=1
  or flag_FANS_0p5U_methylation_gt0=1
  or flag_FANS_1p5U_methylation_gt0=1
  or flag_FANS_5U_methylation_gt0=1
  or flag_FANS_25U_methylation_gt0=1;
 keep chr pos site_type;
run;

proc sort data=meth_data;
  by chr pos site_type;
proc sort data=good_coverage;
  by chr pos site_type;
proc sort data=has_meth;
  by chr pos site_type;
run;

data meth_data2;
  merge good_coverage (in=in1) has_meth (in=in2) meth_data (in=in3);
  by chr pos site_type;
  if in1 and in2 and in3;
run;

proc sort data=region_cat2;
   by chr pos;
proc sort data=meth_data2;
   by chr pos;
run;

data meth_data_regions;
  merge region_cat2 (in=in1) meth_data2 (in=in2);
  by chr pos;
  if in1 and in2 then output meth_data_regions;
run;

/* Count and flag regions with at least 3 testable positions per region */

data flag_testable;
  set meth_data_regions;
  if perc_methyl > 0 then flag_testable=1;
  else flag_testable=0;
  keep chr gene_id_cat region_type_cat strand_cat site_type region_start_cat
       region_stop_cat
       dose units rep pos flag_testable;
run;

/* Per position and site type, count the number of testable observations
   Should be max 6 for non-GC sites, and max 12 for GC sites */
  
proc sort data=flag_testable;
   by chr gene_id_cat region_type_cat strand_cat site_type region_start_cat
       region_stop_cat pos;
proc means data=flag_testable noprint;
   by chr gene_id_cat region_type_cat strand_cat site_type region_start_cat
       region_stop_cat pos;
   var flag_testable;
   output out=flag_testable_by_pos
          sum=;
run;

data flag_testable_by_pos2;
  set flag_testable_by_pos;
  if flag_testable > 0 then flag_testable2=1;
  else flag_testable2=0;
run;

proc freq data=flag_testable_by_pos2;
  tables flag_testable flag_testable2;
run;
 

/* Expand out regions */

data flag_testable_by_pos3;
  set flag_testable_by_pos2;
  length gene_id $15.;
  length region_type $10.;
  length strand $1.;
  format region_start best12.;
  format region_stop best12.;
  do i=1 by 1 while(scan(gene_id_cat,i,"|") ^= "");
      gene_id=scan(gene_id_cat,i,"|");
      region_type=scan(region_type_cat,i,"|");
      strand=scan(strand_cat,i,"|");
      region_start=scan(region_start_cat,i,"|");
      region_stop=scan(region_stop_cat,i,"|");
      output;
      end;
run;

/* Count number of testable positions per region, and flag if less than 3 */

proc sort data=flag_testable_by_pos3;
  by chr gene_id region_type strand region_start region_stop site_type;
proc means data=flag_testable_by_pos3 noprint;
  by chr gene_id region_type strand region_start region_stop site_type;
  var flag_testable2;
  output out=num_testable_pos_by_region sum=num_testable;
run;




data num_testable_pos_by_region2;
  set num_testable_pos_by_region;
  if num_testable >=3 then flag_num_testable_ge3=1;
  else flag_num_testable_ge3=0;
run;



proc freq data=num_testable_pos_by_region2;
  tables site_type*flag_num_testable_ge3;
run;

/*

*/
 

/* Make perm */

data wgbs.num_testable_pos_by_region_fans;
  set num_testable_pos_by_region2;
run;

data wgbsloc.meth_data2region_fans;
  set meth_data_regions;
run;


