libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

/* Subset counts and regions by chromosome for binomial tests and HOMER annotations */

data meth_1 meth_2 meth_3 meth_4 meth_5;
  set wgbsloca.methylation_data_cg_chg_chh;
  if chr="1" then output meth_1;
  if chr="2" then output meth_2;
  if chr="3" then output meth_3;
  if chr="4" then output meth_4;
  if chr="5" then output meth_5;
run;


data acc_1 acc_2 acc_3 acc_4 acc_5;
  set wgbsloca.methylation_data_gc;
  if chr="1" then output acc_1;
  if chr="2" then output acc_2;
  if chr="3" then output acc_3;
  if chr="4" then output acc_4;
  if chr="5" then output acc_5;
run;

data wgbsloca.meth_data_cg_chg_chh_1;  set meth_1; run;
data wgbsloca.meth_data_cg_chg_chh_2;  set meth_2; run;
data wgbsloca.meth_data_cg_chg_chh_3;  set meth_3; run;
data wgbsloca.meth_data_cg_chg_chh_4;  set meth_4; run;
data wgbsloca.meth_data_cg_chg_chh_5;  set meth_5; run;


data wgbsloca.meth_data_gc_1;  set acc_1; run;
data wgbsloca.meth_data_gc_2;  set acc_2; run;
data wgbsloca.meth_data_gc_3;  set acc_3; run;
data wgbsloca.meth_data_gc_4;  set acc_4; run;
data wgbsloca.meth_data_gc_5;  set acc_5; run;


data dmr_1 dmr_2 dmr_3 dmr_4 dmr_5;
   set wgbsloca.cytosine_to_meth_region_ge2;
  if chr="1" then output dmr_1;
  if chr="2" then output dmr_2;
  if chr="3" then output dmr_3;
  if chr="4" then output dmr_4;
  if chr="5" then output dmr_5;
run;

data dar_1 dar_2 dar_3 dar_4 dar_5;
   set wgbsloca.cytosine_to_acc_region_ge2;
  if chr="1" then output dar_1;
  if chr="2" then output dar_2;
  if chr="3" then output dar_3;
  if chr="4" then output dar_4;
  if chr="5" then output dar_5;
run;


data wgbsloca.cytosine2meth_1; set dmr_1; run;
data wgbsloca.cytosine2meth_2; set dmr_2; run;
data wgbsloca.cytosine2meth_3; set dmr_3; run;
data wgbsloca.cytosine2meth_4; set dmr_4; run;
data wgbsloca.cytosine2meth_5; set dmr_5; run;

data wgbsloca.cytosine2acc_1; set dar_1; run;
data wgbsloca.cytosine2acc_2; set dar_2; run;
data wgbsloca.cytosine2acc_3; set dar_3; run;
data wgbsloca.cytosine2acc_4; set dar_4; run;
data wgbsloca.cytosine2acc_5; set dar_5; run;

/* Homer annotations */

data meth_regions;
  set wgbsloca.cytosine_to_meth_region_index;
  keep region_num site_type chr start_pos stop_pos comparison;
run;

data acc_regions;
  set wgbsloca.cytosine_to_acc_region_index;
  keep region_num site_type chr start_pos stop_pos comparison;
run;


data all_Regions;
  set meth_regions acc_regions;
run;

proc sort data=all_regions ;
   by comparison site_type region_num chr start_pos stop_pos;
proc means data=all_regions noprint;
   by comparison site_type region_num chr;
   var start_pos stop_pos;
   output out=all_regions2 min(start_pos)=start_pos max(stop_pos)=stop_pos;
run;

data homer_peaks;
  retain chr start_pos stop_pos region_id score strand;
   length region_id $100.;
   length score $1.;
   length strand $1.;
   set all_regions2;
   region_id=catx(".",comparison,site_type,chr,region_num);
   score=".";
   strand="+";
   keep  chr start_pos stop_pos region_id score strand;
run;


proc export data=homer_peaks
     outfile="/TB14/TB14/sandbox/dtra_sandbox/arabidopsis_regions.bed"
     dbms=tab replace; putnames=no;
run;

