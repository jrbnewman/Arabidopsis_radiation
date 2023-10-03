/* Export regions and sites for annotating with HOMER */

ods listing; ods html close;
libname cold "!PATCON/arabidopsis_wgbs_cold/sas_data";
libname coldloc "/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data";

proc datasets lib=work kill noprint;
run;
quit;


data cg_sites chg_sites chh_sites gc_sites;
  retain chr pos pos_end siteID score strand;
  set coldloc.methylation_data;
  length siteID $200.;
  length strand $1.;
  siteID=catx(":",site_type,chr,pos_end);
  score=".";
  strand="+";
  if site_type="CG" then output cg_sites;
  if site_type="CHG" then output chg_sites;
  if site_type="CHH" then output chh_sites;
  if site_type="GC" then output gc_sites;
  keep chr pos pos_end siteID score strand;
run;

proc sort data=cg_sites nodup;
  by chr pos pos_end siteID;
proc sort data=chg_sites nodup;
  by chr pos pos_end siteID;
proc sort data=chh_sites nodup;
  by chr pos pos_end siteID;
proc sort data=gc_sites nodup;
  by chr pos pos_end siteID;
run;


proc export data=cg_sites
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/CG_sites_for_HOMER.bed"
     dbms=tab replace; putnames=no;
run;

proc export data=chg_sites
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/CHG_sites_for_HOMER.bed"
     dbms=tab replace; putnames=no;
run;

proc export data=chh_sites
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/CHH_sites_for_HOMER.bed"
     dbms=tab replace; putnames=no;
run;

proc export data=gc_sites
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/GC_sites_for_HOMER.bed"
     dbms=tab replace; putnames=no;
run;

/* Create a site-to-window-to-superwindow index */

data pos2sub;
  set cold.cytosine_to_meth_region_index2 cold.cytosine_to_access_region_index2 ;
  keep pos region_num chr site_type ;
  rename region_num=subwindow_num;
run;

data pos2super;
  set cold.cytosine_2_meth_region_index_nd2 cold.cytosine_2_acc_region_index_nd2 ;
  keep pos region_num chr site_type ;
  rename region_num=superwindow_num;
run;


proc sort data=pos2sub;
  by chr pos site_type;
proc sort data=pos2super;
  by chr pos site_type;
run;


data sub2pos2super;
  merge pos2sub (in=in1) pos2super (in=in2);
  by chr pos site_type;
  if in1 and in2;
run;

/* For each window and superwindow I need the start and stop coordinates */

proc sort data=sub2pos2super;
  by site_type chr superwindow_num subwindow_num pos;
proc means data=sub2pos2super noprint;
  by site_type chr superwindow_num subwindow_num;
  var pos;
  output out=subwindow_coord min=window_start max=window_stop;
run;



proc sort data=sub2pos2super;
  by site_type chr superwindow_num pos;
proc means data=sub2pos2super noprint;
  by site_type chr superwindow_num ;
  var pos;
  output out=superwindow_coord min=window_start max=window_stop;
run;


data subwindow_coord2;
  retain chr window_start2 window_stop window_ID score strand;
  set subwindow_coord;
  length window_ID $50.;
  window_ID=catx(":",site_type,chr,window_start,window_stop);
  window_start2=window_Start-1;
  score=".";
  strand="+";
  keep chr window_start2 window_stop window_ID score strand;
run;

data superwindow_coord2;
  retain chr window_start2 window_stop superwindow_ID score strand;
  set superwindow_coord;
  length superwindow_ID $50.;
  superwindow_ID=catx(":",site_type,chr,window_start,window_stop);
  window_start2=window_Start-1;
  score=".";
  strand="+";
  keep chr window_start2 window_stop superwindow_ID score strand;
run;





proc export data=subwindow_coord2
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/all_windows_for_HOMER.bed"
     dbms=tab replace; putnames=no;
run;

proc export data=superwindow_coord2
     outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/all_superwindows_for_HOMER.bed"
     dbms=tab replace; putnames=no;
run;




