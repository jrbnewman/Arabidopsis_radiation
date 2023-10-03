libname wgbsA '!PATCON/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


data dmrs;
   set wgbslocA.cytosine_to_meth_region_ge2_v4;
   length regionID $100.;
   regionID = catx("_",site_type,comparison,chr,region_num);
   keep site_type comparison chr start_pos stop_pos regionID;
run;

proc sort data=dmrs;
  by site_type comparison chr regionID;
proc means data=dmrs noprint;
   by site_type comparison chr regionID;
   var start_pos stop_pos;
   output out=dmr_regions min(start_pos)=start_pos max(stop_pos)=stop_pos;
run;

proc sort data=dmr_regions;
  by chr start_pos stop_pos comparison regionID;
run;



data dars;
   set wgbslocA.cytosine_to_acc_region_ge2_v4;
   length regionID $100.;
   regionID = catx("_",site_type,comparison,chr,region_num);
   keep site_type comparison chr start_pos stop_pos regionID;
run;

proc sort data=dars;
  by site_type comparison chr regionID;
proc means data=dars noprint;
   by site_type comparison chr regionID;
   var start_pos stop_pos;
   output out=dar_regions min(start_pos)=start_pos max(stop_pos)=stop_pos;
run;

proc sort data=dar_regions;
  by chr start_pos stop_pos comparison regionID;
run;



%macro exportSite(minSite);

data dmr_cg_01 dmr_cg_1 
     dmr_chg_01 dmr_chg_1
     dmr_chh_01 dmr_chh_1;
     retain chr start_pos stop_pos regionID score strand;
     set dmr_regions;
     where _FREQ_ >= &minSite.;
     score=0;
     strand="+";
     if site_type = "CG" and comparison = "0Gy_vs_01G" then output dmr_cg_01;
     if site_type = "CG" and comparison = "0Gy_vs_1Gy" then output dmr_cg_1;
     if site_type = "CHG" and comparison = "0Gy_vs_01G" then output dmr_chg_01;
     if site_type = "CHG" and comparison = "0Gy_vs_1Gy" then output dmr_chg_1;
     if site_type = "CHH" and comparison = "0Gy_vs_01G" then output dmr_chh_01;
     if site_type = "CHH" and comparison = "0Gy_vs_1Gy" then output dmr_chh_1;
     keep chr start_pos stop_pos regionID score strand;
run;


data dar_gc_01 dar_gc_1;
     retain chr start_pos stop_pos regionID score strand;
   set dar_regions;
     where _FREQ_ >= &minSite.;
     score=0;
     strand="+";
     if site_type = "GC" and comparison = "01Gy_0Gy" then output dar_gc_01;
     if site_type = "GC" and comparison = "1Gy_0Gy" then output dar_gc_1;
     keep chr start_pos stop_pos regionID score strand;
run;


proc export data=dmr_cg_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_10cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;
proc export data=dmr_cg_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_100cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;

proc export data=dmr_chg_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_10cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;
proc export data=dmr_chg_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_100cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;

proc export data=dmr_chh_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHH_10cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;
proc export data=dmr_chh_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DMR_CHH_100cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;

proc export data=dar_gc_01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_10cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;
proc export data=dar_gc_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_100cGy_minSites_&minSite..bed" dbms=tab replace; putnames=no; run;


%mend;

%exportSite(2);
%exportSite(3);
%exportSite(5);
%exportSite(10);

