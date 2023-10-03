/* Counts for arabidopsis paper */

libname brassMAP '/TB14/TB14/sandbox/dtra_sandbox/brass_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

/* Import bedGRaphs and expression TSVs and average in 100bp and 1kb windows */

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CG_10cGy_0cGy_1h_analyzed.bedGraph"
   out=cg_10cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CG_10cGy_0cGy_72h_analyzed.bedGraph"
   out=cg_10cgy_72h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CG_1p4cGy_0cGy_1h_analyzed.bedGraph"
   out=cg_1p4cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CG_1p4cGy_0cGy_72h_analyzed.bedGraph"
   out=cg_1p4cgy_72h dbms=tab replace; getnames=no;
run;




proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHG_10cGy_0cGy_1h_analyzed.bedGraph"
   out=chg_10cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHG_10cGy_0cGy_72h_analyzed.bedGraph"
   out=chg_10cgy_72h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHG_1p4cGy_0cGy_1h_analyzed.bedGraph"
   out=chg_1p4cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHG_1p4cGy_0cGy_72h_analyzed.bedGraph"
   out=chg_1p4cgy_72h dbms=tab replace; getnames=no;
run;


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHH_10cGy_0cGy_1h_analyzed.bedGraph"
   out=chh_10cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHH_10cGy_0cGy_72h_analyzed.bedGraph"
   out=chh_10cgy_72h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHH_1p4cGy_0cGy_1h_analyzed.bedGraph"
   out=chh_1p4cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_CHH_1p4cGy_0cGy_72h_analyzed.bedGraph"
   out=chh_1p4cgy_72h dbms=tab replace; getnames=no;
run;




proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_GC_10cGy_0cGy_1h_analyzed.bedGraph"
   out=gc_10cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_GC_10cGy_0cGy_72h_analyzed.bedGraph"
   out=gc_10cgy_72h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_GC_1p4cGy_0cGy_1h_analyzed.bedGraph"
   out=gc_1p4cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Brassica_rapa_GC_1p4cGy_0cGy_72h_analyzed.bedGraph"
   out=gc_1p4cgy_72h dbms=tab replace; getnames=no;
run;



/* Macro to reduce data into bins of 100bp or 1000bp */

%macro windowData(inData, outData);

data bin_100_1000;
   set &inData.;
   where var1 ? "A";
   bin_100 = int(var3/100) + 1;
   bin_1000 = int(var3/1000) + 1;
run;


proc sort data=bin_100_1000;
  by var1 bin_100;
proc means data=bin_100_1000 noprint;
  by var1 bin_100;
  var VAR2 VAR3 VAR4;
  output out=mean_bin_100 min(var2)=var2 max(var3)=var3 mean(VAR4)=var4;
run;

proc means data=bin_100_1000 noprint;
  by var1 bin_1000;
  var var2 var3 VAR4;
  output out=mean_bin_1000 min(var2)=var2 max(var3)=var3 mean(VAR4)=var4;
run;


data mean_bin_100_2;
  set mean_bin_100;
*  bin_start = (bin_100 - 1 ) * 100;
*  bin_stop = bin_100  * 100;
run;

data mean_bin_1000_2;
  set mean_bin_1000;
*  bin_start = (bin_1000 - 1 ) * 1000;
*  bin_stop = bin_1000  * 1000;
run;

data export_100;
*  retain VAR1 bin_start bin_stop VAR4;
*  set mean_bin_100_2;
*  keep VAR1 bin_start bin_stop VAR4;
  retain VAR1 var2 var3 VAR4;
  set mean_bin_100_2;
  keep VAR1 var2 var3 VAR4;
run;

data export_1000;
  *retain VAR1 bin_start bin_stop VAR4;
  *set mean_bin_1000_2;
  *keep VAR1 bin_start bin_stop VAR4;
  retain VAR1 var2 var3 VAR4;
  set mean_bin_1000_2;
  keep VAR1 var2 var3 VAR4;
run;

proc export data=export_100
     outfile="/TB14/TB14/sandbox/dtra_sandbox/&outData._100bp_window_2.txt"
     dbms=tab replace; putnames=no;
run;

proc export data=export_1000
     outfile="/TB14/TB14/sandbox/dtra_sandbox/&outData._1kb_window_2.txt"
     dbms=tab replace; putnames=no;
run;

%mend;



%windowData(cg_10cgy_1h, Br_CG_Meth_10cGy_1h);
%windowData(cg_10cgy_72h, Br_CG_Meth_10cGy_72h);
%windowData(cg_1p4cgy_1h, Br_CG_Meth_1p4cGy_1h);
%windowData(cg_1p4cgy_72h, Br_CG_Meth_1p4cGy_72h);


%windowData(chg_10cgy_1h, Br_CHG_Meth_10cGy_1h);
%windowData(chg_10cgy_72h, Br_CHG_Meth_10cGy_72h);
%windowData(chg_1p4cgy_1h, Br_CHG_Meth_1p4cGy_1h);
%windowData(chg_1p4cgy_72h, Br_CHG_Meth_1p4cGy_72h);


%windowData(chh_10cgy_1h, Br_CHH_Meth_10cGy_1h);
%windowData(chh_10cgy_72h, Br_CHH_Meth_10cGy_72h);
%windowData(chh_1p4cgy_1h, Br_CHH_Meth_1p4cGy_1h);
%windowData(chh_1p4cgy_72h, Br_CHH_Meth_1p4cGy_72h);


%windowData(gc_10cgy_1h, Br_GC_Acc_10cGy_1h);
%windowData(gc_10cgy_72h, Br_GC_Acc_10cGy_72h);
%windowData(gc_1p4cgy_1h, Br_GC_Acc_1p4cGy_1h);
%windowData(gc_1p4cgy_72h, Br_GC_Acc_1p4cGy_72h);



/* Set up sig sites */

/* Make tile tracks for DEG, DMR, DAR */

data dar_01_1 dar_1p4_1 dar_01_72 dar_1p4_72;
     set brassMAP.results_by_dar_annot;

     /* FET */
     if comparison="10cGy_0cGy_1h" then output dar_01_1;
     if comparison="1p4cGy_0cGy_1h" then output dar_1p4_1;
     if comparison="10cGy_0cGy_72h" then output dar_01_72;
     if comparison="1p4cGy_0cGy_72h" then output dar_1p4_72;

     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison="10cGy_0cGy_1h" then output dar_01_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison="1p4cGy_0cGy_1h" then output dar_1p4_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison="10cGy_0cGy_72h" then output dar_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison="1p4cGy_0cGy_72h" then output dar_1p4_72;
     keep chr region_start region_stop;
run;

data cg_dmr_01_1 cg_dmr_1p4_1 cg_dmr_01_72 cg_dmr_1p4_72 ;
     set brassMAP.results_by_dmr_annot;
     where site_type="CG";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_10cGy_1h" then output cg_dmr_01_1;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_1p4cGy_1h" then output cg_dmr_1p4_1;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_10cGy_72h" then output cg_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_1p4cGy_72h" then output cg_dmr_1p4_72;

     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_10cGy_1h" then output cg_dmr_01_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_1p4cGy_1h" then output cg_dmr_1p4_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_10cGy_72h" then output cg_dmr_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_1p4cGy_72h" then output cg_dmr_1p4_72;
     keep chr region_start region_stop;
run;

data chg_dmr_01_1 chg_dmr_1p4_1 chg_dmr_01_72 chg_dmr_1p4_72 ;
     set brassMAP.results_by_dmr_annot;
     where site_type="CHG";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_10cGy_1h" then output chg_dmr_01_1;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_1p4cGy_1h" then output chg_dmr_1p4_1;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_10cGy_72h" then output chg_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_1p4cGy_72h" then output chg_dmr_1p4_72;

     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_10cGy_1h" then output chg_dmr_01_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_1p4cGy_1h" then output chg_dmr_1p4_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_10cGy_72h" then output chg_dmr_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_1p4cGy_72h" then output chg_dmr_1p4_72;
     keep chr region_start region_stop;
run;

data chh_dmr_01_1 chh_dmr_1p4_1 chh_dmr_01_72 chh_dmr_1p4_72 ;
     set brassMAP.results_by_dmr_annot;
     where site_type="CHH";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_10cGy_1h" then output chh_dmr_01_1;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_1p4cGy_1h" then output chh_dmr_1p4_1;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_10cGy_72h" then output chh_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0cGy_vs_1p4cGy_72h" then output chh_dmr_1p4_72;

     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_10cGy_1h" then output chh_dmr_01_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_1p4cGy_1h" then output chh_dmr_1p4_1;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_10cGy_72h" then output chh_dmr_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0cGy_vs_1p4cGy_72h" then output chh_dmr_1p4_72;
     keep chr region_start region_stop;
run;






proc sort data=dar_01_72 nodup; by chr region_start region_stop; run;
proc sort data=dar_1p4_72 nodup; by chr region_start region_stop; run;
proc sort data=cg_dmr_01_72 nodup; by chr region_start region_stop; run;
proc sort data=cg_dmr_1p4_72 nodup; by chr region_start region_stop; run;
proc sort data=chg_dmr_01_72 nodup; by chr region_start region_stop; run;
proc sort data=chg_dmr_1p4_72 nodup; by chr region_start region_stop; run;
proc sort data=chh_dmr_01_72 nodup; by chr region_start region_stop; run;
proc sort data=chh_dmr_1p4_72 nodup; by chr region_start region_stop; run;

proc sort data=dar_01_1 nodup; by chr region_start region_stop; run;
proc sort data=dar_1p4_1 nodup; by chr region_start region_stop; run;
proc sort data=cg_dmr_01_1 nodup; by chr region_start region_stop; run;
proc sort data=cg_dmr_1p4_1 nodup; by chr region_start region_stop; run;
proc sort data=chg_dmr_01_1 nodup; by chr region_start region_stop; run;
proc sort data=chg_dmr_1p4_1 nodup; by chr region_start region_stop; run;
proc sort data=chh_dmr_01_1 nodup; by chr region_start region_stop; run;
proc sort data=chh_dmr_1p4_1 nodup; by chr region_start region_stop; run;




/* Macro to subset only sig regions and also to reduce data into bins of 100bp or 1000bp */

%macro windowData(inMeth, inRegion, outData);


data expand_region;
  set &inRegion.;
  where chr ? "A";
  do VAR3 = region_start to region_stop;
  output; end;
  rename chr=VAR1;
run;

proc sort data=expand_region;
   by VAR1 VAR3;
proc sort data=&inMeth.;
   by VAR1 VAR3;
run;

data sig_region;
 merge &inMeth. (in=in1) expand_region (in=in2);
 by VAR1 VAR3;
 if in1 and in2;
run;


data bin_100_1000;
   set sig_region;
   bin_100 = int(var3/100) + 1;
   bin_1000 = int(var3/1000) + 1;
run;


proc sort data=bin_100_1000;
  by var1 bin_100;
proc means data=bin_100_1000 noprint;
  by var1 bin_100;
  var VAR2 VAR3 VAR4;
  output out=mean_bin_100 min(var2)=var2 max(var3)=var3 mean(VAR4)=var4;
run;

proc means data=bin_100_1000 noprint;
  by var1 bin_1000;
  var var2 var3 VAR4;
  output out=mean_bin_1000 min(var2)=var2 max(var3)=var3 mean(VAR4)=var4;
run;


data mean_bin_100_2;
  set mean_bin_100;
*  bin_start = (bin_100 - 1 ) * 100;
*  bin_stop = bin_100  * 100;
run;

data mean_bin_1000_2;
  set mean_bin_1000;
*  bin_start = (bin_1000 - 1 ) * 1000;
*  bin_stop = bin_1000  * 1000;
run;

data export_100;
*  retain VAR1 bin_start bin_stop VAR4;
*  set mean_bin_100_2;
*  keep VAR1 bin_start bin_stop VAR4;
  retain var1 var2 var3 VAR4;
  set mean_bin_100_2;
  keep var1 var2 var3 VAR4;
run;

data export_1000;
  *retain VAR1 bin_start bin_stop VAR4;
  *set mean_bin_1000_2;
  *keep VAR1 bin_start bin_stop VAR4;
  retain var1 var2 var3 VAR4;
  set mean_bin_1000_2;
  keep var1 var2 var3 VAR4;
run;

proc export data=export_100
     outfile="/TB14/TB14/sandbox/dtra_sandbox/&outData._100bp_window_sigOnly.txt"
     dbms=tab replace; putnames=no;
run;

proc export data=export_1000
     outfile="/TB14/TB14/sandbox/dtra_sandbox/&outData._1kb_window_sigOnly.txt"
     dbms=tab replace; putnames=no;
run;

%mend;

%windowData(cg_10cgy_1h, cg_dmr_01_1, Br_CG_Meth_10cGy_1h);
%windowData(chg_10cgy_1h, chg_dmr_01_1, Br_CHG_Meth_10cGy_1h);
%windowData(chh_10cgy_1h, chh_dmr_01_1, Br_CHH_Meth_10cGy_1h);
%windowData(gc_10cgy_1h, dar_01_1, Br_GC_Acc_10cGy_1h);

%windowData(cg_1p4cgy_1h, cg_dmr_1p4_1, Br_CG_Meth_1p4cGy_1h);
%windowData(chg_1p4cgy_1h, chg_dmr_1p4_1, Br_CHG_Meth_1p4cGy_1h);
%windowData(chh_1p4cgy_1h, chh_dmr_1p4_1, Br_CHH_Meth_1p4cGy_1h);
%windowData(gc_1p4cgy_1h, dar_1p4_1, Br_GC_Acc_1p4cGy_1h);

%windowData(cg_10cgy_72h, cg_dmr_01_72, Br_CG_Meth_10cGy_72h);
%windowData(chg_10cgy_72h, chg_dmr_01_72, Br_CHG_Meth_10cGy_72h);
%windowData(chh_10cgy_72h, chh_dmr_01_72, Br_CHH_Meth_10cGy_72h);
%windowData(gc_10cgy_72h, dar_01_72, Br_GC_Acc_10cGy_72h);

%windowData(cg_1p4cgy_72h, cg_dmr_1p4_72, Br_CG_Meth_1p4cGy_72h);
%windowData(chg_1p4cgy_72h, chg_dmr_1p4_72, Br_CHG_Meth_1p4cGy_72h);
%windowData(chh_1p4cgy_72h, chh_dmr_1p4_72, Br_CHH_Meth_1p4cGy_72h);
%windowData(gc_1p4cgy_72h, dar_1p4_72, Br_GC_Acc_1p4cGy_72h);




/* Sig SITES not regions (with a 100bp pad added) */





/* Group DMCs into windows
 Sites must be no more than 100bp away from the next site
 Sites must have the same direction to be included
 Region must have at least 3 sites
 A DMR has at least 3 DMCs (flag 1 and 3)
 */

data dar_01_1 dar_1p4_1 dar_01_72 dar_1p4_72;
retain chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
set brassMAP.results_by_dac_w_meth;
where chr ? "A";
if flag_meth_diff_20perc ne 1 then delete;
if flag_FDR05_CTL_or_TRT ne 1 then delete;
start_pos2=start_pos-50;
stop_pos2=stop_pos+50;
     if comparison="10cGy_0cGy_1h" then output dar_01_1;
     if comparison="1p4cGy_0cGy_1h" then output dar_1p4_1;
     if comparison="10cGy_0cGy_72h" then output dar_01_72;
     if comparison="1p4cGy_0cGy_72h" then output dar_1p4_72;
keep chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
run;




data cg_dmr_01_1 cg_dmr_1p4_1 cg_dmr_01_72 cg_dmr_1p4_72
     chg_dmr_01_1 chg_dmr_1p4_1 chg_dmr_01_72 chg_dmr_1p4_72
     chh_dmr_01_1 chh_dmr_1p4_1 chh_dmr_01_72 chh_dmr_1p4_72;
retain chr start_pos2 stop_pos2 methyl_diff ;
set brassMAP.results_by_dmc_w_meth;
*where chr ? "A";
start_pos2=start_pos-50;
stop_pos2=stop_pos+50;
if site_type="CG" and comparison="0cGy_vs_10cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output cg_dmr_01_1;
if site_type="CG" and comparison="0cGy_vs_1p4cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output cg_dmr_1p4_1;
if site_type="CG" and comparison="0cGy_vs_10cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output cg_dmr_01_72;
if site_type="CG" and comparison="0cGy_vs_1p4cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output cg_dmr_1p4_72;

if site_type="CHG" and comparison="0cGy_vs_10cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chg_dmr_01_1;
if site_type="CHG" and comparison="0cGy_vs_1p4cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chg_dmr_1p4_1;
if site_type="CHG" and comparison="0cGy_vs_10cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chg_dmr_01_72;
if site_type="CHG" and comparison="0cGy_vs_1p4cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chg_dmr_1p4_72;

if site_type="CHH" and comparison="0cGy_vs_10cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chh_dmr_01_1;
if site_type="CHH" and comparison="0cGy_vs_1p4cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chh_dmr_1p4_1;
if site_type="CHH" and comparison="0cGy_vs_10cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chh_dmr_01_72;
if site_type="CHH" and comparison="0cGy_vs_1p4cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chh_dmr_1p4_72;



keep chr start_pos2 stop_pos2 methyl_diff ;
run;


proc export data=dar_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_10cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=dar_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_1p4cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=dar_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_10cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=dar_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_1p4cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;

proc export data=cg_dmr_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_10cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmr_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_1p4cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmr_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_10cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmr_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_1p4cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;


proc export data=chg_dmr_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_10cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmr_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_1p4cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmr_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_10cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmr_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_1p4cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;

proc export data=chh_dmr_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_10cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmr_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_1p4cGy_1h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmr_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_10cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmr_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_1p4cGy_72h_sig_sites.txt" dbms=tab replace; putnames=no;run;


data dar_01_1 dar_1p4_1 dar_01_72 dar_1p4_72;
retain chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
set brassMAP.results_by_dac_w_meth;
where chr ? "A";
if flag_meth_diff_20perc ne 1 then delete;
if flag_FDR05_CTL_or_TRT ne 1 then delete;
start_pos2=start_pos-500;
stop_pos2=stop_pos+500;
     if comparison="10cGy_0cGy_1h" then output dar_01_1;
     if comparison="1p4cGy_0cGy_1h" then output dar_1p4_1;
     if comparison="10cGy_0cGy_72h" then output dar_01_72;
     if comparison="1p4cGy_0cGy_72h" then output dar_1p4_72;
keep chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
run;









data cg_dmc01 cg_dmc1 chg_dmc01 chg_dmc1 chh_dmc01 chh_dmc1;
retain chr start_pos2 stop_pos2 methyl_diff ;
set brassMAP.results_by_dmc_annot;
where chr ? "A";
start_pos2=start_pos-500;
stop_pos2=stop_pos+500;
if site_type="CG" and comparison="0cGy_vs_10cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output cg_dmr_01_1;
if site_type="CG" and comparison="0cGy_vs_1p4cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output cg_dmr_1p4_1;
if site_type="CG" and comparison="0cGy_vs_10cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output cg_dmr_01_72;
if site_type="CG" and comparison="0cGy_vs_1p4cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output cg_dmr_1p4_72;

if site_type="CHG" and comparison="0cGy_vs_10cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chg_dmr_01_1;
if site_type="CHG" and comparison="0cGy_vs_1p4cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chg_dmr_1p4_1;
if site_type="CHG" and comparison="0cGy_vs_10cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chg_dmr_01_72;
if site_type="CHG" and comparison="0cGy_vs_1p4cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chg_dmr_1p4_72;

if site_type="CHH" and comparison="0cGy_vs_10cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chh_dmr_01_1;
if site_type="CHH" and comparison="0cGy_vs_1p4cGy_1h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chh_dmr_1p4_1;
if site_type="CHH" and comparison="0cGy_vs_10cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chh_dmr_01_72;
if site_type="CHH" and comparison="0cGy_vs_1p4cGy_72h" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chh_dmr_1p4_72;
keep chr start_pos2 stop_pos2 methyl_diff ;
run;



proc export data=dar_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_10cGy_1h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=dar_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_1p4cGy_1h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=dar_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_10cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=dar_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_DAC_1p4cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;

proc export data=cg_dmr_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_10cGy_1h_sig_site2s.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmr_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_1p4cGy_1h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmr_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_10cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmr_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CG_DMC_1p4cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;


proc export data=chg_dmr_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_10cGy_1h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmr_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_1p4cGy_1h_sig_site2s.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmr_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_10cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmr_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHG_DMC_1p4cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;

proc export data=chh_dmr_01_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_10cGy_1h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmr_1p4_1 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_1p4cGy_1h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmr_01_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_10cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmr_1p4_72 outfile="/TB14/TB14/sandbox/dtra_sandbox/Br_CHH_DMC_1p4cGy_72h_sig_sites2.txt" dbms=tab replace; putnames=no;run;




