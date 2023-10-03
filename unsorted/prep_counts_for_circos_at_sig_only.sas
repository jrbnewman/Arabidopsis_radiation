/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

/* Import bedGRaphs and expression TSVs and average in 100bp and 1kb windows */

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_CG_01Gy_0Gy_analyzed.bedGraph"
   out=cg_10cgy dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_CG_1Gy_0Gy_analyzed.bedGraph"
   out=cg_100cgy dbms=tab replace; getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_CHG_01Gy_0Gy_analyzed.bedGraph"
   out=chg_10cgy dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_CHG_1Gy_0Gy_analyzed.bedGraph"
   out=chg_100cgy dbms=tab replace; getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_CHH_01Gy_0Gy_analyzed.bedGraph"
   out=chh_10cgy dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_CHH_1Gy_0Gy_analyzed.bedGraph"
   out=chh_100cgy dbms=tab replace; getnames=no;
run;
proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_GC_01Gy_0Gy_analyzed.bedGraph"
   out=gc_10cgy dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/bedGraphs/Arabidopsis_GC_1Gy_0Gy_analyzed.bedGraph"
   out=gc_100cgy dbms=tab replace; getnames=no;
run;


/* Set up sig sites */

/* Make tile tracks for DEG, DMR, DAR */

data dar_01_72 dar_1_72;
     set arabMAP.results_by_dac_annot;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if comparison2="0Gy_vs_01G" then output dar_01_72;
     if comparison2="0Gy_vs_1Gy" then output dar_1_72;

     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison2="0Gy_vs_01G" then output dar_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff_TRT_CTL) >=0.1 and comparison2="0Gy_vs_1Gy" then output dar_1_72;

     keep chr region_start region_stop;
run;

data cg_dmr_01_72 cg_dmr_1_72 ;
     set arabMAP.results_by_dmr_annot;
     where site_type="CG";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2  and comparison="0Gy_vs_01G" then output cg_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0Gy_vs_1Gy" then output cg_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0Gy_vs_01G" then output cg_dmr_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0Gy_vs_1Gy" then output cg_dmr_1_72;
     keep chr region_start region_stop;
run;

data chg_dmr_01_72 chg_dmr_1_72 ;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHG";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2  and comparison="0Gy_vs_01G" then output chg_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0Gy_vs_1Gy" then output chg_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0Gy_vs_01G" then output chg_dmr_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0Gy_vs_1Gy" then output chg_dmr_1_72;
     keep chr region_start region_stop;
run;

data chh_dmr_01_72 chh_dmr_1_72 ;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHH";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2  and comparison="0Gy_vs_01G" then output chh_dmr_01_72;
     if num_sites_FDR05_diff_10perc>=2  and comparison="0Gy_vs_1Gy" then output chh_dmr_1_72;
     /* binomial */
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0Gy_vs_01G" then output chh_dmr_01_72;
     if flag_binomial_FDR05=1 and abs(mean_methyl_diff) >=0.1 and comparison="0Gy_vs_1Gy" then output chh_dmr_1_72;
     keep chr region_start region_stop;
run;





proc sort data=dar_01_72 nodup; by chr region_start region_stop; run;
proc sort data=dar_1_72 nodup; by chr region_start region_stop; run;
proc sort data=cg_dmr_01_72 nodup; by chr region_start region_stop; run;
proc sort data=cg_dmr_1_72 nodup; by chr region_start region_stop; run;
proc sort data=chg_dmr_01_72 nodup; by chr region_start region_stop; run;
proc sort data=chg_dmr_1_72 nodup; by chr region_start region_stop; run;
proc sort data=chh_dmr_01_72 nodup; by chr region_start region_stop; run;
proc sort data=chh_dmr_1_72 nodup; by chr region_start region_stop; run;





/* Macro to subset only sig regions and also to reduce data into bins of 100bp or 1000bp */

%macro windowData(inMeth, inRegion, outData);


data expand_region;
  set &inRegion.;
  VAR1 = chr + 0;
  do VAR3 = region_start to region_stop;
  output; end;
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
   length var1_2 $4.;
   set sig_region;
   var1_2=catt("chr",var1);
   bin_100 = int(var3/100) + 1;
   bin_1000 = int(var3/1000) + 1;
run;


proc sort data=bin_100_1000;
  by var1_2 bin_100;
proc means data=bin_100_1000 noprint;
  by var1_2 bin_100;
  var VAR2 VAR3 VAR4;
  output out=mean_bin_100 min(var2)=var2 max(var3)=var3 mean(VAR4)=var4;
run;

proc means data=bin_100_1000 noprint;
  by var1_2 bin_1000;
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
  retain var1_2 var2 var3 VAR4;
  set mean_bin_100_2;
  keep var1_2 var2 var3 VAR4;
run;

data export_1000;
  *retain VAR1 bin_start bin_stop VAR4;
  *set mean_bin_1000_2;
  *keep VAR1 bin_start bin_stop VAR4;
  retain var1_2 var2 var3 VAR4;
  set mean_bin_1000_2;
  keep var1_2 var2 var3 VAR4;
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

%windowData(cg_10cgy, cg_dmr_01_72, At_CG_Meth_10cGy);
%windowData(chg_10cgy, chg_dmr_01_72, At_CHG_Meth_10cGy);
%windowData(chh_10cgy, chh_dmr_01_72, At_CHH_Meth_10cGy);
%windowData(gc_10cgy, dar_01_72, At_GC_Acc_10cGy);

%windowData(cg_100cgy, cg_dmr_1_72, At_CG_Meth_100cGy);
%windowData(chg_100cgy, chg_dmr_1_72, At_CHG_Meth_100cGy);
%windowData(chh_100cgy, chh_dmr_1_72, At_CHH_Meth_100cGy);
%windowData(gc_100cgy, dar_1_72, At_GC_Acc_100cGy);


/* Sig SITES not regions (with a 100bp pad added) */





/* Group DMCs into windows
 Sites must be no more than 100bp away from the next site
 Sites must have the same direction to be included
 Region must have at least 3 sites
 A DMR has at least 3 DMCs (flag 1 and 3)
 */

data dac01 dac1;
retain chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
set arabMAP.results_by_dac_w_meth;
if flag_meth_diff_20perc ne 1 then delete;
if flag_FDR05_CTL_or_TRT ne 1 then delete;
start_pos2=start_pos-50;
stop_pos2=stop_pos+50;
if comparison="01Gy_0Gy" then output dac01;
if comparison="1Gy_0Gy" then output dac1;
keep chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
run;



data dmc;
set methyl_data_all2;
where flag_meth_diff_10perc=1;
run;



data cg_dmc01 cg_dmc1 chg_dmc01 chg_dmc1 chh_dmc01 chh_dmc1;
retain chr start_pos2 stop_pos2 methyl_diff ;
set arabMAP.results_by_dmc_w_meth;
start_pos2=start_pos-50;
stop_pos2=stop_pos+50;
if site_type="CG" and comparison="0Gy_vs_01G" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output cg_dmc01;
if site_type="CG" and comparison="0Gy_vs_1Gy" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output cg_dmc1;
if site_type="CHG" and comparison="0Gy_vs_01G" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chg_dmc01;
if site_type="CHG" and comparison="0Gy_vs_1Gy" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chg_dmc1;
if site_type="CHH" and comparison="0Gy_vs_01G" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chh_dmc01;
if site_type="CHH" and comparison="0Gy_vs_1Gy" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chh_dmc1;
keep chr start_pos2 stop_pos2 methyl_diff ;
run;


proc export data=dac01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAC_10cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=dac1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAC_100cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;

proc export data=cg_dmc01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_DMC_10cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmc1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_DMC_100cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;

proc export data=chg_dmc01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_DMC_10cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmc1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_DMC_100cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;

proc export data=chh_dmc01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_DMC_10cGy_sig_sites.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmc1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_DMC_100cGy_sig_sites.txt" dbms=tab replace; putnames=no;
run;


data dac01 dac1;
retain chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
set arabMAP.results_by_dac_w_meth;
if flag_meth_diff_20perc ne 1 then delete;
if flag_FDR05_CTL_or_TRT ne 1 then delete;
start_pos2=start_pos-500;
stop_pos2=stop_pos+500;
if comparison="01Gy_0Gy" then output dac01;
if comparison="1Gy_0Gy" then output dac1;
keep chr start_pos2 stop_pos2 methyl_diff_TRT_CTL ;
run;


data cg_dmc01 cg_dmc1 chg_dmc01 chg_dmc1 chh_dmc01 chh_dmc1;
retain chr start_pos2 stop_pos2 methyl_diff ;
set arabMAP.results_by_dmc_annot;
start_pos2=start_pos-500;
stop_pos2=stop_pos+500;
if site_type="CG" and comparison="0Gy_vs_01G" and FET_FDR_P < 0.05  and flag_meth_diff_10perc=1 then output cg_dmc01;
if site_type="CG" and comparison="0Gy_vs_1Gy" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output cg_dmc1;
if site_type="CHG" and comparison="0Gy_vs_01G" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chg_dmc01;
if site_type="CHG" and comparison="0Gy_vs_1Gy" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chg_dmc1;
if site_type="CHH" and comparison="0Gy_vs_01G" and FET_FDR_P < 0.05 and FET_FDR_P ne . and flag_meth_diff_10perc=1 then output chh_dmc01;
if site_type="CHH" and comparison="0Gy_vs_1Gy" and FET_FDR_P < 0.05 and FET_FDR_P ne .  and flag_meth_diff_10perc=1 then output chh_dmc1;
keep chr start_pos2 stop_pos2 methyl_diff ;
run;


proc export data=dac01 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAC_10cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=dac1 outfile="/TB14/TB14/sandbox/dtra_sandbox/DAC_100cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;

proc export data=cg_dmc01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_DMC_10cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=cg_dmc1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CG_DMC_100cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;

proc export data=chg_dmc01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_DMC_10cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chg_dmc1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHG_DMC_100cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;

proc export data=chh_dmc01 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_DMC_10cGy_sig_sites2.txt" dbms=tab replace; putnames=no;run;
proc export data=chh_dmc1 outfile="/TB14/TB14/sandbox/dtra_sandbox/CHH_DMC_100cGy_sig_sites2.txt" dbms=tab replace; putnames=no;
run;

