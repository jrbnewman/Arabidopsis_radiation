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

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/M-1.bedGraph.2"
   out=rna_0cgy_1h dbms=tab replace; getnames=no; 
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/M-3.bedGraph.2"
   out=rna_0cgy_3h dbms=tab replace; getnames=no; 
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/M-24.bedGraph.2"
   out=rna_0cgy_24h dbms=tab replace; getnames=no; 
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/M-72.bedGraph.2"
   out=rna_0cgy_72h dbms=tab replace; getnames=no; 
run;


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/0-1-1.bedGraph.2"
   out=rna_10cgy_1h dbms=tab replace; getnames=no; 
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/0-1-3.bedGraph.2"
   out=rna_10cgy_3h dbms=tab replace; getnames=no; 
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/0-1-24.bedGraph.2"
   out=rna_10cgy_24h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/0-1-72.bedGraph.2"
   out=rna_10cgy_72h dbms=tab replace; getnames=no;
run;


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/1-1.bedGraph.2"
   out=rna_100cgy_1h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/1-3.bedGraph.2"
   out=rna_100cgy_3h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/1-24.bedGraph.2"
   out=rna_100cgy_24h dbms=tab replace; getnames=no;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/1-72.bedGraph.2"
   out=rna_100cgy_72h dbms=tab replace; getnames=no;
run;

/* set up diffs for RNA */

data rna_0cgy_1h_2;  set rna_0cgy_1h;  rename VAR4=mock; run;
data rna_0cgy_3h_2;  set rna_0cgy_3h;  rename VAR4=mock; run;
data rna_0cgy_24h_2;  set rna_0cgy_24h;  rename VAR4=mock; run;
data rna_0cgy_72h_2;  set rna_0cgy_72h;  rename VAR4=mock; run;

data rna_10cgy_1h_2;  set rna_10cgy_1h;  rename VAR4=_10cgy; run;
data rna_10cgy_3h_2;  set rna_10cgy_3h;  rename VAR4=_10cgy; run;
data rna_10cgy_24h_2;  set rna_10cgy_24h;  rename VAR4=_10cgy; run;
data rna_10cgy_72h_2;  set rna_10cgy_72h;  rename VAR4=_10cgy; run;

data rna_100cgy_1h_2;  set rna_100cgy_1h;  rename VAR4=_100cgy; run;
data rna_100cgy_3h_2;  set rna_100cgy_3h;  rename VAR4=_100cgy; run;
data rna_100cgy_24h_2;  set rna_100cgy_24h;  rename VAR4=_100cgy; run;
data rna_100cgy_72h_2;  set rna_100cgy_72h;  rename VAR4=_100cgy; run;

proc sort data=rna_0cgy_1h_2; by var1 var2 var3;
proc sort data=rna_10cgy_1h_2; by var1 var2 var3;
proc sort data=rna_100cgy_1h_2; by var1 var2 var3;

proc sort data=rna_0cgy_3h_2; by var1 var2 var3;
proc sort data=rna_10cgy_3h_2; by var1 var2 var3;
proc sort data=rna_100cgy_3h_2; by var1 var2 var3;

proc sort data=rna_0cgy_24h_2; by var1 var2 var3;
proc sort data=rna_10cgy_24h_2; by var1 var2 var3;
proc sort data=rna_100cgy_24h_2; by var1 var2 var3;

proc sort data=rna_0cgy_72h_2; by var1 var2 var3;
proc sort data=rna_10cgy_72h_2; by var1 var2 var3;
proc sort data=rna_100cgy_72h_2; by var1 var2 var3;
run;

data rna_10cgy_Mock_1h;
     merge rna_10cgy_1h_2 (in=in1) rna_0cgy_1h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _10cgy=0;
     if not in2 then mock=0;
run;

data rna_10cgy_Mock_3h;
     merge rna_10cgy_3h_2 (in=in1) rna_0cgy_3h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _10cgy=0;
     if not in2 then mock=0;
run;

data rna_10cgy_Mock_24h;
     merge rna_10cgy_24h_2 (in=in1) rna_0cgy_24h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _10cgy=0;
     if not in2 then mock=0;
run;

data rna_10cgy_Mock_72h;
     merge rna_10cgy_72h_2 (in=in1) rna_0cgy_72h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _10cgy=0;
     if not in2 then mock=0;
run;


data rna_100cgy_Mock_1h;
     merge rna_100cgy_1h_2 (in=in1) rna_0cgy_1h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _100cgy=0;
     if not in2 then mock=0;
run;

data rna_100cgy_Mock_3h;
     merge rna_100cgy_3h_2 (in=in1) rna_0cgy_3h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _100cgy=0;
     if not in2 then mock=0;
run;

data rna_100cgy_Mock_24h;
     merge rna_100cgy_24h_2 (in=in1) rna_0cgy_24h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _100cgy=0;
     if not in2 then mock=0;
run;

data rna_100cgy_Mock_72h;
     merge rna_100cgy_72h_2 (in=in1) rna_0cgy_72h_2 (in=in2) ;
     by var1 var2 var3;
     if not in1 then _100cgy=0;
     if not in2 then mock=0;
run;


data rna_10cgy_mock_1h_2; set rna_10cgy_mock_1h; var4=log2(_10cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;
data rna_10cgy_mock_3h_2; set rna_10cgy_mock_3h; var4=log2(_10cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;
data rna_10cgy_mock_24h_2; set rna_10cgy_mock_24h; var4=log2(_10cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;
data rna_10cgy_mock_72h_2; set rna_10cgy_mock_72h; var4=log2(_10cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;

data rna_100cgy_mock_1h_2; set rna_100cgy_mock_1h; var4=log2(_100cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;
data rna_100cgy_mock_3h_2; set rna_100cgy_mock_3h; var4=log2(_100cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;
data rna_100cgy_mock_24h_2; set rna_100cgy_mock_24h; var4=log2(_100cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;
data rna_100cgy_mock_72h_2; set rna_100cgy_mock_72h; var4=log2(_100cgy+1)-log2(Mock+1); keep var1 var2 var3 var4; run;


/* Macro to reduce data into bins of 100bp or 1000bp */

%macro windowData(inData, outData);

data bin_100_1000;
   set &inData.;
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



%windowData(rna_10cgy_mock_1h_2, At_RNA_10cGy_1h);
%windowData(rna_10cgy_mock_3h_2,  At_RNA_10cGy_3h);
%windowData(rna_10cgy_mock_24h_2,  At_RNA_10cGy_24h);
%windowData(rna_10cgy_mock_72h_2,  At_RNA_10cGy_72h);
%windowData(rna_100cgy_mock_1h_2,  At_RNA_100cGy_1h);
%windowData(rna_100cgy_mock_3h_2, At_RNA_100cGy_3h);
%windowData(rna_100cgy_mock_24h_2, At_RNA_100cGy_24h);
%windowData(rna_100cgy_mock_72h_2, At_RNA_100cGy_72h);
%windowData(cg_10cgy, At_CG_Meth_10cGy);
%windowData(cg_100cgy, At_CG_Meth_100cGy);
%windowData(chg_10cgy, At_CHG_Meth_10cGy);
%windowData(chg_100cgy, At_CHG_Meth_100cGy);
%windowData(chh_10cgy, At_CHH_Meth_10cGy);
%windowData(chh_100cgy, At_CHH_Meth_100cGy);
%windowData(gc_10cgy, At_GC_Acc_10cGy);
%windowData(gc_100cgy, At_GC_Acc_100cGy);



/* Make tile tracks for DEG, DMR, DAR */



data de_results;
  set arabRNA.arab_results_by_gene;
  keep gene_id flag_Mock_1hr_on flag_Mock_3hr_on flag_Mock_24hr_on flag_Mock_72hr_on
  flag_01gy_1hr_on flag_01gy_3hr_on flag_01gy_24hr_on flag_01gy_72hr_on
  flag_1gy_1hr_on flag_1gy_3hr_on flag_1gy_24hr_on flag_1gy_72hr_on
  mean_cpm_: 
  fdr_: 
  flag_01gy_v_Mock_1h_fdr05 flag_01gy_v_Mock_3h_fdr05
  flag_01gy_v_Mock_24h_fdr05 flag_01gy_v_Mock_72h_fdr05
  flag_1gy_v_Mock_1h_fdr05 flag_1gy_v_Mock_3h_fdr05
  flag_1gy_v_Mock_24h_fdr05 flag_1gy_v_Mock_72h_fdr05 ;
run;

data de_results2;
  set de_results;
  log2fc_01gy_Mock_1h = log2(mean_cpm_01gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_01gy_Mock_3h = log2(mean_cpm_01gy_3h ) - log2(mean_cpm_Mock_3h );
  log2fc_01gy_Mock_24h = log2(mean_cpm_01gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_01gy_Mock_72h = log2(mean_cpm_01gy_72h ) - log2(mean_cpm_Mock_72h );
  log2fc_1gy_Mock_1h = log2(mean_cpm_1gy_1h ) - log2(mean_cpm_Mock_1h );
  log2fc_1gy_Mock_3h = log2(mean_cpm_1gy_3h) - log2(mean_cpm_Mock_3h );
  log2fc_1gy_Mock_24h = log2(mean_cpm_1gy_24h ) - log2(mean_cpm_Mock_24h );
  log2fc_1gy_Mock_72h = log2(mean_cpm_1gy_72h ) - log2(mean_cpm_Mock_72h );
run;


data  up_01_1_fc1 up_01_3_fc1 up_01_24_fc1 up_01_72_fc1
      dn_01_1_fc1 dn_01_3_fc1 dn_01_24_fc1 dn_01_72_fc1
      no_01_1_fc1 no_01_3_fc1 no_01_24_fc1 no_01_72_fc1
      up_1_1_fc1 up_1_3_fc1 up_1_24_fc1 up_1_72_fc1
      dn_1_1_fc1 dn_1_3_fc1 dn_1_24_fc1 dn_1_72_fc1
      no_1_1_fc1 no_1_3_fc1 no_1_24_fc1 no_1_72_fc1
      off_01_1_fc1 off_01_3_fc1 off_01_24_fc1 off_01_72_fc1
      off_1_1_fc1 off_1_3_fc1 off_1_24_fc1 off_1_72_fc1;
     set de_Results2;
if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1_fc1;
else if ((flag_01gy_v_mock_1h_fdr05=1 and abs(log2fc_01gy_Mock_1h) >= 1) or (flag_01gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_01gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1_fc1;
else if (mean_cpm_01gy_1h > 0 and mean_cpm_Mock_1h > 0) then output no_01_1_fc1;
else output off_01_1_fc1;

if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3_fc1;
else if ((flag_01gy_v_mock_3h_fdr05=1 and abs(log2fc_01gy_Mock_3h) >= 1) or (flag_01gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_01gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3_fc1;
else if (mean_cpm_01gy_3h > 0 and mean_cpm_Mock_3h > 0) then output no_01_3_fc1;
else output off_01_3_fc1;

if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24_fc1;
else if ((flag_01gy_v_mock_24h_fdr05=1 and abs(log2fc_01gy_Mock_24h) >= 1) or (flag_01gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_01gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24_fc1;
else if (mean_cpm_01gy_24h > 0 and mean_cpm_Mock_24h > 0) then output no_01_24_fc1;
else output off_01_24_fc1;

if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72_fc1;
else if ((flag_01gy_v_mock_72h_fdr05=1 and abs(log2fc_01gy_Mock_72h) >= 1) or (flag_01gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_01gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72_fc1;
else if (mean_cpm_01gy_72h > 0 and mean_cpm_Mock_72h > 0) then output no_01_72_fc1;
else output off_01_72_fc1;

if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1_fc1;
else if ((flag_1gy_v_mock_1h_fdr05=1 and abs(log2fc_1gy_Mock_1h) >= 1) or (flag_1gy_1hr_on ^= flag_Mock_1hr_on and log2(abs(mean_cpm_1gy_1h - mean_cpm_Mock_1h)) >= 1 ))  and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1_fc1;
else if (mean_cpm_1gy_1h > 0 and mean_cpm_Mock_1h > 0) then output no_1_1_fc1;
else output off_1_1_fc1;

if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3_fc1;
else if ((flag_1gy_v_mock_3h_fdr05=1 and abs(log2fc_1gy_Mock_3h) >= 1) or (flag_1gy_3hr_on ^= flag_Mock_3hr_on and log2(abs(mean_cpm_1gy_3h - mean_cpm_Mock_3h)) >= 1 ))  and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3_fc1;
else if (mean_cpm_1gy_3h > 0 and mean_cpm_Mock_3h > 0) then output no_1_3_fc1;
else output off_1_3_fc1;

if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24_fc1;
else if ((flag_1gy_v_mock_24h_fdr05=1 and abs(log2fc_1gy_Mock_24h) >= 1) or (flag_1gy_24hr_on ^= flag_Mock_24hr_on and log2(abs(mean_cpm_1gy_24h - mean_cpm_Mock_24h)) >= 1 ))  and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24_fc1;
else if (mean_cpm_1gy_24h > 0 and mean_cpm_Mock_24h > 0) then output no_1_24_fc1;
else output off_1_24_fc1;

if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72_fc1;
else if ((flag_1gy_v_mock_72h_fdr05=1 and abs(log2fc_1gy_Mock_72h) >= 1) or (flag_1gy_72hr_on ^= flag_Mock_72hr_on and log2(abs(mean_cpm_1gy_72h - mean_cpm_Mock_72h)) >= 1 ))  and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72_fc1;
else if (mean_cpm_1gy_72h > 0 and mean_cpm_Mock_72h > 0) then output no_1_72_fc1;
else output off_1_72_fc1;

keep gene_id;
run;

/* Get DMRs/DARs by gene  and feature annotation */



data up_dar_01_72_gn up_dar_1_72_gn dn_dar_01_72_gn dn_dar_1_72_gn;
     set arabMAP.results_by_dar_annot;
     if comparison="01Gy_0Gy" then comparison2="0Gy_vs_01G";
     if comparison="1Gy_0Gy" then comparison2="0Gy_vs_1Gy";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_01G" then output up_dar_01_72_gn;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_01G" then output dn_dar_01_72_gn;

     if mean_methyl_diff_TRT_CTL > 0 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72_gn;
     if mean_methyl_diff_TRT_CTL < 0 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72_gn;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_01G" then output up_dar_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_01G" then output dn_dar_01_72_gn;

     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison2="0Gy_vs_1Gy" then output up_dar_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison2="0Gy_vs_1Gy" then output dn_dar_1_72_gn;

     *keep chr region_start region_stop;
run;





data cg_up_dmr_01_72_gn cg_up_dmr_1_72_gn cg_dn_dmr_01_72_gn cg_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CG";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output cg_up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output cg_up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output cg_dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output cg_dn_dmr_1_72_gn;
     keep chr region_start region_stop;
run;



data chg_up_dmr_01_72_gn chg_up_dmr_1_72_gn chg_dn_dmr_01_72_gn chg_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHG";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chg_up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chg_up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chg_dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chg_dn_dmr_1_72_gn;
     keep chr region_start region_stop;
run;


data chh_up_dmr_01_72_gn chh_up_dmr_1_72_gn chh_dn_dmr_01_72_gn chh_dn_dmr_1_72_gn;
     set arabMAP.results_by_dmr_annot;
     where site_type="CHH";
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_gn;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_gn;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_gn;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_gn;
     keep chr region_start region_stop;
run;


proc sort data=up_dar_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=up_dar_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=dn_dar_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=dn_dar_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=cg_up_dmr_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=cg_up_dmr_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=cg_dn_dmr_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=cg_dn_dmr_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chg_up_dmr_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chg_up_dmr_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chg_dn_dmr_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chg_dn_dmr_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chh_up_dmr_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chh_up_dmr_1_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chh_dn_dmr_01_72_gn nodup; by chr region_start region_stop; run;
proc sort data=chh_dn_dmr_1_72_gn nodup; by chr region_start region_stop; run;

proc sort data=up_01_72_fc1 nodup; by gene_id ; run;
proc sort data=up_01_24_fc1 nodup; by gene_id ; run;
proc sort data=up_01_3_fc1 nodup; by gene_id ; run;
proc sort data=up_01_1_fc1 nodup; by gene_id ; run;

proc sort data=dn_01_72_fc1 nodup; by gene_id ; run;
proc sort data=dn_01_24_fc1 nodup; by gene_id ; run;
proc sort data=dn_01_3_fc1 nodup; by gene_id ; run;
proc sort data=dn_01_1_fc1 nodup; by gene_id ; run;

proc sort data=up_1_72_fc1 nodup; by gene_id ; run;
proc sort data=up_1_24_fc1 nodup; by gene_id ; run;
proc sort data=up_1_3_fc1 nodup; by gene_id ; run;
proc sort data=up_1_1_fc1 nodup; by gene_id ; run;

proc sort data=dn_1_72_fc1 nodup; by gene_id ; run;
proc sort data=dn_1_24_fc1 nodup; by gene_id ; run;
proc sort data=dn_1_3_fc1 nodup; by gene_id ; run;
proc sort data=dn_1_1_fc1 nodup; by gene_id ; run;


data tiles_10cgy_rna_1h;
length options $100.;
set up_01_1_fc1 (in=in1) dn_01_1_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_10cgy_rna_3h;
length options $100.;
set up_01_3_fc1 (in=in1) dn_01_3_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_10cgy_rna_24h;
length options $100.;
set up_01_24_fc1 (in=in1) dn_01_24_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_10cgy_rna_72h;
length options $100.;
set up_01_72_fc1 (in=in1) dn_01_72_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;



data tiles_100cgy_rna_1h;
length options $100.;
set up_1_1_fc1 (in=in1) dn_1_1_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_100cgy_rna_3h;
length options $100.;
set up_1_3_fc1 (in=in1) dn_1_3_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_100cgy_rna_24h;
length options $100.;
set up_1_24_fc1 (in=in1) dn_1_24_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_100cgy_rna_72h;
length options $100.;
set up_1_72_fc1 (in=in1) dn_1_72_fc1 (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;


data tiles_10cgy_cg;
length options $100.;
set cg_up_dmr_01_72_gn (in=in1) cg_dn_dmr_01_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_10cgy_chg;
length options $100.;
set chg_up_dmr_01_72_gn (in=in1) chg_dn_dmr_01_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_10cgy_chh;
length options $100.;
set chh_up_dmr_01_72_gn (in=in1) chh_dn_dmr_01_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;

data tiles_10cgy_gc;
length options $100.;
set up_dar_01_72_gn (in=in1) dn_dar_01_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;




data tiles_100cgy_cg;
length options $100.;
set cg_up_dmr_1_72_gn (in=in1) cg_dn_dmr_1_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";run;

data tiles_100cgy_chg;
length options $100.;
set chg_up_dmr_1_72_gn (in=in1) chg_dn_dmr_1_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";run;

data tiles_100cgy_chh;
length options $100.;
set chh_up_dmr_1_72_gn (in=in1) chh_dn_dmr_1_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";run;

data tiles_100cgy_gc;
length options $100.;
set up_dar_1_72_gn (in=in1) dn_dar_1_72_gn (in=in2);
*if in1 then options="color=orange";
*if in2 then options="color=blue";
if in1 then options="color=black";
if in2 then options="color=black";
run;


/* Get gene coordinates */

libname tair '!HOME/concannon/useful_arabidopsis_data/tair10/sas_data';

data exon_coord;
  set tair.tair20_exons_w_info;
  keep chrom start stop gene_id;
run;

proc sort data=exon_coord;
  by chrom gene_id start stop;
proc means data=exon_coord noprint;
  by chrom gene_id;
  var start stop ;
  output out=gene_coord(drop=_TYPE_ _FREQ_) min(start)=gene_start max(stop)=gene_stop;
run;

proc sort data=gene_coord;  by gene_id;
proc sort data=tiles_10cgy_rna_1h; by gene_id;
proc sort data=tiles_10cgy_rna_3h; by gene_id;
proc sort data=tiles_10cgy_rna_24h; by gene_id;
proc sort data=tiles_10cgy_rna_72h; by gene_id;
proc sort data=tiles_100cgy_rna_1h; by gene_id;
proc sort data=tiles_100cgy_rna_3h; by gene_id;
proc sort data=tiles_100cgy_rna_24h; by gene_id;
proc sort data=tiles_100cgy_rna_72h; by gene_id;
run;

data tiles_10cgy_rna_1h_coord;
  merge tiles_10cgy_rna_1h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;

data tiles_10cgy_rna_3h_coord;
  merge tiles_10cgy_rna_3h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;

data tiles_10cgy_rna_24h_coord;
  merge tiles_10cgy_rna_24h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;

data tiles_10cgy_rna_72h_coord;
  merge tiles_10cgy_rna_72h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;



data tiles_100cgy_rna_1h_coord;
  merge tiles_100cgy_rna_1h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;

data tiles_100cgy_rna_3h_coord;
  merge tiles_100cgy_rna_3h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;

data tiles_100cgy_rna_24h_coord;
  merge tiles_100cgy_rna_24h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;

data tiles_100cgy_rna_72h_coord;
  merge tiles_100cgy_rna_72h (in=in1) gene_coord (in=in2);
  by gene_id;
  if in1 and in2;
run;


/* export tiles */

data tiles_10cgy_rna_1h_export;   retain chr2 gene_start gene_stop options;   set tiles_10cgy_rna_1h_coord; length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_10cgy_rna_3h_export;   retain chr2 gene_start gene_stop options;   set tiles_10cgy_rna_3h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_10cgy_rna_24h_export;   retain chr2 gene_start gene_stop options;   set tiles_10cgy_rna_24h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_10cgy_rna_72h_export;   retain chr2 gene_start gene_stop options;   set tiles_10cgy_rna_72h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_10cgy_cg_export; retain chr2 region_start region_stop options; set tiles_10cgy_cg;  length chr2 $6.; chr2=compress(catt("chr",chr)); keep chr2 region_Start region_stop options; run;
data tiles_10cgy_chg_export; retain chr2 region_start region_stop options; set tiles_10cgy_chg; length chr2 $6.; chr2=compress(catt("chr",chr));keep chr2 region_Start region_stop options; run;
data tiles_10cgy_chh_export; retain chr2 region_start region_stop options; set tiles_10cgy_chh; length chr2 $6.; chr2=compress(catt("chr",chr));keep chr2 region_Start region_stop options; run;
data tiles_10cgy_gc_export; retain chr2 region_start region_stop options; set tiles_10cgy_gc; length chr2 $6.; chr2=compress(catt("chr",chr));keep chr2 region_Start region_stop options; run;


data tiles_100cgy_rna_1h_export;   retain chr2 gene_start gene_stop options;   set tiles_100cgy_rna_1h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_100cgy_rna_3h_export;   retain chr2 gene_start gene_stop options;   set tiles_100cgy_rna_3h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_100cgy_rna_24h_export;   retain chr2 gene_start gene_stop options;   set tiles_100cgy_rna_24h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_100cgy_rna_72h_export;   retain chr2 gene_start gene_stop options;   set tiles_100cgy_rna_72h_coord;  length chr2 $6.; chr2=compress(catt("chr",chrom)); keep chr2 gene_start gene_stop options; run;
data tiles_100cgy_cg_export; retain chr2 region_start region_stop options; set tiles_100cgy_cg;length chr2 $6.; chr2=compress(catt("chr",chr)); keep chr2 region_Start region_stop options; run;
data tiles_100cgy_chg_export; retain chr2 region_start region_stop options; set tiles_100cgy_chg;length chr2 $6.; chr2=compress(catt("chr",chr)); keep chr2 region_Start region_stop options; run;
data tiles_100cgy_chh_export; retain chr2 region_start region_stop options; set tiles_100cgy_chh;length chr2 $6.; chr2=compress(catt("chr",chr)); keep chr2 region_Start region_stop options; run;
data tiles_100cgy_gc_export; retain chr2 region_start region_stop options; set tiles_100cgy_gc;length chr2 $6.; chr2=compress(catt("chr",chr));keep chr2 region_Start region_stop options; run;

%macro exportTRNAiles(inData, outFile);

proc sort data=&inData.;
 by chr2 gene_start gene_stop options;
run;


proc export data=&inData. outfile="/TB14/TB14/sandbox/dtra_sandbox/&outFile._2.txt"
     dbms=tab replace;
     putnames=no;
run;

%mend;

%exportTRNAiles(tiles_10cgy_rna_1h_export, At_10cGy_1h_DEG_tiles);
%exportTRNAiles(tiles_10cgy_rna_3h_export, At_10cGy_3h_DEG_tiles);
%exportTRNAiles(tiles_10cgy_rna_24h_export, At_10cGy_24h_DEG_tiles);
%exportTRNAiles(tiles_10cgy_rna_72h_export, At_10cGy_72h_DEG_tiles);

%exportTRNAiles(tiles_100cgy_rna_1h_export, At_100cGy_1h_DEG_tiles);
%exportTRNAiles(tiles_100cgy_rna_3h_export, At_100cGy_3h_DEG_tiles);
%exportTRNAiles(tiles_100cgy_rna_24h_export, At_100cGy_24h_DEG_tiles);
%exportTRNAiles(tiles_100cgy_rna_72h_export, At_100cGy_72h_DEG_tiles);


%macro exportTMAPiles(inData, outFile);

proc sort data=&inData.;
 by chr2 region_start region_stop options;
run;

proc export data=&inData. outfile="/TB14/TB14/sandbox/dtra_sandbox/&outFile._2.txt"
     dbms=tab replace;
     putnames=no;
run;

%mend;
%exportTMAPiles(tiles_10cgy_cg_export, At_10cGy_CG_DMR_tiles);
%exportTMAPiles(tiles_10cgy_chg_export, At_10cGy_CHG_DMR_tiles);
%exportTMAPiles(tiles_10cgy_chh_export, At_10cGy_CHH_DMR_tiles);
%exportTMAPiles(tiles_10cgy_gc_export, At_10cGy_GC_DAR_tiles);

%exportTMAPiles(tiles_100cgy_cg_export, At_100cGy_CG_DMR_tiles);
%exportTMAPiles(tiles_100cgy_chg_export, At_100cGy_CHG_DMR_tiles);
%exportTMAPiles(tiles_100cgy_chh_export, At_100cGy_CHH_DMR_tiles);
%exportTMAPiles(tiles_100cgy_gc_export, At_100cGy_GC_DAR_tiles);






