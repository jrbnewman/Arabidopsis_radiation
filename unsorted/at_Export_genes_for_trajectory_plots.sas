libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

data gene_Data;
  set arabRNA.cpm_norm_counts_by_gene;
  *where gene_id in ("AT3G01830","AT5G13190","AT1G79915","AT3G21330",
            "AT1G78990","AT3G21330","ATMG01370","AT2G30360",
            "AT5G44585","AT5G49690","AT1G72030","AT2G41100",
            "AT5G23510","AT4G34410","AT1G43160","AT3G16770","AT1G63240","AT2G36490",
            "AT3G23240","AT2G01760","AT4G35610","AT3G53600","AT1G61380","AT4G19820");
	    where gene_id in ("AT3G23240","AT5G47220","AT3G50260");
  length group $10.;
  length timepoint $4.;
  if treatment = "Mock" then group="0_Mock";
  else if treatment = "0.1gy" then group="1_10cGy";
  else if treatment = "1gy" then group="2_100cGy";

  if time=1 then timepoint = "01h";
  if time=3 then timepoint = "03h";
  if time=24 then timepoint = "24h";
  if time=72 then timepoint = "72h";
  log_cpm = log(cpm+1);
  keep gene_id replicate cpm group timepoint log_cpm;
run;


/* Export counts */

proc sort data=gene_Data;
  by gene_id group timepoint replicate;
run;

/*
data export_AT3G01830;  set gene_data;   where gene_id="AT3G01830";   drop gene_id;run;
data export_ATGLIP;  set gene_data;   where gene_id="AT5G13190";   drop gene_id;run;
data export_AT1G79915;  set gene_data;   where gene_id="AT1G79915";   drop gene_id;run;
data export_AT3G21330;  set gene_data;   where gene_id="AT3G21330";   drop gene_id;run;
data export_AT1G78990;  set gene_data;   where gene_id="AT1G78990";   drop gene_id;run;
data export_AT3G21330;  set gene_data;   where gene_id="AT3G21330";   drop gene_id;run;
data export_ORF111D;  set gene_data;   where gene_id="ATMG01370";   drop gene_id;run;
data export_SIP4;  set gene_data;   where gene_id="AT2G30360";   drop gene_id;run;
data export_PROSCOOP12;  set gene_data;   where gene_id="AT5G44585";   drop gene_id;run;
data export_UGT91C1;  set gene_data;   where gene_id="AT5G49690";   drop gene_id;run;
data export_GNAT10;  set gene_data;   where gene_id="AT1G72030";   drop gene_id;run;
data export_TCH3;  set gene_data;   where gene_id="AT2G41100";   drop gene_id;run;
data export_AT5G23510;  set gene_data;   where gene_id="AT5G23510";   drop gene_id;run;
data export_RRTF1;  set gene_data;   where gene_id="AT4G34410";   drop gene_id;run;
data export_RAP2_6;  set gene_data;   where gene_id="AT1G43160";   drop gene_id;run;
data export_RAP2_3;  set gene_data;   where gene_id="AT3G16770";   drop gene_id;run;
data export_ERF1;  set gene_data;   where gene_id="AT3G23240";   drop gene_id;run;
data export_ARR14;  set gene_data;   where gene_id="AT2G01760";   drop gene_id;run;
data export_At4g35610;  set gene_data;   where gene_id="AT4G35610";   drop gene_id;run;


data export_ZAT18;  set gene_data;   where gene_id="AT3G53600";   drop gene_id;run;

data export_ZAT_AT1G61380; set gene_data; where gene_id="AT1G61380"; drop gene_id; run;
data export_ZAT_AT4G19820; set gene_data; where gene_id="AT4G19820"; drop gene_id; run;


data export_RMB1_AT1G63240; set gene_data; where gene_id="AT1G63240"; drop gene_id; run;
data export_ROS1_AT2G36490;; set gene_data; where gene_id="AT2G36490"; drop gene_id; run;
*/



data export; set gene_data; where gene_id="AT3G23240"; drop gene_id; run;
proc export data=export outfile="$HOME/concannon/DTRA/new_ERF1_AT3G23240_gene_CPM_counts.csv" dbms=csv replace; run;


data export; set gene_data; where gene_id="AT5G47220"; drop gene_id; run;
proc export data=export outfile="$HOME/concannon/DTRA/new_ERF2_AT5G47220_gene_CPM_counts.csv" dbms=csv replace; run;

data export; set gene_data; where gene_id="AT3G50260"; drop gene_id; run;
proc export data=export outfile="$HOME/concannon/DTRA/new_CEJ1_AT3G50260_gene_CPM_counts.csv" dbms=csv replace; run;
 


/*


proc export data=export_RMB1_AT1G63240 outfile="$HOME/concannon/DTRA/new_RMB1_AT1G63240_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data=export_ROS1_AT2G36490 outfile="$HOME/concannon/DTRA/new_ROS1_AT2G36490_gene_CPM_counts.csv" dbms=csv replace; run;


proc export data=export_ZAT_AT1G61380 outfile="$HOME/concannon/DTRA/new_ZAT_AT1G61380_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data=export_ZAT_AT4G19820 outfile="$HOME/concannon/DTRA/new_ZAT_AT4G19820_gene_CPM_counts.csv" dbms=csv replace; run;

proc export data=export_ZAT18 outfile="$HOME/concannon/DTRA/new_ZAT18_gene_CPM_counts.csv" dbms=csv replace; run;


proc export data=export_AT3G01830 outfile="$HOME/concannon/DTRA/new_AT3G01830_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data=export_ATGLIP outfile="$HOME/concannon/DTRA/new_ATGLIP_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data=export_AT1G79915 outfile="$HOME/concannon/DTRA/new_AT1G79915_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data=export_AT3G21330 outfile="$HOME/concannon/DTRA/new_AT3G21330_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data=export_AT1G78990 outfile="$HOME/concannon/DTRA/new_AT1G78990_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_AT1G78990 outfile="$HOME/concannon/DTRA/new_AT1G78990_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_ORF111D outfile="$HOME/concannon/DTRA/new_ORF111D_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_SIP4 outfile="$HOME/concannon/DTRA/new_SIP4_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_PROSCOOP12 outfile="$HOME/concannon/DTRA/new_PROSCOOP12_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_UGT91C1 outfile="$HOME/concannon/DTRA/new_UGT91C1_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_GNAT10 outfile="$HOME/concannon/DTRA/new_GNAT10_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_TCH3 outfile="$HOME/concannon/DTRA/new_TCH3_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_AT5G23510 outfile="$HOME/concannon/DTRA/new_AT5G23510_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_RRTF1 outfile="$HOME/concannon/DTRA/new_RRTF1_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_RAP2_6 outfile="$HOME/concannon/DTRA/new_RAP2_6_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_RAP2_3 outfile="$HOME/concannon/DTRA/new_RAP2_3_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_ERF1 outfile="$HOME/concannon/DTRA/new_ERF1_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_ARR14 outfile="$HOME/concannon/DTRA/new_ARR14_gene_CPM_counts.csv" dbms=csv replace; run;
proc export data= export_At4g35610 outfile="$HOME/concannon/DTRA/new_At4g35610_gene_CPM_counts.csv" dbms=csv replace; run;

	
*/	
	
	
	



