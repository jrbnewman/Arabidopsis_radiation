/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


data counts;
  set arabRNA.cpm_norm_counts_by_gene;
  length group $12.;
  length timepoint $4.;
  where gene_id = "AT2G40340" or gene_id = "AT3G15500"
     or gene_id = "AT1G43160" or gene_id = "AT3G23240" or gene_id="AT4G34410";
  log_cpm=log(cpm+1);
  if treatment="Mock" then group="Mock";
  if treatment="0.1gy" then group="10cGy";
  if treatment="1gy" then group="100cGy";
  if time = 1 then timepoint = "01h";
  else if time = 3 then timepoint = "03h";
  else if time = 24 then timepoint = "24h";
  else if time = 72 then timepoint = "72h";
  keep replicate timepoint group gene_id cpm log_cpm;
run;	


proc sort data=counts;
  by gene_id timepoint group replicate;
run;


data gene1 gene2 gene3 gene4 gene5;
  set counts;
  if gene_id = "AT2G40340" then output gene1;
  if gene_id = "AT3G15500" then output gene2;
  if gene_id = "AT1G43160" then output gene3;
  if gene_id = "AT3G23240" then output gene4;
  if gene_id = "AT4G34410" then output gene5;
  drop gene_id;
run;


proc export data=gene1 outfile="!HOME/concannon/DTRA/at_rad_exp_DREB2C.csv" 
dbms=csv replace;
run;

proc export data=gene2 outfile="!HOME/concannon/DTRA/at_rad_exp_NAC055.csv" 
dbms=csv replace;
run;

proc export data=gene3 outfile="!HOME/concannon/DTRA/at_rad_exp_RAP2_6.csv" 
dbms=csv replace;
run;

proc export data=gene4 outfile="!HOME/concannon/DTRA/at_rad_exp_ERF1B.csv" 
dbms=csv replace;
run;

proc export data=gene5 outfile="!HOME/concannon/DTRA/at_rad_exp_ERF109.csv" 
dbms=csv replace;
run;


