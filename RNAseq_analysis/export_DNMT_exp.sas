/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


%macro exportCounts(geneId,geneName);

data counts;
  set arabRNA.cpm_norm_counts_by_gene;
  length group $12.;
  length timepoint $4.;
  where gene_id = "&geneID.";
  log_cpm=log(cpm+1);
  if treatment="Mock" then group="Mock";
  if treatment="0.1gy" then group="10cGy";
  if treatment="1gy" then group="100cGy";
  if time = 1 then timepoint = "01h";
  else if time = 3 then timepoint = "03h";
  else if time = 24 then timepoint = "24h";
  else if time = 72 then timepoint = "72h";
  keep replicate timepoint group  cpm log_cpm;
run;	


proc sort data=counts;
  by  timepoint group replicate;
run;




proc export data=counts outfile="!HOME/concannon/DTRA/at_rad_exp_&geneName..csv" 
dbms=csv replace;
run;


%mend;


%exportCounts(AT1G69770,CMT3);
%exportCounts(AT1G80740,CMT1);
%exportCounts(AT3G17310,DRM3);
%exportCounts(AT4G14140,MET2);
%exportCounts(AT4G19020,CMT2);
%exportCounts(AT5G14620,DRM2);
%exportCounts(AT5G49160,MET1);

	
	
	
	
	
	
	
	
	
	
	
	
	
	

