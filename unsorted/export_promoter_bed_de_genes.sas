/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname tair '/nfshome/jrbnewman/concannon/useful_arabidopsis_data/TAIR10/sas_data';

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

data up_01_1 up_01_3 up_01_24 up_01_72
     dn_01_1 dn_01_3 dn_01_24 dn_01_72
     up_1_1  up_1_3  up_1_24  up_1_72
     dn_1_1  dn_1_3  dn_1_24  dn_1_72;
     set de_Results2;
     if (flag_01gy_v_mock_1h_fdr05=1 or (flag_01gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_01gy_1h > mean_cpm_Mock_1h) then output up_01_1;
     if (flag_01gy_v_mock_1h_fdr05=1 or (flag_01gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_01gy_1h < mean_cpm_Mock_1h) then output dn_01_1;
     if (flag_01gy_v_mock_3h_fdr05=1 or (flag_01gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_01gy_3h > mean_cpm_Mock_3h) then output up_01_3;
     if (flag_01gy_v_mock_3h_fdr05=1 or (flag_01gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_01gy_3h < mean_cpm_Mock_3h) then output dn_01_3;
     if (flag_01gy_v_mock_24h_fdr05=1 or (flag_01gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_01gy_24h > mean_cpm_Mock_24h) then output up_01_24;
     if (flag_01gy_v_mock_24h_fdr05=1 or (flag_01gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_01gy_24h < mean_cpm_Mock_24h) then output dn_01_24;
     if (flag_01gy_v_mock_72h_fdr05=1 or (flag_01gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_01gy_72h > mean_cpm_Mock_72h) then output up_01_72;
     if (flag_01gy_v_mock_72h_fdr05=1 or (flag_01gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_01gy_72h < mean_cpm_Mock_72h) then output dn_01_72;

     if (flag_1gy_v_mock_1h_fdr05=1 or (flag_1gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_1gy_1h > mean_cpm_Mock_1h) then output up_1_1;
     if (flag_1gy_v_mock_1h_fdr05=1 or (flag_1gy_1hr_on ^= flag_Mock_1hr_on )) and (mean_cpm_1gy_1h < mean_cpm_Mock_1h) then output dn_1_1;
     if (flag_1gy_v_mock_3h_fdr05=1 or (flag_1gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_1gy_3h > mean_cpm_Mock_3h) then output up_1_3;
     if (flag_1gy_v_mock_3h_fdr05=1 or (flag_1gy_3hr_on ^= flag_Mock_3hr_on )) and (mean_cpm_1gy_3h < mean_cpm_Mock_3h) then output dn_1_3;
     if (flag_1gy_v_mock_24h_fdr05=1 or (flag_1gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_1gy_24h > mean_cpm_Mock_24h) then output up_1_24;
     if (flag_1gy_v_mock_24h_fdr05=1 or (flag_1gy_24hr_on ^= flag_Mock_24hr_on )) and (mean_cpm_1gy_24h < mean_cpm_Mock_24h) then output dn_1_24;
     if (flag_1gy_v_mock_72h_fdr05=1 or (flag_1gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_1gy_72h > mean_cpm_Mock_72h) then output up_1_72;
     if (flag_1gy_v_mock_72h_fdr05=1 or (flag_1gy_72hr_on ^= flag_Mock_72hr_on )) and (mean_cpm_1gy_72h < mean_cpm_Mock_72h) then output dn_1_72;
keep gene_id;
run;


data exons;
  set tair.tair20_exons;
  keep chrom start stop strand gene_id;
run;

proc sort data=exons;
  by chrom strand gene_id start stop;
proc means data=exons noprint;
  by chrom strand gene_id;
  var start stop;
  output out=gene_coord min(start)=start max(stop)=stop;
run;


data promoter;
  retain chrom promoter_start promoter_stop gene_id score strand;
  set gene_coord;
  if strand = "+" then do;
     promoter_start = start - 1000;
     promoter_stop = start + 100;
  end;
  else if strand = "-" then do;
     promoter_stop = stop + 1000;
     promoter_start = stop - 100;
  end;
  score = ".";
  keep gene_id strand promoter_start promoter_stop score chrom;
run;

proc sort data=promoter;
  by gene_id;
proc sort data=up_01_1;
  by gene_id;
proc sort data=up_01_3;
  by gene_id;
proc sort data=up_01_24;
  by gene_id;
proc sort data=up_01_72;
  by gene_id;
proc sort data=dn_01_1;
  by gene_id;
proc sort data=dn_01_3;
  by gene_id;
proc sort data=dn_01_24;
  by gene_id;
proc sort data=dn_01_72;
  by gene_id;
proc sort data=up_1_1;
  by gene_id;
proc sort data=up_1_3;
  by gene_id;
proc sort data=up_1_24;
  by gene_id;
proc sort data=up_1_72;
  by gene_id;
proc sort data=dn_1_1;
  by gene_id;
proc sort data=dn_1_3;
  by gene_id;
proc sort data=dn_1_24;
  by gene_id;
proc sort data=dn_1_72;
  by gene_id;
run;

data up_01_1_promoter; merge promoter (in=in1) up_01_1 (in=in2); by gene_id; if in1 and in2; run;
data up_01_3_promoter; merge promoter (in=in1) up_01_3 (in=in2); by gene_id; if in1 and in2; run;
data up_01_24_promoter; merge promoter (in=in1) up_01_24 (in=in2); by gene_id; if in1 and in2; run;
data up_01_72_promoter; merge promoter (in=in1) up_01_72 (in=in2); by gene_id; if in1 and in2; run;
data dn_01_1_promoter; merge promoter (in=in1) dn_01_1 (in=in2); by gene_id; if in1 and in2; run;
data dn_01_3_promoter; merge promoter (in=in1) dn_01_3 (in=in2); by gene_id; if in1 and in2; run;
data dn_01_24_promoter; merge promoter (in=in1) dn_01_24 (in=in2); by gene_id; if in1 and in2; run;
data dn_01_72_promoter; merge promoter (in=in1) dn_01_72 (in=in2); by gene_id; if in1 and in2; run;

data up_1_1_promoter; merge promoter (in=in1) up_1_1 (in=in2); by gene_id; if in1 and in2; run;
data up_1_3_promoter; merge promoter (in=in1) up_1_3 (in=in2); by gene_id; if in1 and in2; run;
data up_1_24_promoter; merge promoter (in=in1) up_1_24 (in=in2); by gene_id; if in1 and in2; run;
data up_1_72_promoter; merge promoter (in=in1) up_1_72 (in=in2); by gene_id; if in1 and in2; run;
data dn_1_1_promoter; merge promoter (in=in1) dn_1_1 (in=in2); by gene_id; if in1 and in2; run;
data dn_1_3_promoter; merge promoter (in=in1) dn_1_3 (in=in2); by gene_id; if in1 and in2; run;
data dn_1_24_promoter; merge promoter (in=in1) dn_1_24 (in=in2); by gene_id; if in1 and in2; run;
data dn_1_72_promoter; merge promoter (in=in1) dn_1_72 (in=in2); by gene_id; if in1 and in2; run;


data de_01_1_promoter; set up_01_1_promoter dn_01_1_promoter; run;
data de_01_3_promoter; set up_01_3_promoter dn_01_3_promoter; run;
data de_01_24_promoter; set up_01_24_promoter dn_01_24_promoter; run;
data de_01_72_promoter; set up_01_72_promoter dn_01_72_promoter; run;

data de_1_1_promoter; set up_1_1_promoter dn_1_1_promoter; run;
data de_1_3_promoter; set up_1_3_promoter dn_1_3_promoter; run;
data de_1_24_promoter; set up_1_24_promoter dn_1_24_promoter; run;
data de_1_72_promoter; set up_1_72_promoter dn_1_72_promoter; run;


data de_1_promoter; set de_01_1_promoter de_1_1_promoter; run;
data de_3_promoter; set de_01_3_promoter de_1_3_promoter; run;
data de_24_promoter; set de_01_24_promoter de_1_24_promoter; run;
data de_72_promoter; set de_01_72_promoter de_1_72_promoter; run;

data de_any_promoter;  set de_1_promoter de_3_promoter de_24_promoter de_72_promoter; run;

%macro sortExport(datain);

proc sort data=&datain. nodup;
by gene_id;
run;

proc sort data=&datain. ;
by chrom promoter_start promoter_stop;
run;

proc export data=&datain. outfile="/TB14/TB14/sandbox/dtra_sandbox/At_DEG_&datain..bed" 
dbms=tab replace;
putnames=no;
run;

%mend;

%sortExport(up_01_1_promoter); %sortExport(up_01_3_promoter); %sortExport(up_01_24_promoter); %sortExport(up_01_72_promoter);
%sortExport(dn_01_1_promoter); %sortExport(dn_01_3_promoter); %sortExport(dn_01_24_promoter); %sortExport(dn_01_72_promoter);
%sortExport(de_01_1_promoter); %sortExport(de_01_3_promoter); %sortExport(de_01_24_promoter); %sortExport(de_01_72_promoter);

%sortExport(up_1_1_promoter); %sortExport(up_1_3_promoter); %sortExport(up_1_24_promoter); %sortExport(up_1_72_promoter);
%sortExport(dn_1_1_promoter); %sortExport(dn_1_3_promoter); %sortExport(dn_1_24_promoter); %sortExport(dn_1_72_promoter);
%sortExport(de_1_1_promoter); %sortExport(de_1_3_promoter); %sortExport(de_1_24_promoter); %sortExport(de_1_72_promoter);

%sortExport(de_1_promoter); %sortExport(de_3_promoter); %sortExport(de_24_promoter); %sortExport(de_72_promoter);
%sortExport(de_any_promoter);







