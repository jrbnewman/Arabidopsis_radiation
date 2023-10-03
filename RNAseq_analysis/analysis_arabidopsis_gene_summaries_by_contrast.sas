/* Make gene-level summaries for each contrast
I want the following:
(1) GeneID/Symbol
(2) Number of fusions
(3) Number of tested fusions
(4) Concatenation of tested fusions IDs
(5) Concatenation of tested fusions U/D
(6) Number nominally significant total
(7) Number nominally significant and UP
(8) Number nominally significant and DOWN
(9) Average coverage over gene (tested fusions only)
(10) Average estimate over gene
(11) Average LSmeans over gene
(12) Overall gene up/down
(13) Flag if any one fusion in gene is significant
(14) Flag if any one fusion in gene is significant and UP
(15) Flag if any one fusion in gene is significant and DOWN

two sets: (1) no multigene (2) with multigene

*/

/* NON-MULTIGENE SET */

ods listing; ods html close;
libname rs "!PATCON/arabidopsis/sas_data";
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

/* Make gene-to-fusion index -- we can refer back to this later */

* ID genes only comprised of multigene fusions;

data mult nomult;
  set tair.tair20_fusion_si_info;
  if flag_multigene=1 then output mult;
  else output nomult;
  keep primary_Fbgn;
   rename primary_fbgn=gene_id;
run;

proc sort data=mult nodup;
   by gene_id;
proc sort data=nomult nodup;
   by gene_id;
run;

data mult_v_nomult;
  merge mult (in=in1) nomult (in=in2);
  by gene_id;
  if in1 then flag_has_mult=1; else flag_has_mult=0;
  if in2 then flag_no_mult=1; else flag_no_mult=0;
run;

proc freq data=mult_v_nomult;
   tables flag_has_mult*flag_no_mult;
run;

data gene2fus;
  set tair.tair20_fusion_si_info;
  if flag_multigene=0;
  keep fusion_id primary_Fbgn;
  rename primary_fbgn=gene_id;
run;

proc sort data=gene2fus nodup;
   by fusion_id gene_id;
run;

proc freq data=gene2fus noprint;
   tables gene_id / out=check;
proc sort data=check;
   by descending count;
run;

%macro geneSummary(contrast,conlabel);

/* Get number of tested fusions */

data tested_fus;
  set rs.arab_fus_cntrs_constr_non;
  where label=&conlabel.;
  keep fusion_id;
run;

proc sort data=tested_fus;
   by fusion_id;
proc sort data=gene2fus;
   by fusion_id;
run;

data tested_fus2gene;
  merge tested_fus (in=in1) gene2fus (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc freq data=tested_fus2gene noprint;
  table gene_id / out=tested_fus_per_gene;
run;

data tested_fus_per_gene2;
  set tested_fus_per_gene;
  keep gene_id count;
  rename count=num_tested_exons;
run;

/* Concatenation of tested fusions IDs -- I want this in genomic order */

data fus2coord;
  set tair.tair20_fusion_si_info_unique;
  keep chr fusion_start fusion_stop fusion_id;
run;

proc sort data=tested_fus2gene;
  by fusion_id;
proc sort data=fus2coord;
  by fusion_id;
run;

data fus2gene_w_coord;
  merge fus2coord (in=in1) tested_fus2gene (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus2gene_w_coord;
   by gene_id fusion_id;
proc freq data=fus2gene_w_coord noprint;
   tables gene_id / out=fus_count;
proc sort data=fus_count;
   by descending count;
run;

%local fusCount;

data _null_;
  set fus_count (firstobs=1 obs=1);
  call symputx('fusCount',count);
  stop;
run;

%put &fusCount.;


proc sort data=fus2gene_w_coord;
  by gene_id fusion_start fusion_stop ;
run;

data fus_id_cat;
   array fusions[&fusCount.] $10.;
   retain fusions1-fusions&fusCount. ; 
   set fus2gene_w_coord;
   by gene_id;
   if first.gene_id then do;
       call missing(of fusions1-fusions&fusCount.);
       records=0;
       end;
   records + 1;
   fusions[records]=fusion_id;
   if last.gene_id then output;
run;

data fus_id_cat2;
   set fus_id_cat;
   length exonic_region_id_cat $1000.;
   exonic_region_id_cat = catx("|", OF fusions1-fusions&fusCount.);
   keep gene_id exonic_region_id_cat records ; 
   rename records = num_exonic_regions_tested;
run;

/* Concatenation of tested fusions U/D */

data signs;
  set rs.arab_sign_by_contrast;
  keep fusion_id sign_&contrast.;
run;

proc sort data=signs;
  by fusion_id;
proc sort data=fus2gene_w_coord;
  by fusion_id;
run;

data fus2gene_w_sign;
  merge fus2gene_w_coord (in=in1) signs (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=fus2gene_w_sign;
  by gene_id fusion_start fusion_stop ;
run;

data fus_sign_cat;
   array fusions[&fusCount.] $1.;
   retain fusions1-fusions&fusCount. ; 
   set fus2gene_w_sign;
   by gene_id;
   if first.gene_id then do;
       call missing(of fusions1-fusions&fusCount.);
       records=0;
       end;
   records + 1;
   fusions[records]=sign_&contrast.;
   if last.gene_id then output;
run;

data fus_sign_cat2;
   set fus_sign_cat;
   length exon_sign_cat $100.;
   exon_sign_cat = catt(OF fusions1-fusions&fusCount.);
   keep gene_id exon_sign_cat ; 
run;

/* Number nominally significant total */

data sig_total;
   set rs.arab_fus_cntrs_constr_non;
   where label=&conlabel.;
   if ProbF ne . and  ProbF < 0.05;
   keep fusion_id;
run;

proc sort data=sig_total;
   by fusion_id;
proc sort data=tested_fus2gene;
   by fusion_id;
run;

data sig_fus2gene;
  merge tested_fus2gene (in=in1) sig_total (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc freq data=sig_fus2gene noprint;
  tables gene_id / out=sig_count;
run;

data sig_count2;
  set sig_count;
  keep gene_id count;
  rename count=num_exons_nom_p05;
run;

/* Number nominally significant and UP, DOWN */

data fus_up_down;
  set rs.arab_sign_by_contrast;
  keep fusion_id sign_&contrast.; 
run;

proc sort data=fus_up_down;
   by fusion_id;
proc sort data=tested_fus2gene;
   by fusion_id;
run;

data sign_fus2gene;
  merge tested_fus2gene (in=in1) fus_up_down (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc freq data=sign_fus2gene noprint;
  where sign_&contrast.="U";
  tables gene_id / out=up_count;
run;

proc freq data=sign_fus2gene noprint;
  where sign_&contrast.="D";
  tables gene_id / out=down_count;
run;

data up_count2;
  set up_count;
  keep gene_id count;
  rename count=num_exons_up_reg;
run;

data down_count2;
  set down_count;
  keep gene_id count;
  rename count=num_exons_down_reg;
run;

/* Average estimate over gene */

data estimates;
   set rs.arab_fus_cntrs_estim_non;
   where label=&conlabel.;
   keep fusion_id estimate;
run;

proc sort data=estimates;
   by fusion_id;
proc sort data=tested_fus2gene;
   by fusion_id;
run;

data est_fus2gene;
  merge estimates (in=in1) tested_fus2gene (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc sort data=est_fus2gene;
  by gene_id;
proc means data=est_fus2gene noprint;
  by gene_id;
  var estimate;
  output out=mean_est_by_gene mean=mean_estimate;
run;

/* Overall gene up/down */
data mean_est_by_gene2;
  set mean_est_by_gene;
  length overall_sign $1.;
  if mean_estimate < 0 then overall_sign="D";
  else if mean_estimate > 0 then overall_sign="U";
  else overall_sign="N";
  keep gene_id mean_estimate overall_sign;
run;

/* Merge all together */

proc sort data=tested_fus_per_gene2;
   by gene_id;
proc sort data=fus_id_cat2;
   by gene_id;
proc sort data=fus_sign_cat2;
   by gene_id;
proc sort data=sig_count2;
   by gene_id;
proc sort data=up_count2;
   by gene_id;
proc sort data=down_count2;
   by gene_id;
proc sort data=mean_est_by_gene2;
   by gene_id;
run;

data gene_summary;
   merge tested_fus_per_gene2 (in=in1) fus_id_cat2 (in=in2) fus_sign_cat2 (in=in3)
         sig_count2 (in=in4) up_count2 (in=in5) down_count2 (in=in6) mean_est_by_gene2 (in=in7);
   by gene_id;
   if not in4 then num_exons_nom_p05=0;
   if not in5 then num_exons_up_reg=0;
   if not in6 then num_exons_down_reg=0;
   if in1;
run;

data flag_genes;
  set gene_summary;
  if num_exons_nom_p05 > 0 then flag_gene_sig_exons=1; else flag_gene_sig_exons=0;
  if num_exons_up_reg > 0 then flag_gene_up_exons=1; else flag_gene_up_exons=0;
  if num_exons_down_reg > 0 then flag_gene_down_exons=1; else flag_gene_down_exons=0;
run;

data sym2gene;
  set tair.gene2symbol_ens;
run;

proc sort data=sym2gene nodup;
   by gene_id;
proc sort data=flag_genes;
   by gene_id;
run;

data gene_summary_w_sym;
  merge sym2gene (in=in1) flag_genes (in=in2);
  by gene_id;
  if in2;
run;

data rs.gene_sum_&contrast.;
  set gene_summary_w_sym;
run;

proc export data=gene_summary_w_sym
     outfile="!PATCON/arabidopsis/analysis_output/gene_summary_contrast_&contrast._nomulti.csv"
     dbms=csv replace;
run;

%mend;



%geneSummary(01gy_v_Mock_1h,"0.1gy-Mock: 1h");
%geneSummary(1gy_v_Mock_1h,"1gy-Mock: 1h");
%geneSummary(1gy_v_01gy_1h,"1gy-0.1gy: 1h");

%geneSummary(01gy_v_Mock_24h,"0.1gy-Mock: 24h");
%geneSummary(1gy_v_Mock_24h,"1gy-Mock: 24h");
%geneSummary(1gy_v_01gy_24h,"1gy-0.1gy: 24h");

%geneSummary(01gy_v_Mock_3h,"0.1gy-Mock: 3h");
%geneSummary(1gy_v_Mock_3h,"1gy-Mock: 3h");
%geneSummary(1gy_v_01gy_3h,"1gy-0.1gy: 3h");

%geneSummary(01gy_M_72h,"0.1gy-Mock: 72h");
%geneSummary(1gy_M_72h,"1gy-Mock: 72h");
%geneSummary(1gy_v_01gy_72h,"1gy-0.1gy: 72h");

%geneSummary(Mock_3h_v_1h,"Mock: 3h-1h");
%geneSummary(Mock_24h_v_1h,"Mock: 24h-1h");
%geneSummary(Mock_72h_v_1h,"Mock: 72h-1h");
%geneSummary(Mock_24h_v_3h,"Mock: 24h-3h");
%geneSummary(Mock_72h_v_3h,"Mock: 72h-3h");
%geneSummary(Mock_72h_v_24h,"Mock: 72h-24h");

%geneSummary(01gy_3h_v_1h,"0.1gy: 3h-1h");
%geneSummary(01gy_24h_v_1h,"0.1gy: 24h-1h");
%geneSummary(01gy_72h_v_1h,"0.1gy: 72h-1h");
%geneSummary(01gy_24h_v_3h,"0.1gy: 24h-3h");
%geneSummary(01gy_72h_v_3h,"0.1gy: 72h-3h");
%geneSummary(01gy_72h_v_24h,"0.1gy: 72h-24h");

%geneSummary(1gy_3h_v_1h,"1gy: 3h-1h");
%geneSummary(1gy_24h_v_1h,"1gy: 24h-1h");
%geneSummary(1gy_72h_v_1h,"1gy: 72h-1h");
%geneSummary(1gy_24h_v_3h,"1gy: 24h-3h");
%geneSummary(1gy_72h_v_3h,"1gy: 72h-3h");
%geneSummary(1gy_72h_v_24h,"1gy: 72h-24h");

%geneSummary(01gy_v_Mock_3h_sub_1h,"0.1gy 3h-1h = Mock 3h-1h");
%geneSummary(01gy_v_Mock_24h_sub_1h,"0.1gy 24h-1h = Mock 24h-1h");
%geneSummary(01gy_v_Mock_72h_sub_1h,"0.1gy 72h-1h = Mock 72h-1h");

%geneSummary(1gy_v_Mock_3h_sub_1h,"1gy 3h-1h = Mock 3h-1h");
%geneSummary(1gy_v_Mock_24h_sub_1h,"1gy 24h-1h = Mock 24h-1h");
%geneSummary(1gy_v_Mock_72h_sub_1h,"1gy 72h-1h = Mock 72h-1h");

12345678901324567890123456789012
gene_summary_01gy_v_Mock_72h_sub_1h
