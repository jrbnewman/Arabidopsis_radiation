libname brassRNA "!HOME/concannon/DTRA/brassica/sas_data";

/* Export counts for expressed DNMT genes 
   format should be in long format so I can group things properly 
*/

%macro exportData(geneID);

data counts;
  set brassRNA.cpm_norm_counts_by_gene;
  where gene_id = "&geneID.";
  if treatment="Mock" then order=1;
  if treatment="1.4cGy" then order=2;
  if treatment="0.1Gy" then order=3;
  if treatment="1Gy" then order=4;
  keep sample_id treatment time gene_id cpm order rep;
run;

proc sort data=counts;
  by order time rep;
run;

proc export data=counts outfile="!HOME/concannon/DTRA/Br_DNMT_cpm_counts_&geneID..csv"
dbms=csv replace;
run;



%mend;


%exportData(BraA02g06530R_gene);
%exportData(BraA09g52360R_gene);
%exportData(BraA10g23200R_gene);
%exportData(BraA10g23210R_gene);
%exportData(BraA05g28670R_gene);
%exportData(BraA01g03230R_gene);

