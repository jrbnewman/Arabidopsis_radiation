proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/downloaded_files/TAIR10_TE.gtf"
out=gtf dbms=tab replace;
getnames=no;
guessingrows=max;
run;


data gtf2;
  set gtf;
  length gene $20.;
  length transcript $20.;
  gene=scan(scan(var9,1,";"),2," ");
  transcript=scan(scan(var9,2,";"),2," ");
run;

data gtf3;
  set gtf2;
  length gene2 $50.;
  gene2=compress(catt('"',tranwrd(gene,'"',''),"/",tranwrd(transcript,'"',''),'"'));
run;

data gtf4;
  set gtf3;
  length VAR9_2 $255.;
  VAR9_2=cat("gene_id ",compress(gene2),'; ', scan(var9,2,";"), '; ',scan(var9,3,";"), '; ',scan(var9,4,";"));
run;

data gtf_export;
  set gtf4;
  drop VAR9 gene transcript gene2;
run;

data bed_export;
  set gtf4;
  length gene3 $50.;
  gene3=compress(tranwrd(gene2,'"',''));
  keep VAR1 VAR4 VAR5 gene3;
run;

proc export data=gtf_export outfile="!HOME/concannon/useful_arabidopsis_data/TAIR10/output/TE_fragments.gtf"
    dbms=tab replace;
    putnames=no;
run;

proc export data=bed_export outfile="!HOME/concannon/useful_arabidopsis_data/TAIR10/output/TE_fragments.bed"
    dbms=tab replace;
    putnames=no;
run;



data _null_;
set gtf_export;
file "!HOME/concannon/useful_arabidopsis_data/TAIR10/output/TE_fragments.gtf" dlm='09'x ;
put VAR1 VAR2 VAR3 VAR4 VAR5 VAR6 VAR7 VAR8 VAR9_2  ;
run;

 

