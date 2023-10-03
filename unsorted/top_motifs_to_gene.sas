ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname tair10 

proc import datafile="!HOME/concannon/DTRA/top_motifs_genomewide.txt"
   out=motif_hits dbms=tab replace;
   getnames=no;
   guessingrows=1000000;
run;

proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/output/tair20_genes.bed" out=gene_coord dbms=tab replace;
 getnames=no;
guessingrows=max;
run;


data gene_1kb;
  set gene_coord;
  if var6 = "-" then do;
    promoter_start=VAR3+1000;
    promoter_stop=VAR3+1;
    downstream_start=VAR2-1;
    downstream_stop=VAR2-1000;
    end;
  else do;
    promoter_start=VAR2-1000;
    promoter_stop=VAR2-1;
    downstream_start=VAR3+1;
    downstream_stop=VAR3+1000;
    end;
  rename var1=chr var2=genebody_start var3=genebody_stop var6=strand var4=gene_id;
  drop VAR5;
run;

data gene_expand;
  set gene_1kb;
  if strand="-" then do;
      do pos = downstream_stop to promoter_start ;
      if pos >= downstream_stop and pos <= downstream_start then flag_downstream=1;
      else flag_downstream=0;
      if pos >= promoter_stop and pos <= promoter_start then flag_promoter=1;
      else flag_promoter=0;
      output;
    end;
    end;
    else do;
      do pos = promoter_start to downstream_stop;
      if pos <= downstream_stop and pos >= downstream_start then flag_downstream=1;
      else flag_downstream=0;
      if pos <= promoter_stop and pos >= promoter_start then flag_promoter=1;
      else flag_promoter=0;
      output;
      end;
  end;
run;

data motif_hits2;
  set motif_hits;
  length chr $2.;
  chr=compress(scan(VAR2, 1, " "));
  rename var1=motif_id var7=motif_sequence var3=pos var4=end_pos var6=score;
  drop var2;
run;

proc sort data=motif_hits2;
  by chr pos end_pos;
proc freq data=motif_hits2 noprint;
  tables chr*pos*end_pos / out=ctabs_cnt;
proc sort data=ctabs_cnt;
  by descending count;
run;


data _null_;
   set ctabs_cnt (obs=1 firstobs=1);
   call symputx("numMotif",count);
   stop;
run;


%put &numMotif.;


data motif_cat;
   array motif[&numMotif.] $30.;
   array scoreCat[&numMotif.] $30.;
   array seq[&numMotif.] $30.;
   retain motif1-motif&numMotif. ;
   retain scoreCat1-scoreCat&numMotif. ;
   retain seq1-seq&numMotif. ;
   set motif_hits2;
   by chr pos end_pos;
   if first.end_pos then do;
       call missing(of motif1-motif&numMotif. );
       records=0;
       end;
   records + 1;
   motif[records]=motif_id;
   scoreCat[records]=score;
   seq[records]=motif_sequence;
   if last.end_pos then output;
run;

data motif_cat2;
  set motif_cat;
  length motif_id_cat $100.;
  length motif_score_cat $100.;
  length motif_seq_cat $100.;
  motif_id_cat = catx("|", OF motif1-motif&numMotif. );
  motif_score_cat = catx("|", OF scoreCat1-scoreCat&numMotif. );
  motif_seq_cat = catx("|", OF seq1-seq&numMotif. );
  drop motif_id score motif_sequence records;
run;






proc sort data=motif_cat2;
  by chr pos;
proc sort data=gene_expand;
  by chr pos;
run;

data gene_expand_motif1 no_gene no_motif;
  merge gene_expand (in=in1) motif_cat2 (in=in2);
  by chr pos;
  if in1 and in2 then output gene_expand_motif1;
  else if in1 then output no_motif;
  else output no_gene;
run;


data no_gene2;
  set no_gene;
  keep chr pos end_pos motif_id_cat motif_score_cat motif_seq_cat;
  rename pos=start_pos end_pos=pos;
run;

proc sort data=no_gene2;
  by chr pos;
proc sort data=no_motif;
  by chr pos;
run;

data gene_expand_motif2;
  merge no_motif (in=in1) no_gene2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data gene2motif;
  set gene_expand_motif1 gene_expand_motif2;
  keep chr genebody_start genebody_stop gene_id strand promoter_start promoter_stop downstream_start downstream_stop
  flag_downstream flag_promoter motif_id_cat motif_score_cat motif_seq_cat;
run;

data gene2motif2;
  set gene2motif;
  length motif_id $32.;
  format score best12.;
  length motif_sequence $15.;
  do i = 1 by 1 while(scan(motif_id_cat,i,"|") ne "");
    motif_id=compress(scan(motif_id_cat,i,"|"));
    score=scan(motif_score_cat,i,"|") + 0;
    motif_sequence=compress(scan(motif_seq_cat,i,"|"));
    output;
  end;
  drop motif_id_cat motif_score_cat motif_seq_cat;
run;

data gene2motif3;
  set gene2motif2;
  length motif_family $32.;
  motif_family=compress(scan(motif_id,1,"-"));
run;

data all_motifs;
  set gene2motif3;
  keep gene_id motif_family;
run;

data promoter_motifs;
  set gene2motif3;
  where flag_promoter=1;
  keep gene_id motif_family;
run;

data downstream_motifs;
  set gene2motif3;
  where flag_downstream=1;
  keep gene_id motif_family;
run;


data genebody_motifs;
  set gene2motif3;
  where flag_downstream=0 and flag_promoter=0;
  keep gene_id motif_family;
run;


proc sort data=all_motifs nodup;
  by gene_id motif_family;
proc sort data=promoter_motifs nodup;
  by gene_id motif_family;
proc sort data=downstream_motifs nodup;
  by gene_id motif_family;
proc sort data=genebody_motifs nodup;
  by gene_id motif_family;
run;



%macro catMotif(location);

proc sort data=&location._motifs;
  by gene_id;
proc freq data=&location._motifs noprint;
  tables gene_id / out=ctabs_cnt;
proc sort data=ctabs_cnt;
  by descending count;
run;


data _null_;
   set ctabs_cnt (obs=1 firstobs=1);
   call symputx("numMotif",count);
   stop;
run;


%put &numMotif.;


data &location._motif_cat;
   array motif[&numMotif.] $30.;
   retain motif1-motif&numMotif. ;
   set &location._motifs;
   by gene_id;
   if first.gene_id then do;
       call missing(of motif1-motif&numMotif. );
       records=0;
       end;
   records + 1;
   motif[records]=motif_family;
   if last.gene_id then output;
run;

data &location._motif_cat2;
  set &location._motif_cat;
  length motif_id_cat_&location. $100.;
  motif_id_cat_&location. = catx("|", OF motif1-motif&numMotif. );
  drop motif_family motif1-motif&numMotif. records;
run;

%mend;



%catMotif(all);
%catMotif(promoter);
%catMotif(downstream);
%catMotif(genebody);


proc sort data=all_motif_cat2;
  by gene_id;
proc sort data=promoter_motif_cat2;
  by gene_id;
proc sort data=downstream_motif_cat2;
  by gene_id;
proc sort data=genebody_motif_cat2;
  by gene_id;
run;

data top_motif_flags_by_gene;
   merge all_motif_cat2 promoter_motif_cat2 downstream_motif_cat2 genebody_motif_cat2;
  by gene_id;
run;




data arabRNA.top_motif_flags_by_gene;
  set top_motif_flags_by_gene;
run;

  





  
