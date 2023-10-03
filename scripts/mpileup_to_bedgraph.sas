ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


%macro importMpile(group,sample1,sample2,sample3);


   data WORK.REP1    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile "/blue/concannon/share/jnewman/mingqi_arab/bwa_mem_aln_pe/mpileups_genome_tair10/&sample1.-&group.-1.tsv" delimiter='09'x MISSOVER DSD lrecl=32767 ;
      informat VAR1 $200. ;
      informat VAR2 best32. ;
      informat VAR3 best32. ;
      format VAR1 $200. ;
      format VAR2 best12. ;
      format VAR3 best12. ;
   input
               VAR1 $
               VAR2
               VAR3
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

   data WORK.REP2    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile "/blue/concannon/share/jnewman/mingqi_arab/bwa_mem_aln_pe/mpileups_genome_tair10/&sample2.-&group.-2.tsv" delimiter='09'x MISSOVER DSD lrecl=32767 ;
      informat VAR1 $200. ;
      informat VAR2 best32. ;
      informat VAR3 best32. ;
      format VAR1 $200. ;
      format VAR2 best12. ;
      format VAR3 best12. ;
   input
               VAR1 $
               VAR2
               VAR3
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;


   data WORK.REP3    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile "/blue/concannon/share/jnewman/mingqi_arab/bwa_mem_aln_pe/mpileups_genome_tair10/&sample3.-&group.-3.tsv" delimiter='09'x MISSOVER DSD lrecl=32767 ;
      informat VAR1 $200. ;
      informat VAR2 best32. ;
      informat VAR3 best32. ;
      format VAR1 $200. ;
      format VAR2 best12. ;
      format VAR3 best12. ;
   input
               VAR1 $
               VAR2
               VAR3
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;



data stack;
  retain VAR1 start VAR2 VAR3;
  length VAR1 $200.;
  set rep1 rep2 rep3;
  start = VAR2 - 1;
  keep VAR1 start VAR2 VAR3;
run;

proc sort data=stack;
  by VAR1 start VAR2 ;
run;

proc means data=stack noprint;
  by var1 start var2;
  var var3;
  output out=mean_stack(drop=_TYPE_ _FREQ_) mean=;
run;

data header;
length VAR1 $200.;
VAR1="track name=&group. description=&group. name=bedGraph";
output;
VAR1="#chrom chromStart chromEnd value";
output;
run;


data header_stack;
  set header mean_stack;
run;

proc export data=header_stack
     outfile="/blue/concannon/share/jnewman/mingqi_arab/bwa_mem_aln_pe/mpileups_genome_tair10/&group..bedGraph"
     dbms=tab replace; putnames=no;
run;

%mend;

%importMpile(M-1,1,2,3);
%importMpile(0-1-1,4,5,6);
%importMpile(1-1,7,8,9);

%importMpile(M-3,10,11,12);
%importMpile(0-1-3,13,14,15);
%importMpile(1-3,16,17,18);

%importMpile(M-24,19,20,21);
%importMpile(0-1-24,22,23,24);
%importMpile(1-24,25,26,27);

%importMpile(M-72,28,29,30);
%importMpile(0-1-72,31,32,33);
%importMpile(1-72,34,35,36);


