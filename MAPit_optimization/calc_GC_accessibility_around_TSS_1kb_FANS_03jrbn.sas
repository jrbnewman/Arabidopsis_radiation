/* For 22C and 4C, show and plot the average level of GC accessibility
   (calculated as 100U - 0U ) +/- 1kb around TSS for all genes 

doing this now for transcripts. I think what I want to do is pull in the transcript TSS sites
and then compare them against the HOMER annotations to see which looks "better"
let Mingqi decide which he likes better
idc
*/

libname cold '!PATCON/arabidopsis_wgbs_cold/sas_data';
libname coldloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

ods listing;
ods html close;

proc datasets lib=work kill noprint;
run;
quit;


* GC methylation data;
data gc_data;
  set coldloc.methylation_data;
  where site_type="GC" ;
run;

proc sort data=gc_data;
  by chr pos_end condition;
proc transpose data=gc_data out=gc_sbys;
  by chr pos_end;
  id condition;
  var perc_methyl;
run;


proc transpose data=gc_data out=gc_sbys10;
  where total_C >= 10;
  by chr pos_end;
  id condition;
  var perc_methyl;
run;

data gc_sbys_2;
  set gc_sbys;
  if _22C_0U ne . and _22C_100U ne . then _22C_100U_0U=_22C_100U-_22C_0U;
  else _22C_100U_0U=.;
  if _4C_0U ne . and _4C_100U ne . then _4C_100U_0U=_4C_100U-_4C_0U;
  else _4C_100U_0U=.;
  if _22C_0U ne . and _22C_100U ne . and _4C_0U ne . and _4C_100U ne . then do;
        _22C_100U_0U_common=_22C_100U - _22C_0U;
        _4C_100U_0U_common=_4C_100U - _4C_0U;
  end;
  rename pos_end=pos;
run;


data gc_sbys10_2;
  set gc_sbys10;
  if _22C_0U ne . and _22C_100U ne . then _22C_100U_0U=_22C_100U-_22C_0U;
  else _22C_100U_0U=.;
  if _4C_0U ne . and _4C_100U ne . then _4C_100U_0U=_4C_100U-_4C_0U;
  else _4C_100U_0U=.;
  if _22C_0U ne . and _22C_100U ne . and _4C_0U ne . and _4C_100U ne . then do;
        _22C_100U_0U_common=_22C_100U - _22C_0U;
        _4C_100U_0U_common=_4C_100U - _4C_0U;
  end;
  rename pos_end=pos;
run;

/* Get TSSs and regions from Ensembl TAIR10/Araport11 GTF annotation */

data exon2xs;
  set tair.tair20_exons_w_info;
  length transcript_id2 $50.;
  do i = 1 by 1 while(scan(transcript_id,i,"|") ^= "");
    transcript_id2 = compress(scan(transcript_id,i,"|"));
    output;
    end;
  keep chrom start stop strand exon_id transcript_id2;
  rename transcript_id2=transcript_id;
run;

*ID TSS;

proc sort data=exon2xs nodup;
  by transcript_id chrom start stop ;
run;

data tss_gene_plus;
  set exon2xs;
  by transcript_id;
  if first.transcript_id and strand="+" then output;
  keep transcript_id chrom start strand;
  rename chrom=chr start=tss;
run;

proc sort data=exon2xs;
  by transcript_id chrom descending stop descending start;
run;

data tss_gene_minus;
  set exon2xs;
  by transcript_id;
  if first.transcript_id and strand="-" then output;
  keep transcript_id chrom stop strand;
  rename chrom=chr stop=tss;
run;

data tss_gene;
  set tss_gene_plus tss_gene_minus;
    tss_1kb_start=tss-1000;
    tss_1kb_stop=tss+1000;
run;


proc sort data=tss_gene;
  by transcript_id;
run;

data tss_pos_index;
  set tss_gene;
  by transcript_id ;
  do pos = tss_1kb_start to tss_1kb_stop ;
  pos_abs=pos-tss;
  output;
  end;
 run;

data tss_pos_index2;
  set tss_pos_index;
  if strand="-" then do;
       pos_abs2=pos_abs * -1;
    output;
    end;
  else do;
    pos_abs2=pos_abs;
    output;
    end;
  drop pos_abs;
  rename pos_abs2=pos_abs;
run;


proc sort data=tss_pos_index2;
  by chr pos;
proc sort data=gc_sbys_2;
  by chr pos;
proc sort data=gc_sbys10_2;
  by chr pos;
run;

data tss_pos_gc_all;
  merge gc_sbys_2 (in=in1) tss_pos_index2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

data tss_pos_gc_10;
  merge gc_sbys10_2 (in=in1) tss_pos_index2 (in=in2);
  by chr pos;
  if in1 and in2;
run;

/* Also going to group every 10 positions */
data tss_pos_gc_all_2;
  set tss_pos_gc_all;
  grouped_pos=int(pos_abs/10) * 10;
run;

data tss_pos_gc_10_2;
  set tss_pos_gc_10;
  grouped_pos=int(pos_abs/10) * 10;
run;

/* For each, calculate the mean and SD */

proc sort data=tss_pos_gc_all_2;
  by pos_abs grouped_pos;
proc sort data=tss_pos_gc_10_2;
  by pos_abs grouped_pos;
run;

proc means data=tss_pos_gc_all_2 noprint;
  by pos_abs grouped_pos;
  var  _22C_100U _22C_100U_0U _4C_100U_0U _22C_100U_0U_common _4C_100U_0U_common FANS_0p5U FANS_1p5U FANS_25U FANS_5U;
  output out=mean_diff_tss_all_1000
  mean(_22C_100U)=mean_22C_100U  stddev(_22C_100U)=sd_22C_100U
  mean(_22C_100U_0U)=mean_22C_100U_0U  stddev(_22C_100U_0U)=sd_22C_100U_0U
  mean(_4C_100U_0U)=mean_4C_100U_0U  stddev(_4C_100U_0U)=sd_4C_100U_0U
  mean(_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(FANS_0p5U)=mean_FANS_0p5U  stddev(FANS_0p5U)=sd_FANS_0p5U
  mean(FANS_1p5U)=mean_FANS_1p5U   stddev(FANS_1p5U)=sd_FANS_1p5U
  mean(FANS_25U)=mean_FANS_25U   stddev(FANS_25U)=sd_FANS_25U
  mean(FANS_5U)=mean_FANS_5U  stddev(FANS_5U)=sd_FANS_5U;
run;

proc means data=tss_pos_gc_10_2 noprint;
  by pos_abs grouped_pos;
  var  _22C_100U _22C_100U_0U _4C_100U_0U _22C_100U_0U_common _4C_100U_0U_common FANS_0p5U FANS_1p5U FANS_25U FANS_5U;
  output out=mean_diff_tss_10X_1000
  mean(_22C_100U)=mean_22C_100U  stddev(_22C_100U)=sd_22C_100U
  mean(_22C_100U_0U)=mean_22C_100U_0U  stddev(_22C_100U_0U)=sd_22C_100U_0U
  mean(_4C_100U_0U)=mean_4C_100U_0U  stddev(_4C_100U_0U)=sd_4C_100U_0U
  mean(_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(FANS_0p5U)=mean_FANS_0p5U  stddev(FANS_0p5U)=sd_FANS_0p5U
  mean(FANS_1p5U)=mean_FANS_1p5U   stddev(FANS_1p5U)=sd_FANS_1p5U
  mean(FANS_25U)=mean_FANS_25U   stddev(FANS_25U)=sd_FANS_25U
  mean(FANS_5U)=mean_FANS_5U  stddev(FANS_5U)=sd_FANS_5U;
run;



proc sort data=mean_diff_tss_all_1000;
  by grouped_pos;
proc sort data=mean_diff_tss_10X_1000;
  by grouped_pos;
run;

proc means data=mean_diff_tss_all_1000 noprint;
  by grouped_pos;
  var mean_22C_100U mean_22C_100U_0U mean_4C_100U_0U mean_22C_100U_0U_common mean_4C_100U_0U_common mean_FANS_0p5U mean_FANS_1p5U mean_FANS_25U mean_FANS_5U;
  output out=mean_diff_tss_all_100
  mean(mean_22C_100U)=mean_22C_100U stddev(mean_22C_100U)=sd_22C_100U
  mean(mean_22C_100U_0U)=mean_22C_100U_0U  stddev(mean_22C_100U_0U)=sd_22C_100U_0U
  mean(mean_4C_100U_0U)=mean_4C_100U_0U  stddev(mean_4C_100U_0U)=sd_4C_100U_0U
  mean(mean_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(mean_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(mean_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(mean_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(mean_FANS_0p5U)=mean_FANS_0p5U   stddev(mean_FANS_0p5U)=sd_FANS_0p5U
  mean(mean_FANS_1p5U)=mean_FANS_1p5U   stddev(mean_FANS_1p5U)=sd_FANS_1p5U
  mean(mean_FANS_25U)=mean_FANS_25U   stddev(mean_FANS_25U)=sd_FANS_25U
  mean(mean_FANS_5U)=mean_FANS_5U   stddev(mean_FANS_5U)=sd_FANS_5U;
run;

proc means data=mean_diff_tss_10X_1000 noprint;
  by grouped_pos;
  var mean_22C_100U mean_22C_100U_0U mean_4C_100U_0U mean_22C_100U_0U_common mean_4C_100U_0U_common mean_FANS_0p5U mean_FANS_1p5U mean_FANS_25U mean_FANS_5U;;
  output out=mean_diff_tss_10X_100
  mean(mean_22C_100U)=mean_22C_100U stddev(mean_22C_100U)=sd_22C_100U
  mean(mean_22C_100U_0U)=mean_22C_100U_0U  stddev(mean_22C_100U_0U)=sd_22C_100U_0U
  mean(mean_4C_100U_0U)=mean_4C_100U_0U  stddev(mean_4C_100U_0U)=sd_4C_100U_0U
  mean(mean_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(mean_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(mean_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(mean_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(mean_FANS_0p5U)=mean_FANS_0p5U   stddev(mean_FANS_0p5U)=sd_FANS_0p5U
  mean(mean_FANS_1p5U)=mean_FANS_1p5U   stddev(mean_FANS_1p5U)=sd_FANS_1p5U
  mean(mean_FANS_25U)=mean_FANS_25U   stddev(mean_FANS_25U)=sd_FANS_25U
  mean(mean_FANS_5U)=mean_FANS_5U   stddev(mean_FANS_5U)=sd_FANS_5U;
run;

data mean_diff_tss_all_1000_2;
  set mean_diff_tss_all_1000;
  drop grouped_pos;
  rename pos_abs=pos;
run;

data mean_diff_tss_10X_1000_2;
  set mean_diff_tss_10X_1000;
  drop grouped_pos;
  rename pos_abs=pos;
run;



data mean_diff_tss_all_100_2;
  set mean_diff_tss_all_100;
  rename grouped_pos=pos;
run;

data mean_diff_tss_10X_100_2;
  set mean_diff_tss_10X_100;
  rename grouped_pos=pos;
run;


/* Export and plot in python */

proc export data=mean_diff_tss_all_1000_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_nobin_transcript_TSS.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_10X_1000_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_nobin_transcript_TSS.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_all_100_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_binned_transcript_TSS.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_10X_100_2
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_binned_transcript_TSS.csv"
   dbms=csv
   replace;
run;

/* Do the same as above, using HOMER annotations */

/* Get TSSs and regions from Ensembl TAIR10/Araport11 GTF annotation */



    data GC_annot   ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!PATCON/arabidopsis_wgbs_cold/analysis_output/HOMER_annotation/GC_sites_for_HOMER_annotated2.txt" delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat regionID $50.;
       informat Chr $3. ;
       informat Start best32. ;
       informat End best32. ;
       informat Strand $1. ;
       informat Peak_Score best32. ;
       informat Focus_Ratio_Region_Size $23. ;
       informat Annotation $100. ;
       informat Detailed_Annotation $100. ;
       informat Distance_to_TSS best32. ;
       informat Nearest_PromoterID $18. ;
       informat Entrez_ID best32. ;
       informat Nearest_Unigene $15. ;
       informat Nearest_Refseq $14. ;
       informat Nearest_Ensembl $15. ;
       informat Gene_Name $15. ;
       informat Gene_Alias $500. ;
       informat Gene_Description $200. ;
       informat Gene_Type $14. ;
        format regionID $50. ;
        format Chr $3. ;
        format Start best12. ;
        format End best12. ;
        format Strand $1. ;
        format Peak_Score best12. ;
        format Focus_Ratio_Region_Size $23. ;
        format Annotation $100. ;
        format Detailed_Annotation $100. ;
        format Distance_to_TSS best12. ;
        format Nearest_PromoterID $18. ;
        format Entrez_ID best32. ;
        format Nearest_Unigene $15. ;
        format Nearest_Refseq $14. ;
        format Nearest_Ensembl $15. ;
        format Gene_Name $15. ;
        format Gene_Alias $500. ;
        format Gene_Description $200. ;
        format Gene_Type $14. ;
     input
                 regionID $
                 Chr $
                 Start
                 End
                 Strand $
                 Peak_Score
                 Focus_Ratio_Region_Size $
                 Annotation $
                 Detailed_Annotation $
                 Distance_to_TSS
                 Nearest_PromoterID $
                 Entrez_ID
                 Nearest_Unigene $
                 Nearest_Refseq $
                 Nearest_Ensembl $
                 Gene_Name $
                 Gene_Alias $
                 Gene_Description $
                 Gene_Type $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;



data site2promoter;
  set GC_annot;
  length gene_id $20.;
  if abs(distance_to_tss) > 1000 then delete;
  if count(nearest_promoterID, "-T1") > 0 then gene_ID=compress(upcase(tranwrd(Nearest_PromoterID,"-T1","")));
  else gene_ID=compress(upcase(scan(Nearest_PromoterID,1,".")));
  keep regionID gene_ID chr start strand distance_to_tss nearest_promoterID;
  rename start=pos nearest_promoterID=transcript_id;
run;

/* Get strand for each transcript so I know what sites should be up/downstream of TSS */

data xs_strand;
  set tair.tair20_exons_w_info;
  keep strand gene_id;
  rename strand=xs_strand;
run;

proc sort data=xs_strand nodup;
  by gene_id;
proc sort data=site2promoter;
    by gene_id;
run;

data site2promoter2 no_strand no_site;
  merge site2promoter (in=in1) xs_strand (in=in2);
  by gene_id;
  if in1 and in2 then output site2promoter2;
  else if in1 then output no_strand;
  else output no_site;
run;

/* distance to TSS is relative to strand of transcript, so I don't need to flip anything!!! */

proc sort data=site2promoter nodup;
  by chr pos;
proc sort data=gc_sbys_2;
  by chr pos;
proc sort data=gc_sbys10_2;
  by chr pos;
run;

data tss_pos_gc_allH;
  merge gc_sbys_2 (in=in1) site2promoter (in=in2);
  by chr pos;
  if in1 and in2;
run;

data tss_pos_gc_10H;
  merge gc_sbys10_2 (in=in1) site2promoter (in=in2);
  by chr pos;
  if in1 and in2;
run;

/* Also going to group every 10 positions */
data tss_pos_gc_all_2H;
  set tss_pos_gc_allH;
  grouped_pos=int(distance_to_tss/10) * 10;
run;

data tss_pos_gc_10_2H;
  set tss_pos_gc_10H;
  grouped_pos=int(distance_to_tss/10) * 10;
run;

/* For each, calculate the mean and SD */

proc sort data=tss_pos_gc_all_2H;
  by distance_to_tss grouped_pos;
proc sort data=tss_pos_gc_10_2H;
  by distance_to_tss grouped_pos;
run;


proc means data=tss_pos_gc_all_2H noprint;
  by distance_to_tss grouped_pos;
  var  _22C_100U _22C_100U_0U _4C_100U_0U _22C_100U_0U_common _4C_100U_0U_common FANS_0p5U FANS_1p5U FANS_25U FANS_5U;
  output out=mean_diff_tss_all_1000H
  mean(_22C_100U)=mean_22C_100U  stddev(_22C_100U)=sd_22C_100U
  mean(_22C_100U_0U)=mean_22C_100U_0U  stddev(_22C_100U_0U)=sd_22C_100U_0U
  mean(_4C_100U_0U)=mean_4C_100U_0U  stddev(_4C_100U_0U)=sd_4C_100U_0U
  mean(_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(FANS_0p5U)=mean_FANS_0p5U  stddev(FANS_0p5U)=sd_FANS_0p5U
  mean(FANS_1p5U)=mean_FANS_1p5U   stddev(FANS_1p5U)=sd_FANS_1p5U
  mean(FANS_25U)=mean_FANS_25U   stddev(FANS_25U)=sd_FANS_25U
  mean(FANS_5U)=mean_FANS_5U  stddev(FANS_5U)=sd_FANS_5U;
run;

proc means data=tss_pos_gc_10_2H noprint;
  by distance_to_tss grouped_pos;
  var  _22C_100U _22C_100U_0U _4C_100U_0U _22C_100U_0U_common _4C_100U_0U_common FANS_0p5U FANS_1p5U FANS_25U FANS_5U;
  output out=mean_diff_tss_10X_1000H
  mean(_22C_100U)=mean_22C_100U  stddev(_22C_100U)=sd_22C_100U
  mean(_22C_100U_0U)=mean_22C_100U_0U  stddev(_22C_100U_0U)=sd_22C_100U_0U
  mean(_4C_100U_0U)=mean_4C_100U_0U  stddev(_4C_100U_0U)=sd_4C_100U_0U
  mean(_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(FANS_0p5U)=mean_FANS_0p5U  stddev(FANS_0p5U)=sd_FANS_0p5U
  mean(FANS_1p5U)=mean_FANS_1p5U   stddev(FANS_1p5U)=sd_FANS_1p5U
  mean(FANS_25U)=mean_FANS_25U   stddev(FANS_25U)=sd_FANS_25U
  mean(FANS_5U)=mean_FANS_5U  stddev(FANS_5U)=sd_FANS_5U;
run;



proc sort data=mean_diff_tss_all_1000H;
  by grouped_pos;
proc sort data=mean_diff_tss_10X_1000H;
  by grouped_pos;
run;

proc means data=mean_diff_tss_all_1000H noprint;
  by grouped_pos;
  var mean_22C_100U mean_22C_100U_0U mean_4C_100U_0U mean_22C_100U_0U_common mean_4C_100U_0U_common mean_FANS_0p5U mean_FANS_1p5U mean_FANS_25U mean_FANS_5U;
  output out=mean_diff_tss_all_100H
  mean(mean_22C_100U)=mean_22C_100U stddev(mean_22C_100U)=sd_22C_100U
  mean(mean_22C_100U_0U)=mean_22C_100U_0U  stddev(mean_22C_100U_0U)=sd_22C_100U_0U
  mean(mean_4C_100U_0U)=mean_4C_100U_0U  stddev(mean_4C_100U_0U)=sd_4C_100U_0U
  mean(mean_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(mean_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(mean_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(mean_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(mean_FANS_0p5U)=mean_FANS_0p5U   stddev(mean_FANS_0p5U)=sd_FANS_0p5U
  mean(mean_FANS_1p5U)=mean_FANS_1p5U   stddev(mean_FANS_1p5U)=sd_FANS_1p5U
  mean(mean_FANS_25U)=mean_FANS_25U   stddev(mean_FANS_25U)=sd_FANS_25U
  mean(mean_FANS_5U)=mean_FANS_5U   stddev(mean_FANS_5U)=sd_FANS_5U;
run;

proc means data=mean_diff_tss_10X_1000H noprint;
  by grouped_pos;
  var mean_22C_100U mean_22C_100U_0U mean_4C_100U_0U mean_22C_100U_0U_common mean_4C_100U_0U_common mean_FANS_0p5U mean_FANS_1p5U mean_FANS_25U mean_FANS_5U;;
  output out=mean_diff_tss_10X_100H
  mean(mean_22C_100U)=mean_22C_100U stddev(mean_22C_100U)=sd_22C_100U
  mean(mean_22C_100U_0U)=mean_22C_100U_0U  stddev(mean_22C_100U_0U)=sd_22C_100U_0U
  mean(mean_4C_100U_0U)=mean_4C_100U_0U  stddev(mean_4C_100U_0U)=sd_4C_100U_0U
  mean(mean_22C_100U_0U_common)=mean_22C_100U_0U_common  stddev(mean_22C_100U_0U_common)=sd_22C_100U_0U_common
  mean(mean_4C_100U_0U_common)=mean_4C_100U_0U_common  stddev(mean_4C_100U_0U_common)=sd_4C_100U_0U_common
  mean(mean_FANS_0p5U)=mean_FANS_0p5U   stddev(mean_FANS_0p5U)=sd_FANS_0p5U
  mean(mean_FANS_1p5U)=mean_FANS_1p5U   stddev(mean_FANS_1p5U)=sd_FANS_1p5U
  mean(mean_FANS_25U)=mean_FANS_25U   stddev(mean_FANS_25U)=sd_FANS_25U
  mean(mean_FANS_5U)=mean_FANS_5U   stddev(mean_FANS_5U)=sd_FANS_5U;
run;

data mean_diff_tss_all_1000_2H;
  set mean_diff_tss_all_1000H;
  drop grouped_pos;
  rename distance_to_tss=pos;
run;

data mean_diff_tss_10X_1000_2H;
  set mean_diff_tss_10X_1000H;
  drop grouped_pos;
  rename distance_to_tss=pos;
run;



data mean_diff_tss_all_100_2H;
  set mean_diff_tss_all_100H;
  rename grouped_pos=pos;
run;

data mean_diff_tss_10X_100_2H;
  set mean_diff_tss_10X_100H;
  rename grouped_pos=pos;
run;


/* Export and plot in python */

proc export data=mean_diff_tss_all_1000_2H
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_nobin_HOMER.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_10X_1000_2H
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_nobin_HOMER.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_all_100_2H
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_all_sites_binned_HOMER.csv"
   dbms=csv
   replace;
run;

proc export data=mean_diff_tss_10X_100_2H
   outfile="!PATCON/arabidopsis_wgbs_cold/analysis_output/FANS_GC_accessibility_1kb_TSS_10X_sites_binned_HOMER.csv"
   dbms=csv
   replace;
run;


