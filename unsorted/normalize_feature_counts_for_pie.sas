/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DMRs_min_5_sites_for_HOMER_annotation.txt"
   out=dmr_annot dbms=tab replace;
   guessingrows=max;
run;

proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER_annotation.txt"
   out=dar_annot dbms=tab replace;
   guessingrows=max;
run;


data dmr_annot2;
  set dmr_annot;
  length comparison $12.;
  length site_type $4.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  length feature $20.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'|'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'|'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'|'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'|')) + 0;
  feature = scan(annotation, 1, " ");
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;

data dar_annot2;
  set dar_annot;
  length comparison $12.;
  length site_type $4.;
  length chrom $3.;
  format region_num best12.;
  length geneID $100.;
  comparison=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,1,'|'));
  site_type=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,2,'|'));
  chrom=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,3,'|'));
  region_num=compress(scan(PeakID__cmd_annotatePeaks_pl__bl,4,'|')) + 0;
  feature = scan(annotation, 1, " ");
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"intron") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(tranwrd(upcase(scan(Nearest_PromoterID,1,".")),"-T1",""));  
   drop PeakID__cmd_annotatePeaks_pl__bl chr;
   rename chrom=chr;
run;


data results_by_dmr;
   set wgbsA.results_by_dmr_5sites;
run;

data results_by_dar;
   set wgbsA.results_by_dar_5sites;
run;

proc freq data=results_by_dar;
  tables comparison*flag_fdr05;
run;


proc sort data=dmr_annot2;
  by comparison site_type chr  region_num;
proc sort data=dar_annot2;
  by comparison site_type chr  region_num;
proc sort data=results_by_dmr;
  by comparison site_type chr  region_num;
proc sort data=results_by_dar;
  by comparison site_type chr  region_num;
run;

data dmr_w_annot;
  merge dmr_annot2 (in=in1) results_by_dmr (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;

data dar_w_annot;
  merge dar_annot2 (in=in1) results_by_dar (in=in2);
  by comparison site_type chr region_num;
  if in1 and in2;
run;


data dmr_w_annot2;
  set dmr_w_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;


data dar_w_annot2;
  set dar_w_annot;
  if mean_methyl_diff < 0 then flag_direction=-1;
  else if  mean_methyl_diff > 0 then flag_direction=1;
  else mean_methyl_diff = 0;
run;

/* Count DMRs */

data all_01_up all_01_dn all_1_up all_1_dn
     cg_01_up cg_01_dn  cg_1_up cg_1_dn
     chg_01_up chg_01_dn  chg_1_up chg_1_dn
     chh_01_up chh_01_dn  chh_1_up chh_1_dn;
     set dmr_w_annot2;
     /* all DMRs by site type and comparison */
     if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output all_01_up;
     if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output all_01_dn;
     if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output all_1_up;
     if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output all_1_dn;

     if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output all_01_up;
     if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output all_01_dn;
     if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output all_1_up;
     if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output all_1_dn;
     if site_type="CG" then do;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output cg_01_up;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output cg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output cg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output cg_1_dn;

         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output cg_01_up;
         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output cg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output cg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output cg_1_dn;
     end;
     if site_type="CHG" then do;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chg_01_up;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chg_1_dn;

         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output chg_01_up;
         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output chg_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output chg_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output chg_1_dn;
     end;
     if site_type="CHH" then do;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chh_01_up;
         if comparison = "0Gy_vs_01G" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chh_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=1 then output chh_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_DMC10_ge2 = 1 and flag_direction=-1 then output chh_1_dn;

         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=1 then output chh_01_up;
         if comparison = "0Gy_vs_01G" and flag_fdr05=1 and flag_direction=-1 then output chh_01_dn;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=1 then output chh_1_up;
         if comparison = "0Gy_vs_1Gy" and flag_fdr05=1 and flag_direction=-1 then output chh_1_dn;
     end;
run;

proc sort data=all_01_up nodup; by _all_;
proc sort data=all_01_dn nodup; by _all_;
proc sort data=all_1_up nodup; by _all_;
proc sort data=all_1_dn nodup; by _all_;
proc sort data=cg_01_up nodup; by _all_;
proc sort data=cg_01_dn nodup; by _all_;
proc sort data=cg_1_up nodup; by _all_;
proc sort data=cg_1_dn nodup; by _all_;
proc sort data=chg_01_up nodup; by _all_;
proc sort data=chg_01_dn nodup; by _all_;
proc sort data=chg_1_up nodup; by _all_;
proc sort data=chg_1_dn nodup; by _all_;
proc sort data=chh_01_up nodup; by _all_;
proc sort data=chh_01_dn nodup; by _all_;
proc sort data=chh_1_up nodup; by _all_;
proc sort data=chh_1_dn nodup; by _all_;
run;


/* Count DARs */

data gc_01_up gc_01_dn gc_1_up gc_1_dn;
     set dar_w_annot2;
     /* all DMRs by site type and comparison */
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output gc_01_up;
     if comparison = "01Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output gc_01_dn;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=1 then output gc_1_up;
     if comparison = "1Gy_0Gy" and flag_DAC10_ge2 = 1 and flag_direction=-1 then output gc_1_dn;

     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output gc_01_up;
     if comparison = "01Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output gc_01_dn;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=1 then output gc_1_up;
     if comparison = "1Gy_0Gy" and flag_fdr05=1 and flag_direction=-1 then output gc_1_dn;
run;


proc sort data=all_01_up nodup; by feature;
proc sort data=all_01_dn nodup; by feature;
proc sort data=all_1_up nodup; by feature;
proc sort data=all_1_dn nodup; by feature;
proc sort data=cg_01_up nodup; by feature;
proc sort data=cg_01_dn nodup; by feature;
proc sort data=cg_1_up nodup; by feature;
proc sort data=cg_1_dn nodup; by feature;
proc sort data=chg_01_up nodup; by feature;
proc sort data=chg_01_dn nodup; by feature;
proc sort data=chg_1_up nodup; by feature;
proc sort data=chg_1_dn nodup; by feature;
proc sort data=chh_01_up nodup; by feature;
proc sort data=chh_01_dn nodup; by feature;
proc sort data=chh_1_up nodup; by feature;
proc sort data=chh_1_dn nodup; by feature;
proc sort data=gc_01_up nodup; by feature;
proc sort data=gc_01_dn nodup; by feature;
proc sort data=gc_1_up nodup; by feature;
proc sort data=gc_1_dn nodup; by feature;
run;

%macro normMotif(input);

proc sort data=&input. nodup; by feature;
run;


proc means data=&input. noprint;
   by feature;
   var _FREQ_;
   output out=&input._cnt sum=num_motifs;
run;

data &input._norm;
  set &input._cnt;
  norm_size = _FREQ_ / num_motifs;
run;
 
proc print data=&input._norm;
run;


%mend;
%normMotif(gc_01_up); %normMotif(gc_01_dn); %normMotif(gc_1_up); %normMotif(gc_1_dn);

%normMotif(cg_01_up); %normMotif(cg_01_dn); %normMotif(cg_1_up); %normMotif(cg_1_dn);
%normMotif(chg_01_up); %normMotif(chg_01_dn); %normMotif(chg_1_up); %normMotif(chg_1_dn);

%normMotif(chh_01_up); %normMotif(chh_01_dn); %normMotif(chh_1_up); %normMotif(chh_1_dn);


/*

%normMotif(gc_01_up);


                                     num_      norm_
feature         _TYPE_    _FREQ_    motifs      size

Intergenic         0       2382      14694    0.16211
TTS                0        251       1471    0.17063
exon               0        525       3056    0.17179
intron             0        161        933    0.17256
promoter-TSS       0        252       1458    0.17284



%normMotif(gc_01_dn);

                                     num_      norm_
feature         _TYPE_    _FREQ_    motifs      size

Intergenic         0         8        41      0.19512
TTS                0         3        45      0.06667
exon               0         6        50      0.12000
promoter-TSS       0         1         5      0.20000


%normMotif(gc_1_up);

                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0        568      3371     0.16850
  TTS                0         39       222     0.17568
  exon               0         84       474     0.17722
  intron             0         31       182     0.17033
  promoter-TSS       0         50       291     0.17182


%normMotif(gc_1_dn);

                                      num_      norm_
 feature         _TYPE_    _FREQ_    motifs      size

 Intergenic         0        144       835     0.17246
 TTS                0         16       101     0.15842
 exon               0         24       129     0.18605
 intron             0         10        56     0.17857
 promoter-TSS       0         14        75     0.18667



%normMotif(cg_01_up);


                                      num_      norm_
 feature         _TYPE_    _FREQ_    motifs      size

 Intergenic         0        472      3067     0.15390
 TTS                0        244      1778     0.13723
 exon               0        443      2961     0.14961
 intron             0         69       403     0.17122
 promoter-TSS       0        269      1685     0.15964




%normMotif(cg_01_dn);

                                      num_      norm_
 feature         _TYPE_    _FREQ_    motifs      size

 Intergenic         0       1509      9729     0.15510
 TTS                0        453      2988     0.15161
 exon               0        493      3152     0.15641
 intron             0        122       762     0.16010
 promoter-TSS       0        557      4134     0.13474



%normMotif(cg_1_up);


                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0        36        354     0.10169
  TTS                0         8        105     0.07619
  exon               0        16        155     0.10323
  intron             0         2         17     0.11765
  promoter-TSS       0        26        226     0.11504


%normMotif(cg_1_dn);

                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0        141      1373     0.10269
  TTS                0         50       503     0.09940
  exon               0         48       473     0.10148
  intron             0          4        38     0.10526
  promoter-TSS       0         66       775     0.08516

%normMotif(chg_01_up);

                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0        852      5657     0.15061
  TTS                0         85       877     0.09692
  exon               0         98       981     0.09990
  intron             0         20       146     0.13699
  promoter-TSS       0         78       555     0.14054


%normMotif(chg_01_dn);

                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0       6298      44925    0.14019
  TTS                0        531       3759    0.14126
  exon               0        264       1838    0.14363
  intron             0        152       1063    0.14299
  promoter-TSS       0        816       5519    0.14785

%normMotif(chg_1_up);


                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0        427      3965     0.10769
  TTS                0         42       473     0.08879
  exon               0         27       268     0.10075
  intron             0         13       125     0.10400
  promoter-TSS       0         64       611     0.10475


%normMotif(chg_1_dn);


                                        num_      norm_
   feature         _TYPE_    _FREQ_    motifs      size

   Intergenic         0        337      3250     0.10369
   TTS                0         27       243     0.11111
   exon               0         22       240     0.09167
   intron             0          6        55     0.10909
   promoter-TSS       0         31       283     0.10954

%normMotif(chh_01_up);

                                      num_      norm_
 feature         _TYPE_    _FREQ_    motifs      size

 Intergenic         0       4748      33525    0.14163
 TTS                0        665       6263    0.10618
 exon               0        896       8430    0.10629
 intron             0        202       1570    0.12866
 promoter-TSS       0        937       6940    0.13501


%normMotif(chh_01_dn);

                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0       33687    330913    0.10180
  TTS                0        4828     52687    0.09164
  exon               0        1789     15545    0.11509
  intron             0        1633     15116    0.10803
  promoter-TSS       0        8477     95238    0.08901



%normMotif(chh_1_up);


                                       num_      norm_
  feature         _TYPE_    _FREQ_    motifs      size

  Intergenic         0       19138    186512    0.10261
  TTS                0        2602     25943    0.10030
  exon               0         962      9183    0.10476
  intron             0         647      6236    0.10375
  promoter-TSS       0        4573     45209    0.10115



%normMotif(chh_1_dn);

                                      num_      norm_
 feature         _TYPE_    _FREQ_    motifs      size

 Intergenic         0       10773     88394    0.12187
 TTS                0        1527     13378    0.11414
 exon               0         597      5399    0.11058
 intron             0         486      4139    0.11742
 promoter-TSS       0        2631     22241    0.11830



*/


