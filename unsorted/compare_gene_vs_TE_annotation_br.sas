/* Counts for arabidopsis paper */

libname brassRNA '!HOME/concannon/DTRA/brassica/sas_data';
libname brassMAP '/TB14/TB14/sandbox/dtra_sandbox/brass_wgbs';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

%macro geneTE(conditDMR,conditDAR,name);


data cg_up_&name. chg_up_&name. chh_up_&name.
     cg_dn_&name. chg_dn_&name. chh_dn_&name.;
     set brassMAP.results_by_dmr_annot;
     length feature $32.;
     feature = scan(annotation, 1, " ");
     /* FET */
     if site_type="CG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output cg_up_&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output cg_dn_&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output cg_up_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output cg_dn_&name.;
     end;


     /* FET */
     if site_type="CHG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output chg_up_&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output chg_dn_&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output chg_up_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output chg_dn_&name.;
     end;

     /* FET */
     if site_type="CHH" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output chh_up_&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output chh_dn_&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output chh_up_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output chh_dn_&name.;
     end;


     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
run;

proc sort data=cg_up_&name. nodup; by _all_; run;
proc sort data=cg_dn_&name. nodup; by _all_; run;
proc sort data=chg_up_&name. nodup; by _all_; run;
proc sort data=chg_dn_&name. nodup; by _all_; run;
proc sort data=chh_up_&name. nodup; by _all_; run;
proc sort data=chh_dn_&name. nodup; by _all_; run;




data dar_up_&name. dar_dn_&name. ;
     set brassMAP.results_by_dar_annot;
     length feature $32.;

     feature = scan(annotation, 1, " ");

     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison="&conditDAR." then output dar_up_&name.;
     if mean_methyl_diff_TRT_CTL < 0 and comparison="&conditDAR." then output dar_dn_&name.;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison="&conditDAR." then output dar_up_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison="&conditDAR." then output dar_dn_&name.;
     comparison2="&conditDMR.";
     keep comparison2 site_type chr region_Start region_stop feature distance_to_tss geneID;
     rename comparison2=comparison;
run;

proc sort data=dar_up_&name. nodup; by _all_; run;
proc sort data=dar_dn_&name. nodup; by _all_; run;



data cg_up_TE_&name. chg_up_TE_&name. chh_up_TE_&name.
     cg_dn_TE_&name. chg_dn_TE_&name. chh_dn_TE_&name.;
     set brassMAP.results_by_dmr_annot_TE;
     length feature_TE $32.;
     feature_TE = scan(annotation, 1, " ");
     /* FET */
     if site_type="CG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output cg_up_TE_&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output cg_dn_TE_&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output cg_up_TE_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output cg_dn_TE_&name.;
     end;


     /* FET */
     if site_type="CHG" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output chg_up_TE_&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output chg_dn_TE_&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output chg_up_TE_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output chg_dn_TE_&name.;
     end;

     /* FET */
     if site_type="CHH" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output chh_up_TE_&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output chh_dn_TE_&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output chh_up_TE_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output chh_dn_TE_&name.;
     end;


     keep comparison site_type chr region_Start region_stop feature_TE distance_to_tss geneID;
run;

proc sort data=cg_up_TE_&name. nodup; by _all_; run;
proc sort data=cg_dn_TE_&name. nodup; by _all_; run;
proc sort data=chg_up_TE_&name. nodup; by _all_; run;
proc sort data=chg_dn_TE_&name. nodup; by _all_; run;
proc sort data=chh_up_TE_&name. nodup; by _all_; run;
proc sort data=chh_dn_TE_&name. nodup; by _all_; run;




data dar_up_TE_&name. dar_dn_TE_&name. ;
     set brassMAP.results_by_dar_annot_TE;
     length feature_TE $32.;

     feature_TE = scan(annotation, 1, " ");

     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison="&conditDAR." then output dar_up_TE_&name.;
     if mean_methyl_diff_TRT_CTL < 0 and comparison="&conditDAR." then output dar_dn_TE_&name.;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison="&conditDAR." then output dar_up_TE_&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison="&conditDAR." then output dar_dn_TE_&name.;
     comparison2="&conditDMR.";
     keep comparison2 site_type chr region_Start region_stop feature_TE distance_to_tss geneID;
     rename comparison2=comparison;

run;

proc sort data=dar_up_TE_&name. nodup; by _all_; run;
proc sort data=dar_dn_TE_&name. nodup; by _all_; run;



%macro mergeData(inData1,inData2);
proc sort data=&inData1.;
  by comparison site_type chr region_start region_stop;
proc sort data=&inData2.;
  by comparison site_type chr region_start region_stop;
run;


data &inData1._gn_vs_te;
  merge &inData1. (in=in1) &inData2. (in=in2);
  by comparison site_type chr region_start region_stop;
  if in1 and in2;
run;

data flag_&inData1._gn_vs_te;
  set &inData1._gn_vs_te;
  if feature ne "Intergenic" then flag_has_gene=1; else flag_has_gene=0;
  if feature="exon" or feature="intron" then flag_gene_body=1; else flag_gene_body=0;

  if feature_TE ne "Intergenic" then flag_has_TE=1; else flag_has_TE=0;
  if feature_TE="exon" or feature="intron" then flag_TE_body=1; else flag_TE_body=0;
run;

proc freq data=flag_&inData1._gn_vs_te noprint;
  table flag_has_gene*flag_has_TE / out=ctab1;
  table flag_has_gene*flag_gene_body*flag_has_TE*flag_TE_body / out=ctab2;
  table flag_has_gene*flag_has_TE*feature*feature_TE / out=ctab3;
run;

proc print data=ctab1;
proc print data=ctab2;
proc print data=ctab3;
run;


%mend;

%mergeData(cg_up_&name., cg_up_TE_&name. );
%mergeData(cg_dn_&name., cg_dn_TE_&name. );
%mergeData(chg_up_&name., chg_up_TE_&name. );
%mergeData(chg_dn_&name., chg_dn_TE_&name. );
%mergeData(chh_up_&name., chh_up_TE_&name. );
%mergeData(chh_dn_&name., chh_dn_TE_&name. );
%mergeData(dar_up_&name., dar_up_TE_&name. );
%mergeData(dar_dn_&name., dar_dn_TE_&name. );

%mend;


%geneTE(0cGy_vs_10cGy_1h,10cGy_0cGy_1h,10_1);


/*

10cGy 1h CG up:


  flag_      flag_
has_gene    has_TE    feature         feature_TE      COUNT

    0          0      Intergenic      Intergenic        80
    0          1      Intergenic      TTS                5
    0          1      Intergenic      exon               9
    0          1      Intergenic      promoter-TSS       9
    1          0      TTS             Intergenic        13
    1          0      exon            Intergenic        15
    1          0      intron          Intergenic         4
    1          0      promoter-TSS    Intergenic        19
    1          1      TTS             promoter-TSS       2
    1          1      exon            TTS                1
    1          1      exon            exon               1
    1          1      intron          TTS                2
    1          1      intron          promoter-TSS       1
    1          1      promoter-TSS    TTS                2
    1          1      promoter-TSS    promoter-TSS       1

10cGy 1h CG dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic        34
     0          1      Intergenic      TTS                1
     0          1      Intergenic      exon               3
     0          1      Intergenic      promoter-TSS       8
     1          0      TTS             Intergenic         3
     1          0      exon            Intergenic         5
     1          0      intron          Intergenic         4
     1          0      promoter-TSS    Intergenic         5
     1          1      TTS             promoter-TSS       1
     1          1      exon            TTS                1
     1          1      intron          TTS                1
     1          1      promoter-TSS    TTS                1
     1          1      promoter-TSS    exon               1
     1          1      promoter-TSS    promoter-TSS       1


10cGy 1h CHG up:

  flag_      flag_
has_gene    has_TE    feature     feature_TE     COUNT

    1          0      intron     Intergenic        1
    1          1      intron     promoter-TSS      1


10cGy 1h CHG dn:


    flag_      flag_
  has_gene    has_TE    feature         feature_TE    COUNT

      0          0      Intergenic      Intergenic      5
      1          0      TTS             Intergenic      1
      1          0      exon            Intergenic      2
      1          0      promoter-TSS    Intergenic      1


10cGy 1h CHH up:


   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       2115
     0          1      Intergenic      TTS               266
     0          1      Intergenic      exon              262
     0          1      Intergenic      promoter-TSS      491
     1          0      TTS             Intergenic        461
     1          0      exon            Intergenic         51
     1          0      intron          Intergenic        331
     1          0      promoter-TSS    Intergenic        799
     1          1      TTS             TTS                47
     1          1      TTS             exon               40
     1          1      TTS             promoter-TSS       63
     1          1      exon            TTS                 2
     1          1      exon            exon                9
     1          1      exon            promoter-TSS       11
     1          1      intron          TTS                43
     1          1      intron          exon               46
     1          1      intron          promoter-TSS       63
     1          1      promoter-TSS    TTS                67
     1          1      promoter-TSS    exon               66
     1          1      promoter-TSS    promoter-TSS      126


10cGy 1h CHH dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       2873
     0          1      Intergenic      TTS               350
     0          1      Intergenic      exon              350
     0          1      Intergenic      promoter-TSS      648
     1          0      TTS             Intergenic        512
     1          0      exon            Intergenic         92
     1          0      intron          Intergenic        394
     1          0      promoter-TSS    Intergenic        930
     1          1      TTS             TTS                60
     1          1      TTS             exon               35
     1          1      TTS             promoter-TSS       67
     1          1      exon            TTS                 8
     1          1      exon            exon                7
     1          1      exon            promoter-TSS       11
     1          1      intron          TTS                41
     1          1      intron          exon               51
     1          1      intron          promoter-TSS      101
     1          1      promoter-TSS    TTS                73
     1          1      promoter-TSS    exon               95
     1          1      promoter-TSS    promoter-TSS      157


10cGy 1h DAR up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       9876
     0          1      Intergenic      TTS               997
     0          1      Intergenic      exon             1363
     0          1      Intergenic      promoter-TSS     1861
     1          0      TTS             Intergenic       2303
     1          0      exon            Intergenic       4394
     1          0      intron          Intergenic       1776
     1          0      promoter-TSS    Intergenic       3231
     1          1      TTS             TTS               193
     1          1      TTS             exon              128
     1          1      TTS             promoter-TSS      320
     1          1      exon            TTS               205
     1          1      exon            exon               96
     1          1      exon            promoter-TSS      304
     1          1      intron          TTS               132
     1          1      intron          exon              194
     1          1      intron          promoter-TSS      262
     1          1      promoter-TSS    TTS               265
     1          1      promoter-TSS    exon              205
     1          1      promoter-TSS    promoter-TSS      486

                          The SAS System

10cGy 1h DAR dn:

    flag_      flag_
  has_gene    has_TE    feature         feature_TE      COUNT

      0          0      Intergenic      Intergenic       6609
      0          1      Intergenic      TTS               626
      0          1      Intergenic      exon              952
      0          1      Intergenic      promoter-TSS     1263
      1          0      TTS             Intergenic       1282
      1          0      exon            Intergenic       1854
      1          0      intron          Intergenic       1013
      1          0      promoter-TSS    Intergenic       1745
      1          1      TTS             TTS               101
      1          1      TTS             exon              102
      1          1      TTS             promoter-TSS      199
      1          1      exon            TTS                95
      1          1      exon            exon               40
      1          1      exon            promoter-TSS      127
      1          1      intron          TTS                82
      1          1      intron          exon              153
      1          1      intron          promoter-TSS      152
      1          1      promoter-TSS    TTS               146
      1          1      promoter-TSS    exon              139
      1          1      promoter-TSS    promoter-TSS      257





*/

%geneTE(0cGy_vs_10cGy_72h,10cGy_0cGy_72h,10_72);


/*

10cGy 72h CG up:

    flag_      flag_
  has_gene    has_TE    feature         feature_TE      COUNT

      0          0      Intergenic      Intergenic       230
      0          1      Intergenic      TTS               16
      0          1      Intergenic      exon              30
      0          1      Intergenic      promoter-TSS      32
      1          0      TTS             Intergenic        33
      1          0      exon            Intergenic        87
      1          0      intron          Intergenic        28
      1          0      promoter-TSS    Intergenic        31
      1          1      TTS             TTS                2
      1          1      TTS             exon               5
      1          1      TTS             promoter-TSS       3
      1          1      exon            TTS                3
      1          1      exon            exon               2
      1          1      exon            promoter-TSS       5
      1          1      intron          TTS                5
      1          1      intron          exon               4
      1          1      intron          promoter-TSS       1
      1          1      promoter-TSS    TTS                4
      1          1      promoter-TSS    exon               4
      1          1      promoter-TSS    promoter-TSS      11


10cGy 72h CG dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       134
     0          1      Intergenic      TTS               11
     0          1      Intergenic      exon               4
     0          1      Intergenic      promoter-TSS      27
     1          0      TTS             Intergenic       114
     1          0      exon            Intergenic       563
     1          0      intron          Intergenic       173
     1          0      promoter-TSS    Intergenic        65
     1          1      TTS             TTS                5
     1          1      TTS             exon               3
     1          1      TTS             promoter-TSS       6
     1          1      exon            TTS               15
     1          1      exon            exon              10
     1          1      exon            promoter-TSS      20
     1          1      intron          TTS                5
     1          1      intron          exon               3
     1          1      intron          promoter-TSS      10
     1          1      promoter-TSS    TTS                6
     1          1      promoter-TSS    exon               3
     1          1      promoter-TSS    promoter-TSS      16


10cGy 72h CHG up:

  flag_      flag_
has_gene    has_TE    feature         feature_TE      COUNT

    0          0      Intergenic      Intergenic        59
    0          1      Intergenic      TTS                6
    0          1      Intergenic      exon              16
    0          1      Intergenic      promoter-TSS       5
    1          0      TTS             Intergenic         9
    1          0      exon            Intergenic         9
    1          0      intron          Intergenic         4
    1          0      promoter-TSS    Intergenic         5
    1          1      TTS             exon               1
    1          1      TTS             promoter-TSS       1
    1          1      exon            TTS                1
    1          1      exon            exon               2
    1          1      exon            promoter-TSS       1
    1          1      intron          TTS                1
    1          1      intron          exon               2
    1          1      intron          promoter-TSS       1
    1          1      promoter-TSS    TTS                1
    1          1      promoter-TSS    exon               1
    1          1      promoter-TSS    promoter-TSS       2

                         The SAS System

10cGy 72h CHG dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic        28
     0          1      Intergenic      TTS                4
     0          1      Intergenic      exon               2
     0          1      Intergenic      promoter-TSS       9
     1          0      TTS             Intergenic         1
     1          0      exon            Intergenic         6
     1          0      intron          Intergenic         8
     1          0      promoter-TSS    Intergenic         6
     1          1      TTS             promoter-TSS       1
     1          1      exon            TTS                2
     1          1      exon            promoter-TSS       1
     1          1      intron          TTS                1
     1          1      intron          exon               2
     1          1      intron          promoter-TSS       1
     1          1      promoter-TSS    TTS                2
     1          1      promoter-TSS    promoter-TSS       1

                          The SAS System


10cGy 72h CHH up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       3027
     0          1      Intergenic      TTS               376
     0          1      Intergenic      exon              379
     0          1      Intergenic      promoter-TSS      671
     1          0      TTS             Intergenic        612
     1          0      exon            Intergenic         99
     1          0      intron          Intergenic        473
     1          0      promoter-TSS    Intergenic       1111
     1          1      TTS             TTS                54
     1          1      TTS             exon               62
     1          1      TTS             promoter-TSS      108
     1          1      exon            TTS                 6
     1          1      exon            exon               12
     1          1      exon            promoter-TSS       12
     1          1      intron          TTS                35
     1          1      intron          exon               73
     1          1      intron          promoter-TSS      112
     1          1      promoter-TSS    TTS               104
     1          1      promoter-TSS    exon              121
     1          1      promoter-TSS    promoter-TSS      184

                          The SAS System

10cGy 72h CHH dn:

    flag_      flag_
  has_gene    has_TE    feature         feature_TE      COUNT

      0          0      Intergenic      Intergenic       4615
      0          1      Intergenic      TTS               463
      0          1      Intergenic      exon              537
      0          1      Intergenic      promoter-TSS      980
      1          0      TTS             Intergenic        805
      1          0      exon            Intergenic        208
      1          0      intron          Intergenic        635
      1          0      promoter-TSS    Intergenic       1210
      1          1      TTS             TTS                87
      1          1      TTS             exon               79
      1          1      TTS             promoter-TSS      134
      1          1      exon            TTS                29
      1          1      exon            exon               16
      1          1      exon            promoter-TSS       18
      1          1      intron          TTS                62
      1          1      intron          exon               92
      1          1      intron          promoter-TSS      135
      1          1      promoter-TSS    TTS               113
      1          1      promoter-TSS    exon              112
      1          1      promoter-TSS    promoter-TSS      216

                           The SAS System


10cGy 72h DAR up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       9365
     0          1      Intergenic      TTS               926
     0          1      Intergenic      exon             1240
     0          1      Intergenic      promoter-TSS     1727
     1          0      TTS             Intergenic       2749
     1          0      exon            Intergenic       6306
     1          0      intron          Intergenic       1983
     1          0      promoter-TSS    Intergenic       3265
     1          1      TTS             TTS               192
     1          1      TTS             exon              145
     1          1      TTS             promoter-TSS      318
     1          1      exon            TTS               283
     1          1      exon            exon              106
     1          1      exon            promoter-TSS      388
     1          1      intron          TTS               140
     1          1      intron          exon              196
     1          1      intron          promoter-TSS      273
     1          1      promoter-TSS    TTS               240
     1          1      promoter-TSS    exon              169
     1          1      promoter-TSS    promoter-TSS      450

                          The SAS System

10cGy 72h DAR dn:


   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       260
     0          1      Intergenic      TTS               26
     0          1      Intergenic      exon              47
     0          1      Intergenic      promoter-TSS      44
     1          0      TTS             Intergenic        44
     1          0      exon            Intergenic        30
     1          0      intron          Intergenic        38
     1          0      promoter-TSS    Intergenic        64
     1          1      TTS             TTS                4
     1          1      TTS             exon               4
     1          1      TTS             promoter-TSS       7
     1          1      exon            TTS                3
     1          1      exon            promoter-TSS       1
     1          1      intron          TTS                2
     1          1      intron          exon               2
     1          1      intron          promoter-TSS       5
     1          1      promoter-TSS    TTS                3
     1          1      promoter-TSS    exon               8
     1          1      promoter-TSS    promoter-TSS      14



*/

%geneTE(0cGy_vs_1p4cGy_1h,1p4cGy_0cGy_1h,14_1);


/*

1.4cGy 1h CG up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic        88
     0          1      Intergenic      TTS                5
     0          1      Intergenic      exon               4
     0          1      Intergenic      promoter-TSS       9
     1          0      TTS             Intergenic        12
     1          0      exon            Intergenic        29
     1          0      intron          Intergenic         6
     1          0      promoter-TSS    Intergenic        18
     1          1      TTS             promoter-TSS       3
     1          1      exon            TTS                1
     1          1      exon            exon               3
     1          1      exon            promoter-TSS       2
     1          1      intron          exon               2
     1          1      intron          promoter-TSS       1
     1          1      promoter-TSS    promoter-TSS       3


1.4cGy 1h CG dn:

  flag_      flag_
has_gene    has_TE    feature         feature_TE      COUNT

    0          0      Intergenic      Intergenic        26
    0          1      Intergenic      TTS                3
    0          1      Intergenic      exon               2
    0          1      Intergenic      promoter-TSS       3
    1          0      TTS             Intergenic         7
    1          0      exon            Intergenic         3
    1          0      intron          Intergenic         4
    1          0      promoter-TSS    Intergenic         3
    1          1      TTS             exon               1
    1          1      TTS             promoter-TSS       3
    1          1      intron          exon               1
    1          1      promoter-TSS    exon               1
    1          1      promoter-TSS    promoter-TSS       1

                         The SAS System


1.4cGy 1h CHG up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic        6
     0          1      Intergenic      TTS               1
     0          1      Intergenic      exon              1
     0          1      Intergenic      promoter-TSS      1
     1          0      intron          Intergenic        2
     1          0      promoter-TSS    Intergenic        1


1.4cGy 1h CHG dn:

   flag_      flag_
 has_gene    has_TE     feature      feature_TE      COUNT

     0          0      Intergenic    Intergenic        3
     0          1      Intergenic    exon              2
     0          1      Intergenic    promoter-TSS      1


1.4cGy 1h CHH up:


    flag_      flag_
  has_gene    has_TE    feature         feature_TE      COUNT

      0          0      Intergenic      Intergenic       2040
      0          1      Intergenic      TTS               275
      0          1      Intergenic      exon              249
      0          1      Intergenic      promoter-TSS      470
      1          0      TTS             Intergenic        428
      1          0      exon            Intergenic         66
      1          0      intron          Intergenic        328
      1          0      promoter-TSS    Intergenic        769
      1          1      TTS             TTS                49
      1          1      TTS             exon               34
      1          1      TTS             promoter-TSS       75
      1          1      exon            TTS                 6
      1          1      exon            exon                4
      1          1      exon            promoter-TSS        7
      1          1      intron          TTS                40
      1          1      intron          exon               35
      1          1      intron          promoter-TSS       65
      1          1      promoter-TSS    TTS                64
      1          1      promoter-TSS    exon               75
      1          1      promoter-TSS    promoter-TSS      128

1.4cGy 1h CHH dn:

  flag_      flag_
has_gene    has_TE    feature         feature_TE      COUNT

    0          0      Intergenic      Intergenic       2592
    0          1      Intergenic      TTS               300
    0          1      Intergenic      exon              315
    0          1      Intergenic      promoter-TSS      530
    1          0      TTS             Intergenic        493
    1          0      exon            Intergenic         86
    1          0      intron          Intergenic        345
    1          0      promoter-TSS    Intergenic        796
    1          1      TTS             TTS                49
    1          1      TTS             exon               40
    1          1      TTS             promoter-TSS       83
    1          1      exon            TTS                 2
    1          1      exon            exon               14
    1          1      exon            promoter-TSS       16
    1          1      intron          TTS                32
    1          1      intron          exon               54
    1          1      intron          promoter-TSS       79
    1          1      promoter-TSS    TTS                68
    1          1      promoter-TSS    exon               74
    1          1      promoter-TSS    promoter-TSS      130


1.4cGy 1h DAR up:

  flag_      flag_
has_gene    has_TE    feature         feature_TE      COUNT

    0          0      Intergenic      Intergenic      11299
    0          1      Intergenic      TTS              1177
    0          1      Intergenic      exon             1373
    0          1      Intergenic      promoter-TSS     2211
    1          0      TTS             Intergenic       4512
    1          0      exon            Intergenic      10013
    1          0      intron          Intergenic       3260
    1          0      promoter-TSS    Intergenic       6250
    1          1      TTS             TTS               315
    1          1      TTS             exon              191
    1          1      TTS             promoter-TSS      481
    1          1      exon            TTS               430
    1          1      exon            exon              147
    1          1      exon            promoter-TSS      609
    1          1      intron          TTS               214
    1          1      intron          exon              196
    1          1      intron          promoter-TSS      359
    1          1      promoter-TSS    TTS               471
    1          1      promoter-TSS    exon              239
    1          1      promoter-TSS    promoter-TSS      718


1.4cGy 1h DAR dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       753
     0          1      Intergenic      TTS               62
     0          1      Intergenic      exon             121
     0          1      Intergenic      promoter-TSS     128
     1          0      TTS             Intergenic        99
     1          0      exon            Intergenic        52
     1          0      intron          Intergenic        82
     1          0      promoter-TSS    Intergenic       174
     1          1      TTS             TTS               12
     1          1      TTS             exon              16
     1          1      TTS             promoter-TSS      21
     1          1      exon            TTS                5
     1          1      exon            exon               6
     1          1      exon            promoter-TSS       9
     1          1      intron          TTS                7
     1          1      intron          exon              18
     1          1      intron          promoter-TSS      12
     1          1      promoter-TSS    TTS               16
     1          1      promoter-TSS    exon              19
     1          1      promoter-TSS    promoter-TSS      26



*/

%geneTE(0cGy_vs_1p4cGy_72h,1p4cGy_0cGy_72h,14_72);


/*

1.4cGy 72h CG up:

    flag_      flag_
  has_gene    has_TE    feature         feature_TE      COUNT

      0          0      Intergenic      Intergenic       189
      0          1      Intergenic      TTS               15
      0          1      Intergenic      exon              32
      0          1      Intergenic      promoter-TSS      27
      1          0      TTS             Intergenic         9
      1          0      exon            Intergenic        19
      1          0      intron          Intergenic        10
      1          0      promoter-TSS    Intergenic        19
      1          1      TTS             promoter-TSS       2
      1          1      exon            exon               2
      1          1      exon            promoter-TSS       2
      1          1      intron          exon               2
      1          1      intron          promoter-TSS       3
      1          1      promoter-TSS    TTS                1
      1          1      promoter-TSS    exon               4
      1          1      promoter-TSS    promoter-TSS       4


1.4cGy 72h CG dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic        37
     0          1      Intergenic      TTS                7
     0          1      Intergenic      exon               3
     0          1      Intergenic      promoter-TSS      12
     1          0      TTS             Intergenic        17
     1          0      exon            Intergenic        34
     1          0      intron          Intergenic        26
     1          0      promoter-TSS    Intergenic        19
     1          1      TTS             exon               2
     1          1      TTS             promoter-TSS       1
     1          1      exon            TTS                3
     1          1      exon            exon               1
     1          1      exon            promoter-TSS       6
     1          1      intron          TTS                2
     1          1      intron          exon               1
     1          1      intron          promoter-TSS       3
     1          1      promoter-TSS    TTS                3
     1          1      promoter-TSS    exon               4
     1          1      promoter-TSS    promoter-TSS       6

                          The SAS System

1.4cGy 72h CHG up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       116
     0          1      Intergenic      TTS                6
     0          1      Intergenic      exon              23
     0          1      Intergenic      promoter-TSS      16
     1          0      TTS             Intergenic        18
     1          0      exon            Intergenic        14
     1          0      intron          Intergenic        14
     1          0      promoter-TSS    Intergenic         9
     1          1      TTS             TTS                1
     1          1      TTS             exon               3
     1          1      TTS             promoter-TSS       3
     1          1      exon            promoter-TSS       2
     1          1      intron          TTS                1
     1          1      intron          exon               6
     1          1      intron          promoter-TSS       2
     1          1      promoter-TSS    TTS                3
     1          1      promoter-TSS    exon               2
     1          1      promoter-TSS    promoter-TSS       5


1.4cGy 72h CHG dn:

  flag_      flag_
has_gene    has_TE    feature         feature_TE      COUNT

    0          0      Intergenic      Intergenic        5
    0          1      Intergenic      TTS               1
    0          1      Intergenic      exon              3
    0          1      Intergenic      promoter-TSS      1
    1          0      TTS             Intergenic        1
    1          0      exon            Intergenic        1
    1          0      intron          Intergenic        1
    1          0      promoter-TSS    Intergenic        5
    1          1      promoter-TSS    TTS               1


1.4cGy 72h CHH up:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       5474
     0          1      Intergenic      TTS               600
     0          1      Intergenic      exon              627
     0          1      Intergenic      promoter-TSS     1192
     1          0      TTS             Intergenic       1033
     1          0      exon            Intergenic        458
     1          0      intron          Intergenic        891
     1          0      promoter-TSS    Intergenic       1610
     1          1      TTS             TTS                88
     1          1      TTS             exon               87
     1          1      TTS             promoter-TSS      180
     1          1      exon            TTS                23
     1          1      exon            exon               18
     1          1      exon            promoter-TSS       43
     1          1      intron          TTS                87
     1          1      intron          exon               89
     1          1      intron          promoter-TSS      176
     1          1      promoter-TSS    TTS               141
     1          1      promoter-TSS    exon              146
     1          1      promoter-TSS    promoter-TSS      262

                          The SAS System


1.4cGy 72h CHH dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       1137
     0          1      Intergenic      TTS               135
     0          1      Intergenic      exon              135
     0          1      Intergenic      promoter-TSS      258
     1          0      TTS             Intergenic        235
     1          0      exon            Intergenic         42
     1          0      intron          Intergenic        149
     1          0      promoter-TSS    Intergenic        357
     1          1      TTS             TTS                25
     1          1      TTS             exon               17
     1          1      TTS             promoter-TSS       43
     1          1      exon            TTS                 2
     1          1      exon            exon                1
     1          1      exon            promoter-TSS        6
     1          1      intron          TTS                24
     1          1      intron          exon               15
     1          1      intron          promoter-TSS       40
     1          1      promoter-TSS    TTS                41
     1          1      promoter-TSS    exon               30
     1          1      promoter-TSS    promoter-TSS       59

1.4cGy 72h DAR up:

    flag_      flag_
  has_gene    has_TE    feature         feature_TE      COUNT

      0          0      Intergenic      Intergenic      10072
      0          1      Intergenic      TTS               990
      0          1      Intergenic      exon             1208
      0          1      Intergenic      promoter-TSS     1910
      1          0      TTS             Intergenic       4237
      1          0      exon            Intergenic      10921
      1          0      intron          Intergenic       3162
      1          0      promoter-TSS    Intergenic       4938
      1          1      TTS             TTS               290
      1          1      TTS             exon              142
      1          1      TTS             promoter-TSS      397
      1          1      exon            TTS               486
      1          1      exon            exon              133
      1          1      exon            promoter-TSS      624
      1          1      intron          TTS               180
      1          1      intron          exon              222
      1          1      intron          promoter-TSS      306
      1          1      promoter-TSS    TTS               373
      1          1      promoter-TSS    exon              198
      1          1      promoter-TSS    promoter-TSS      597

                           The SAS System

1.4cGy 72h DAR dn:

   flag_      flag_
 has_gene    has_TE    feature         feature_TE      COUNT

     0          0      Intergenic      Intergenic       873
     0          1      Intergenic      TTS               98
     0          1      Intergenic      exon             149
     0          1      Intergenic      promoter-TSS     146
     1          0      TTS             Intergenic       117
     1          0      exon            Intergenic       123
     1          0      intron          Intergenic       109
     1          0      promoter-TSS    Intergenic       180
     1          1      TTS             TTS               11
     1          1      TTS             exon              15
     1          1      TTS             promoter-TSS      25
     1          1      exon            TTS                8
     1          1      exon            exon               8
     1          1      exon            promoter-TSS      14
     1          1      intron          TTS               12
     1          1      intron          exon              14
     1          1      intron          promoter-TSS      17
     1          1      promoter-TSS    TTS               23
     1          1      promoter-TSS    exon              19
     1          1      promoter-TSS    promoter-TSS      25


*/



%macro compDMRDARTE(siteType, conditDMR, conditDAR, name);


data dmr_&siteType._&name.;
     set brassMAP.results_by_dmr_annot_te;
     length feature $32.;
     where site_type="&siteType.";
     feature = scan(annotation, 1, " ");
     /* FET */
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="&conditDMR." then output dmr_&siteType._&name.;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="&conditDMR." then output dmr_&siteType._&name.;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="&conditDMR." then output dmr_&siteType._&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="&conditDMR." then output dmr_&siteType._&name.;

     keep comparison site_type chr region_Start region_stop feature;
run;




data dar_&siteType._&name.;
     set brassMAP.results_by_dar_annot;
     comparison2="&conditDMR.";
     /* FET */
     if mean_methyl_diff_TRT_CTL > 0 and comparison="&conditDAR." then output dar_&siteType._&name.;
     if mean_methyl_diff_TRT_CTL < 0 and comparison="&conditDAR." then output dar_&siteType._&name.;

     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL >=0.1 and comparison="&conditDAR." then output dar_&siteType._&name.;
     if flag_binomial_FDR05=1 and mean_methyl_diff_TRT_CTL <=-0.1 and comparison="&conditDAR." then output dar_&siteType._&name.;

     keep comparison2 site_type chr region_Start region_stop ;
     rename comparison2=comparison;
run;




data stack_dar_&siteType._&name. ;
  set dmr_&siteType._&name. (in=in1) dar_&siteType._&name. (in=in2);
  length comp $25.;
  if in1 then comp="DMR_TE_&name.";
  if in2 then comp="DAR_&name.";
  keep comp site_type chr region_start region_stop feature;
run;

proc sort data=stack_dar_&siteType._&name. nodup;
  by   chr  region_start region_stop comp ;
run;


data stack_dar_&siteType._&name._super;
  retain superregion_num;
  set stack_dar_&siteType._&name.;
  by  chr region_start region_stop;
  prev_start=lag1(region_Start);
  prev_stop=lag1(region_Stop)-1;
  if first.chr then superregion_num=1;
  else do;
    if prev_stop >= region_start then superregion_num=superregion_num;
    else superregion_num = superregion_num + 1;
    end;
run;

data stack_dar_&siteType._&name._super1;
   set stack_dar_&siteType._&name._super;
   flag_present=1;
   keep flag_present superregion_num comp chr;
run;

proc sort data=stack_dar_&siteType._&name._super1 nodup;
  by chr superregion_num comp flag_present;
run;

proc transpose data=stack_dar_&siteType._&name._super1 out=stack_dar_&siteType._&name._sbys;
  by chr superregion_num;
  id comp;
  var flag_present;
run;

data stack_dar_&siteType._&name._sbys2;
  set stack_dar_&siteType._&name._sbys;
  if  DMR_TE_&name.=. then DMR_TE_&name.=0;
  if DAR_&name.=. then DAR_&name.=0;
run;

data superregion_te_annot;
  set stack_dar_&siteType._&name._super;
  if feature = "" then delete;
  if feature = "Intergenic" then delete;
  keep chr superregion_num ;
run;

proc sort data=superregion_te_annot nodup;
  by chr superregion_num;
proc sort data=stack_dar_&siteType._&name._sbys2;
  by chr superregion_num;
run;

data stack_dar_&siteType._&name._sbys3;
   merge stack_dar_&siteType._&name._sbys2 (in=in1) superregion_te_annot (in=in2);
  by chr superregion_num;
  if in2 then flag_te=1; else flag_te=0;
  if in1 then output;
run;

proc freq data=stack_dar_&siteType._&name._sbys3 noprint;
  tables flag_te*DMR_TE_&name.*DAR_&name. / out=ctabs_&siteType._dar_te_&name.;
run;

proc print data=ctabs_&siteType._dar_te_&name.;
run;


proc freq data=stack_dar_&siteType._&name._sbys3 ;
  where DMR_TE_&name.=1;
  tables flag_te*DAR_&name. / chisq;
run;

%mend;

%compDMRDARTE(CG, 0cGy_vs_10cGy_1h, 10cGy_0cGy_1h, 10_1);

/*

             DMR_TE_
  flag_te      10_1     DAR_10_1    COUNT

     0          0           1       51571
     0          1           0         135
     0          1           1          47
     1          1           0          32
     1          1           1          19


      Statistics for Table of flag_te by DAR_10_1

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      2.5638    0.1093
 Likelihood Ratio Chi-Square    1      2.4688    0.1161
 Continuity Adj. Chi-Square     1      2.0317    0.1540
 Mantel-Haenszel Chi-Square     1      2.5528    0.1101
 Phi Coefficient                       0.1049
 Contingency Coefficient               0.1043
 Cramer's V                            0.1049


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)       135
           Left-sided Pr <= F          0.9602
           Right-sided Pr >= F         0.0788

           Table Probability (P)       0.0389
           Two-sided Pr <= P           0.1166

                   Sample Size = 233


*/


%compDMRDARTE(CHG, 0cGy_vs_10cGy_1h, 10cGy_0cGy_1h, 10_1);

/*

            DMR_TE_
 flag_te      10_1     DAR_10_1    COUNT

    0          0           1       51638
    0          1           0           7
    0          1           1           3
    1          1           0           1


       Statistics for Table of flag_te by DAR_10_1

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      0.4125    0.5207
  Likelihood Ratio Chi-Square    1      0.6737    0.4118
  Continuity Adj. Chi-Square     1      0.0000    1.0000
  Mantel-Haenszel Chi-Square     1      0.3750    0.5403
  Phi Coefficient                      -0.1936
  Contingency Coefficient               0.1901
  Cramer's V                           -0.1936

   WARNING: 75% of the cells have expected counts less
            than 5. Chi-Square may not be a valid test.


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)         7
            Left-sided Pr <= F          0.7273
            Right-sided Pr >= F         1.0000

            Table Probability (P)       0.7273
            Two-sided Pr <= P           1.0000

                     Sample Size = 11
*/

%compDMRDARTE(CHH, 0cGy_vs_10cGy_1h, 10cGy_0cGy_1h, 10_1);


/*
              DMR_TE_
   flag_te      10_1     DAR_10_1    COUNT

      0          0           1       50039
      0          1           0        7416
      0          1           1        1131
      1          1           0        3206
      1          1           1         447

                  The SAS System


      Statistics for Table of flag_te by DAR_10_1

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      2.2553    0.1332
 Likelihood Ratio Chi-Square    1      2.2758    0.1314
 Continuity Adj. Chi-Square     1      2.1677    0.1409
 Mantel-Haenszel Chi-Square     1      2.2551    0.1332
 Phi Coefficient                      -0.0136
 Contingency Coefficient               0.0136
 Cramer's V                           -0.0136


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      7416
           Left-sided Pr <= F          0.0700
           Right-sided Pr >= F         0.9376

           Table Probability (P)       0.0076
           Two-sided Pr <= P           0.1408

                  Sample Size = 12200

*/

%compDMRDARTE(CG, 0cGy_vs_10cGy_72h, 10cGy_0cGy_72h, 10_72);

/*
             DMR_TE_     DAR_
  flag_te     10_72     10_72    COUNT    PERCENT

     0          0         1      35590    95.3695
     0          1         0       1342     3.5961
     0          1         1        116     0.3108
     1          1         0        246     0.6592
     1          1         1         24     0.0643


    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      0.2662    0.6059
    Likelihood Ratio Chi-Square    1      0.2601    0.6100
    Continuity Adj. Chi-Square     1      0.1557    0.6932
    Mantel-Haenszel Chi-Square     1      0.2661    0.6060
    Phi Coefficient                       0.0124
    Contingency Coefficient               0.0124
    Cramer's V                            0.0124


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      1342
              Left-sided Pr <= F          0.7425
              Right-sided Pr >= F         0.3394

              Table Probability (P)       0.0819
              Two-sided Pr <= P           0.6270

                      Sample Size = 1728

*/

%compDMRDARTE(CHG, 0cGy_vs_10cGy_72h, 10cGy_0cGy_72h, 10_72);

/*
             DMR_TE_     DAR_
  flag_te     10_72     10_72    COUNT

     0          0         1      35705
     0          1         0        119
     0          1         1         16
     1          1         0         56
     1          1         1         11



    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      0.8062    0.3692
    Likelihood Ratio Chi-Square    1      0.7842    0.3759
    Continuity Adj. Chi-Square     1      0.4601    0.4976
    Mantel-Haenszel Chi-Square     1      0.8022    0.3704
    Phi Coefficient                       0.0632
    Contingency Coefficient               0.0630
    Cramer's V                            0.0632


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)       119
              Left-sided Pr <= F          0.8674
              Right-sided Pr >= F         0.2459

              Table Probability (P)       0.1133
              Two-sided Pr <= P           0.3855

                      Sample Size = 202

*/

%compDMRDARTE(CHH, 0cGy_vs_10cGy_72h, 10cGy_0cGy_72h, 10_72);

/*
            DMR_TE_     DAR_
 flag_te     10_72     10_72    COUNT

    0          0         1      34361
    0          1         0      11838
    0          1         1        950
    1          1         0       4963
    1          1         1        415

              The SAS System

       Statistics for Table of flag_te by DAR_10_72

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      0.4512    0.5018
  Likelihood Ratio Chi-Square    1      0.4491    0.5028
  Continuity Adj. Chi-Square     1      0.4107    0.5216
  Mantel-Haenszel Chi-Square     1      0.4512    0.5018
  Phi Coefficient                       0.0050
  Contingency Coefficient               0.0050
  Cramer's V                            0.0050


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     11838
            Left-sided Pr <= F          0.7594
            Right-sided Pr >= F         0.2601

            Table Probability (P)       0.0195
            Two-sided Pr <= P           0.4978

                   Sample Size = 18166


*/




%compDMRDARTE(CG, 0cGy_vs_1p4cGy_1h, 1p4cGy_0cGy_1h, 14_1);

/*

             DMR_TE_
  flag_te      14_1     DAR_14_1    COUNT

     0          0           1       49142
     0          1           0         161
     0          1           1          35
     1          1           0          35
     1          1           1          13

                 The SAS System

     Statistics for Table of flag_te by DAR_14_1

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      2.0770    0.1495
Likelihood Ratio Chi-Square    1      1.9552    0.1620
Continuity Adj. Chi-Square     1      1.5342    0.2155
Mantel-Haenszel Chi-Square     1      2.0685    0.1504
Phi Coefficient                       0.0923
Contingency Coefficient               0.0919
Cramer's V                            0.0923


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)       161
          Left-sided Pr <= F          0.9464
          Right-sided Pr >= F         0.1096

          Table Probability (P)       0.0560
          Two-sided Pr <= P           0.1593

                  Sample Size = 244
*/

%compDMRDARTE(CHG, 0cGy_vs_1p4cGy_1h, 1p4cGy_0cGy_1h, 14_1);

/*

            DMR_TE_
 flag_te      14_1     DAR_14_1    COUNT

    0          0           1       49188
    0          1           0          10
    0          1           1           2
    1          1           0           5
    1          1           1           1

                The SAS System


          Statistics for Table of flag_te by DAR_14_1

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      0.0000    1.0000
     Likelihood Ratio Chi-Square    1      0.0000    1.0000
     Continuity Adj. Chi-Square     1      0.0000    1.0000
     Mantel-Haenszel Chi-Square     1      0.0000    1.0000
     Phi Coefficient                       0.0000
     Contingency Coefficient               0.0000
     Cramer's V                            0.0000

      WARNING: 50% of the cells have expected counts less
               than 5. Chi-Square may not be a valid test.


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)        10
               Left-sided Pr <= F          0.7549
               Right-sided Pr >= F         0.7304

               Table Probability (P)       0.4853
               Two-sided Pr <= P           1.0000

                        Sample Size = 18
*/

%compDMRDARTE(CHH, 0cGy_vs_1p4cGy_1h, 1p4cGy_0cGy_1h, 14_1);

/*
            DMR_TE_
 flag_te      14_1     DAR_14_1    COUNT

    0          0           1       48247
    0          1           0        7314
    0          1           1         623
    1          1           0        3041
    1          1           1         317


         Statistics for Table of flag_te by DAR_14_1

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      7.8270    0.0051
    Likelihood Ratio Chi-Square    1      7.6584    0.0057
    Continuity Adj. Chi-Square     1      7.6198    0.0058
    Mantel-Haenszel Chi-Square     1      7.8263    0.0051
    Phi Coefficient                       0.0263
    Contingency Coefficient               0.0263
    Cramer's V                            0.0263


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      7314
              Left-sided Pr <= F          0.9975
              Right-sided Pr >= F         0.0031

              Table Probability (P)       0.0006
              Two-sided Pr <= P           0.0058

                     Sample Size = 11295
*/


%compDMRDARTE(CG, 0cGy_vs_1p4cGy_72h, 1p4cGy_0cGy_72h, 14_72);

/*
             DMR_TE_     DAR_
  flag_te     14_72     14_72    COUNT

     0          0         1      47409
     0          1         0        330
     0          1         1         49
     1          1         0        121
     1          1         1         27

      Statistics for Table of flag_te by DAR_14_72

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      2.4358    0.1186
 Likelihood Ratio Chi-Square    1      2.3431    0.1258
 Continuity Adj. Chi-Square     1      2.0242    0.1548
 Mantel-Haenszel Chi-Square     1      2.4312    0.1189
 Phi Coefficient                       0.0680
 Contingency Coefficient               0.0678
 Cramer's V                            0.0680


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)       330
           Left-sided Pr <= F          0.9532
           Right-sided Pr >= F         0.0792

           Table Probability (P)       0.0324
           Two-sided Pr <= P           0.1297

                   Sample Size = 527
*/

%compDMRDARTE(CHG, 0cGy_vs_1p4cGy_72h, 1p4cGy_0cGy_72h, 14_72);

/*
            DMR_TE_     DAR_
 flag_te     14_72     14_72    COUNT    PERCENT

    0          0         1      47448    99.4488
    0          1         0        156     0.3270
    0          1         1         28     0.0587
    1          1         0         70     0.1467
    1          1         1          9     0.0189

              The SAS System


       Statistics for Table of flag_te by DAR_14_72

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      0.6689    0.4134
  Likelihood Ratio Chi-Square    1      0.6921    0.4055
  Continuity Adj. Chi-Square     1      0.3899    0.5324
  Mantel-Haenszel Chi-Square     1      0.6663    0.4143
  Phi Coefficient                      -0.0504
  Contingency Coefficient               0.0504
  Cramer's V                           -0.0504


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)       156
            Left-sided Pr <= F          0.2703
            Right-sided Pr >= F         0.8442

            Table Probability (P)       0.1145
            Two-sided Pr <= P           0.5618

                    Sample Size = 263
*/

%compDMRDARTE(CHH, 0cGy_vs_1p4cGy_72h, 1p4cGy_0cGy_72h, 14_72);

/*

            DMR_TE_     DAR_
 flag_te     14_72     14_72    COUNT    PERCENT

    0          0         1      46642    74.4901
    0          1         0      10800    17.2483
    0          1         1        584     0.9327
    1          1         0       4334     6.9217
    1          1         1        255     0.4073

       Statistics for Table of flag_te by DAR_14_72

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      1.1969    0.2739
  Likelihood Ratio Chi-Square    1      1.1843    0.2765
  Continuity Adj. Chi-Square     1      1.1127    0.2915
  Mantel-Haenszel Chi-Square     1      1.1968    0.2740
  Phi Coefficient                       0.0087
  Contingency Coefficient               0.0087
  Cramer's V                            0.0087


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     10800
            Left-sided Pr <= F          0.8712
            Right-sided Pr >= F         0.1459

            Table Probability (P)       0.0170
            Two-sided Pr <= P           0.2727

                   Sample Size = 15973

*/



