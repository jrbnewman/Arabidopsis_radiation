
data chh_up_dmr_01_72_te chh_up_dmr_1_72_te chh_dn_dmr_01_72_te chh_dn_dmr_1_72_te;
     set arabMAP.results_by_dmr_annot_TE;
     length feature $32.;

     feature = scan(annotation, 1, " ");
     /* FET */
     /* FET */
     if site_type="CHH" then do;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff > 0 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_te;
     if num_sites_FDR05_diff_10perc>=2 and mean_methyl_diff < 0 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_te;
     /* binomial */
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_01G" then output chh_up_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff >=0.1 and comparison="0Gy_vs_1Gy" then output chh_up_dmr_1_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_01G" then output chh_dn_dmr_01_72_te;
     if flag_binomial_FDR05=1 and mean_methyl_diff <=-0.1 and comparison="0Gy_vs_1Gy" then output chh_dn_dmr_1_72_te;
     end;


     keep comparison site_type chr region_Start region_stop feature distance_to_tss geneID;
     rename feature=feature_TE distance_to_tss=distance_to_tss_TE geneID=geneID_TE;
run;


proc sort data=chh_up_dmr_01_72_te nodup; by _all_; run;
proc sort data=chh_dn_dmr_01_72_te nodup; by _all_; run;
proc sort data=chh_up_dmr_1_72_te nodup; by _all_; run;
proc sort data=chh_dn_dmr_1_72_te nodup; by _all_; run;

proc import datafile="!HOME/concannon/useful_arabidopsis_data/TAIR10/output/TE_fragments.bed"
  out=TE_frag dmbs=tab replace;
  getnames=no;
run;

data te_frag2;
  set te_frag;
  length geneID_TE $20.;
  length TE_family $20.;
  geneID_TE=compress(scan(VAR4,2,"/"));
  TE_family=compress(scan(VAR4,1,"/"));
  keep geneID_TE TE_family;
run;


proc sort data=te_frag2;
  by geneID_TE;
proc sort data=chh_up_dmr_01_72_te;
  by geneID_TE;
proc sort data=chh_dn_dmr_01_72_te;
  by geneID_TE;
proc sort data=chh_up_dmr_1_72_te;
  by geneID_TE;
proc sort data=chh_dn_dmr_1_72_te;
  by geneID_TE;
run;

data chh_up_dmr_01_72_te2;   merge chh_up_dmr_01_72_te (in=in1) te_frag2 (in=in2); by geneID_TE; if in1 and in2; run;
data chh_dn_dmr_01_72_te2;   merge chh_dn_dmr_01_72_te (in=in1) te_frag2 (in=in2); by geneID_TE; if in1 and in2; run;
data chh_up_dmr_1_72_te2;   merge chh_up_dmr_1_72_te (in=in1) te_frag2 (in=in2); by geneID_TE; if in1 and in2; run;
data chh_dn_dmr_1_72_te2;   merge chh_dn_dmr_1_72_te (in=in1) te_frag2 (in=in2); by geneID_TE; if in1 and in2; run;


proc freq data=chh_up_dmr_01_72_te2 noprint;
tables TE_family / out=ctabs1;
run;

proc freq data=chh_dn_dmr_01_72_te2 noprint;
tables TE_family / out=ctabs2;
run;

proc freq data=chh_up_dmr_1_72_te2 noprint;
tables TE_family / out=ctabs3;
run;

proc freq data=chh_dn_dmr_1_72_te2 noprint;
tables TE_family / out=ctabs4;
run;


proc print data=ctabs1; where PERCENT > 1; run;
proc print data=ctabs2; where PERCENT > 1; run;
proc print data=ctabs3; where PERCENT > 1; run;
proc print data=ctabs4; where PERCENT > 1; run;

/*

CHH hyper DMR 10cGy:
TE_family      COUNT    PERCENT

ATDNA12T3_2      21     1.82768
ATGP1            20     1.74064
ATGP10           23     2.00174
ATHILA           18     1.56658
ATHILA0_I        19     1.65361
ATHILA2          63     5.48303
ATHILA3          28     2.43690
ATHILA4          24     2.08877
ATHILA4A         22     1.91471
ATHILA6A         39     3.39426
ATHILA6B         13     1.13142
ATLANTYS1        31     2.69800
ATLANTYS2        18     1.56658
ATLINE1_1        14     1.21845
ATLINE1_6        15     1.30548
ATLINEIII        14     1.21845
ATREP10D         21     1.82768
ATREP3           22     1.91471
ATREP4           20     1.74064
ATREP5           12     1.04439
BRODYAGA2        15     1.30548
HELITRONY3       28     2.43690
TA11             20     1.74064


CHH hypo DMR 10cGy:
 TE_family      COUNT    PERCENT

 ATDNA12T3_2      360    1.44006
 ATGP1            360    1.44006
 ATHILA           392    1.56806
 ATHILA2         1221    4.88420
 ATHILA3          541    2.16409
 ATHILA4          377    1.50806
 ATHILA4A         448    1.79207
 ATHILA4C         363    1.45206
 ATHILA6A         635    2.54010
 ATHILA6B         307    1.22805
 ATLANTYS1        286    1.14405
 ATLANTYS2        258    1.03204
 ATLINE1A         280    1.12004
 ATREP10D         548    2.19209
 ATREP11          263    1.05204
 ATREP15          340    1.36005
 ATREP3           667    2.66811
 ATREP5           252    1.00804
 BRODYAGA2        281    1.12404
 HELITRONY1D      267    1.06804
 HELITRONY3       778    3.11212
 TA11             429    1.71607
 VANDAL3          271    1.08404



CHH hyper DMR 100cGy:
 TE_family      COUNT    PERCENT

 ARNOLDY1        120     1.16088
 ARNOLDY2        152     1.47045
 ATDNA12T3_2     104     1.00609
 ATGP1           267     2.58295
 ATHILA          179     1.73164
 ATHILA2         553     5.34971
 ATHILA3         134     1.29631
 ATHILA4         147     1.42208
 ATHILA4A        131     1.26729
 ATHILA4C        112     1.08349
 ATHILA6A        248     2.39915
 ATLANTYS1       125     1.20925
 ATLANTYS2       118     1.14153
 ATLINE1A        109     1.05446
 ATREP10D        181     1.75099
 ATREP15         116     1.12218
 ATREP18         119     1.15120
 ATREP3          212     2.05089
 BRODYAGA2       114     1.10283
 HELITRONY1D     106     1.02544
 HELITRONY3      340     3.28916
 TA11            170     1.64458
 VANDAL3         191     1.84773


CHH hypo DMR 100cGy:

 TE_family      COUNT    PERCENT

 ATDNA12T3_2      60     1.70843
 ATGP1            50     1.42369
 ATHILA           48     1.36674
 ATHILA0_I        46     1.30979
 ATHILA2         188     5.35308
 ATHILA3          91     2.59112
 ATHILA4          46     1.30979
 ATHILA4A         60     1.70843
 ATHILA4C         72     2.05011
 ATHILA6A         94     2.67654
 ATHILA6B         44     1.25285
 ATLANTYS1        41     1.16743
 ATLINE1A         47     1.33827
 ATREP10D         94     2.67654
 ATREP15          46     1.30979
 ATREP3           87     2.47722
 ATREP5           41     1.16743
 BRODYAGA2        45     1.28132
 HELITRONY3      126     3.58770
 TA11             39     1.11048




*/



