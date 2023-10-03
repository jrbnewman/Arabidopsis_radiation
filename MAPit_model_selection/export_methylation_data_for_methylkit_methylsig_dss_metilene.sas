libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';
libname wgbslocA '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';

*libname wgbslocA '/blue/concannon/share/jnewman/mingqi_arab/sas_data';
ods listing;
ods html close;
options nosyntaxcheck;



proc datasets lib=work kill noprint;
run;
quit;


/* Prep count matrices for methylKit, dss, metilene, DSS */

%macro prepData(trt);

data counts;
   set wgbslocA.methylation_data_cg_chg_chh;
   if chr = "Mt" then delete;
   if chr = "Pt" then delete;
   where units="0U" and treatment="&trt.";
run;  

data cov;
   set wgbslocA.flag_coverage_10x_cg_chg_chh;
   where flag_&trt._coverage_ge_10x=1;
   keep site_type chr start_pos stop_pos;
run;

proc sort data=cov;
  by site_type chr start_pos stop_pos;
proc sort data=counts;
  by site_type chr start_pos stop_pos;
run;

data counts_&trt.;
   merge counts (in=in1) cov (in=in2);
   by site_type chr start_pos stop_pos;
   if in1 and in2;
run;

%mend;

%prepData(0Gy);
%prepData(01Gy);
%prepData(1Gy);

/* MethylKit format
for each rep:
(1) chr.pos 
(2) chr
(3) position
(4) strand (F/R)
(5) Coverage
(6) Freq methyl
(7) Freq non-methyl
*/

data methylkit_0gy_R1_CG methylkit_0gy_R2_CG
     methylkit_0gy_R1_CHG methylkit_0gy_R2_CHG
     methylkit_0gy_R1_CHH methylkit_0gy_R2_CHH;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_0gy;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(methyl_C / total_C * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if site_type="CG" and rep = 1 then output methylkit_0gy_R1_CG;
   if site_type="CG" and rep = 2 then output methylkit_0gy_R2_CG;
   if site_type="CHG" and rep = 1 then output methylkit_0gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output methylkit_0gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output methylkit_0gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output methylkit_0gy_R2_CHH;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;


data methylkit_01gy_R1_CG methylkit_01gy_R2_CG
     methylkit_01gy_R1_CHG methylkit_01gy_R2_CHG
     methylkit_01gy_R1_CHH methylkit_01gy_R2_CHH;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_01gy;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(methyl_C / total_C * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if site_type="CG" and rep = 1 then output methylkit_01gy_R1_CG;
   if site_type="CG" and rep = 2 then output methylkit_01gy_R2_CG;
   if site_type="CHG" and rep = 1 then output methylkit_01gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output methylkit_01gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output methylkit_01gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output methylkit_01gy_R2_CHH;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;



data methylkit_1gy_R1_CG methylkit_1gy_R2_CG
     methylkit_1gy_R1_CHG methylkit_1gy_R2_CHG
     methylkit_1gy_R1_CHH methylkit_1gy_R2_CHH;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_01gy;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(methyl_C / total_C * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if site_type="CG" and rep = 1 then output methylkit_1gy_R1_CG;
   if site_type="CG" and rep = 2 then output methylkit_1gy_R2_CG;
   if site_type="CHG" and rep = 1 then output methylkit_1gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output methylkit_1gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output methylkit_1gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output methylkit_1gy_R2_CHH;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;

proc export data=methylkit_0gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r1_CG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_0gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r2_CG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_01gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r1_CG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_01gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r2_CG.txt"
     dbms=tab replace;
run;
proc export data=methylkit_1gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r1_CG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_1gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r2_CG.txt"
     dbms=tab replace;
run;


proc export data=methylkit_0gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r1_CHG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_0gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r2_CHG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_01gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r1_CHG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_01gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r2_CHG.txt"
     dbms=tab replace;
run;
proc export data=methylkit_1gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r1_CHG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_1gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r2_CHG.txt"
     dbms=tab replace;
run;

proc export data=methylkit_0gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r1_CHH.txt"
     dbms=tab replace;
run;

proc export data=methylkit_0gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r2_CHH.txt"
     dbms=tab replace;
run;

proc export data=methylkit_01gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r1_CHH.txt"
     dbms=tab replace;
run;

proc export data=methylkit_01gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r2_CHH.txt"
     dbms=tab replace;
run;
proc export data=methylkit_1gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r1_CHH.txt"
     dbms=tab replace;
run;

proc export data=methylkit_1gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r2_CHH.txt"
     dbms=tab replace;
run;


/* methylSig format
The coverage output looks like this (tab-delimited; 1-based genomic coords):
<chromosome>
<start position>
<end position>
<methylation percentage>
<count methylated>
<count unmethylated>


 */


data methylsig_0gy_R1_CG methylsig_0gy_R2_CG
     methylsig_0gy_R1_CHG methylsig_0gy_R2_CHG
     methylsig_0gy_R1_CHH methylsig_0gy_R2_CHH;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_0gy;
   perc_methyl = round(methyl_C / total_C * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if site_type="CG" and rep = 1 then output methylsig_0gy_R1_CG;
   if site_type="CG" and rep = 2 then output methylsig_0gy_R2_CG;
   if site_type="CHG" and rep = 1 then output methylsig_0gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output methylsig_0gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output methylsig_0gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output methylsig_0gy_R2_CHH;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;



data methylsig_01gy_R1_CG methylsig_01gy_R2_CG
     methylsig_01gy_R1_CHG methylsig_01gy_R2_CHG
     methylsig_01gy_R1_CHH methylsig_01gy_R2_CHH;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_01gy;
   perc_methyl = round(methyl_C / total_C * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if site_type="CG" and rep = 1 then output methylsig_01gy_R1_CG;
   if site_type="CG" and rep = 2 then output methylsig_01gy_R2_CG;
   if site_type="CHG" and rep = 1 then output methylsig_01gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output methylsig_01gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output methylsig_01gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output methylsig_01gy_R2_CHH;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;


data methylsig_1gy_R1_CG methylsig_1gy_R2_CG
     methylsig_1gy_R1_CHG methylsig_1gy_R2_CHG
     methylsig_1gy_R1_CHH methylsig_1gy_R2_CHH;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_1gy;
   perc_methyl = round(methyl_C / total_C * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if site_type="CG" and rep = 1 then output methylsig_1gy_R1_CG;
   if site_type="CG" and rep = 2 then output methylsig_1gy_R2_CG;
   if site_type="CHG" and rep = 1 then output methylsig_1gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output methylsig_1gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output methylsig_1gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output methylsig_1gy_R2_CHH;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;


proc export data=methylsig_0gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r1_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_0gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r2_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_01gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r1_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_01gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r2_CG.txt"
     dbms=tab replace;
     putnames=no;
run;
proc export data=methylsig_1gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r1_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_1gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r2_CG.txt"
     dbms=tab replace;
     putnames=no;
run;


proc export data=methylsig_0gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r1_CHG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_0gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r2_CHG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_01gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r1_CHG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_01gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r2_CHG.txt"
     dbms=tab replace;
     putnames=no;
run;
proc export data=methylsig_1gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r1_CHG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_1gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r2_CHG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_0gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r1_CHH.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_0gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r2_CHH.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_01gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r1_CHH.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_01gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r2_CHH.txt"
     dbms=tab replace;
     putnames=no;
run;
proc export data=methylsig_1gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r1_CHH.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=methylsig_1gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r2_CHH.txt"
     dbms=tab replace;
     putnames=no;
run;





data dss_0gy_R1_CG dss_0gy_R2_CG
     dss_0gy_R1_CHG dss_0gy_R2_CHG
     dss_0gy_R1_CHH dss_0gy_R2_CHH;
   retain chr stop_pos total_C methyl_C;
   set counts_0gy;
   if site_type="CG" and rep = 1 then output dss_0gy_R1_CG;
   if site_type="CG" and rep = 2 then output dss_0gy_R2_CG;
   if site_type="CHG" and rep = 1 then output dss_0gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output dss_0gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output dss_0gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output dss_0gy_R2_CHH;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;



data dss_01gy_R1_CG dss_01gy_R2_CG
     dss_01gy_R1_CHG dss_01gy_R2_CHG
     dss_01gy_R1_CHH dss_01gy_R2_CHH;
   retain  chr stop_pos total_C methyl_C;
   set counts_01gy;
   if site_type="CG" and rep = 1 then output dss_01gy_R1_CG;
   if site_type="CG" and rep = 2 then output dss_01gy_R2_CG;
   if site_type="CHG" and rep = 1 then output dss_01gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output dss_01gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output dss_01gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output dss_01gy_R2_CHH;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;


data dss_1gy_R1_CG dss_1gy_R2_CG
     dss_1gy_R1_CHG dss_1gy_R2_CHG
     dss_1gy_R1_CHH dss_1gy_R2_CHH;
   retain  chr stop_pos total_C methyl_C;
   set counts_1gy;
   if site_type="CG" and rep = 1 then output dss_1gy_R1_CG;
   if site_type="CG" and rep = 2 then output dss_1gy_R2_CG;
   if site_type="CHG" and rep = 1 then output dss_1gy_R1_CHG;
   if site_type="CHG" and rep = 2 then output dss_1gy_R2_CHG;
   if site_type="CHH" and rep = 1 then output dss_1gy_R1_CHH;
   if site_type="CHH" and rep = 2 then output dss_1gy_R2_CHH;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;


proc export data=dss_0gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r1_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=dss_0gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r2_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=dss_01gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_r1_CG.txt"
     dbms=tab replace;
     putnames=no;
run;

proc export data=dss_01gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_r2_CG.txt"
     dbms=tab replace;
run;
proc export data=dss_1gy_r1_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r1_CG.txt"
     dbms=tab replace;
run;

proc export data=dss_1gy_r2_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r2_CG.txt"
     dbms=tab replace;
run;


proc export data=dss_0gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r1_CHG.txt"
     dbms=tab replace;
run;

proc export data=dss_0gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r2_CHG.txt"
     dbms=tab replace;
run;

proc export data=dss_01gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_r1_CHG.txt"
     dbms=tab replace;
run;

proc export data=dss_01gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_r2_CHG.txt"
     dbms=tab replace;
run;
proc export data=dss_1gy_r1_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r1_CHG.txt"
     dbms=tab replace;
run;

proc export data=dss_1gy_r2_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r2_CHG.txt"
     dbms=tab replace;
run;

proc export data=dss_0gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r1_CHH.txt"
     dbms=tab replace;
run;

proc export data=dss_0gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r2_CHH.txt"
     dbms=tab replace;
run;

proc export data=dss_01gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_r1_CHH.txt"
     dbms=tab replace;
run;

proc export data=dss_01gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_r2_CHH.txt"
     dbms=tab replace;
run;
proc export data=dss_1gy_r1_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r1_CHH.txt"
     dbms=tab replace;
run;

proc export data=dss_1gy_r2_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r2_CHH.txt"
     dbms=tab replace;
run;

/* Format for metiline:
chr
pos
g1_r1 methyl
g1_r2 methyl
g2_r1 methyl
g2_r2 methyl
*/


data met_0gy;
  set counts_0Gy;
  keep site_Type chr stop_pos perc_methyl rep;
  rename stop_pos=pos;
run;

proc sort data=met_0gy;
  by site_type chr pos rep;
proc transpose data=met_0gy out=met_0gy_sbys(rename=(_1=C0gy_R1 _2=C0gy_R2) drop=_NAME_);
  by site_type chr pos;
  id rep;
  var perc_methyl;
run;


data met_01gy;
  set counts_01Gy;
  keep site_Type chr stop_pos perc_methyl rep;
  rename stop_pos=pos;
run;

proc sort data=met_01gy;
  by site_type chr pos rep;
proc transpose data=met_01gy out=met_01gy_sbys(rename=(_1=C01gy_R1 _2=C01gy_R2) drop=_NAME_);
  by site_type chr pos;
  id rep;
  var perc_methyl;
run;


data met_1gy;
  set counts_1Gy;
  keep site_Type chr stop_pos perc_methyl rep;
  rename stop_pos=pos;
run;

proc sort data=met_1gy;
  by site_type chr pos rep;
proc transpose data=met_1gy out=met_1gy_sbys(rename=(_1=C1gy_R1 _2=C1gy_R2) drop=_NAME_);
  by site_type chr pos;
  id rep;
  var perc_methyl;
run;

proc sort data=met_0gy_sbys;
  by site_type chr pos;
proc sort data=met_01gy_sbys;
  by site_type chr pos;
proc sort data=met_1gy_sbys;
  by site_type chr pos;
run;

data met_0gy_vs_01gy;
  merge met_0gy_sbys (in=in1) met_01gy_sbys (in=in2);
  by site_type chr pos;
  if in1 and in2;
run;

data met_0gy_vs_1gy;
  merge met_0gy_sbys (in=in1) met_1gy_sbys (in=in2);
  by site_type chr pos;
  if in1 and in2;
run;

data met_0gy_vs_01gy_CG met_0gy_vs_01gy_CHG met_0gy_vs_01gy_CHH;
  set met_0gy_vs_01gy;
  if site_type = "CG" then output met_0gy_vs_01gy_CG;
   if site_type = "CHG" then output met_0gy_vs_01gy_CHG;
   if site_type = "CHH" then output met_0gy_vs_01gy_CHH;
  drop site_type;
 run;



data met_0gy_vs_1gy_CG met_0gy_vs_1gy_CHG met_0gy_vs_1gy_CHH;
  set met_0gy_vs_1gy;
  if site_type = "CG" then output met_0gy_vs_1gy_CG;
   if site_type = "CHG" then output met_0gy_vs_1gy_CHG;
   if site_type = "CHH" then output met_0gy_vs_1gy_CHH;
  drop site_type;
 run;


proc export data=met_0gy_vs_01gy_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_CG.txt"
     dbms=tab replace;
run;


proc export data=met_0gy_vs_01gy_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_CHG.txt"
     dbms=tab replace;
run;


proc export data=met_0gy_vs_01gy_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_CHH.txt"
     dbms=tab replace;
run;

proc export data=met_0gy_vs_1gy_CG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_CG.txt"
     dbms=tab replace;
run;


proc export data=met_0gy_vs_1gy_CHG
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_CHG.txt"
     dbms=tab replace;
run;


proc export data=met_0gy_vs_1gy_CHH
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_CHH.txt"
     dbms=tab replace;
run;



/***********  GC DATA ******************/



%macro prepData(trt,units);

data counts;
set wgbslocA.methylation_data_gc;
if chr = "Mt" then delete;
if chr = "Pt" then delete;
where units="&units." and treatment="&trt.";
run;

data cov;
set wgbslocA.flag_coverage_10x_gc;
where flag_&trt._&units._coverage_ge_10x=1;
keep site_type chr start_pos stop_pos;
run;

proc sort data=cov;
by site_type chr start_pos stop_pos;
proc sort data=counts;
by site_type chr start_pos stop_pos;
run;

data counts_&trt._&units.;
merge counts (in=in1) cov (in=in2);
by site_type chr start_pos stop_pos;
if in1 and in2;
run;

%mend;


%prepData(0Gy,0U);
%prepData(01Gy,0U);
%prepData(1Gy,0U);
%prepData(0Gy,100U);
%prepData(01Gy,100U);
%prepData(1Gy,100U);


/* MethylKit format
for each rep:
(1) chr.pos 
(2) chr
(3) position
(4) strand (F/R)
(5) Coverage
(6) Freq methyl
(7) Freq non-methyl
*/

data methylkit_0gy_0U_R1 methylkit_0gy_0U_R2;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_0gy_0U;
   where perc_methyl_norm ne .;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(perc_methyl_norm * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if rep = 1 then output methylkit_0gy_0U_R1;
   if rep = 2 then output methylkit_0gy_0U_R2;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;

data methylkit_0gy_100U_R1 methylkit_0gy_100U_R2;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_0gy_100U;
   where perc_methyl_norm ne .;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(perc_methyl_norm * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if rep = 1 then output methylkit_0gy_100U_R1;
   if rep = 2 then output methylkit_0gy_100U_R2;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;


data methylkit_01gy_0U_R1 methylkit_01gy_0U_R2;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_01gy_0U;
   where perc_methyl_norm ne .;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(perc_methyl_norm * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if rep = 1 then output methylkit_01gy_0U_R1;
   if rep = 2 then output methylkit_01gy_0U_R2;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;

data methylkit_01gy_100U_R1 methylkit_01gy_100U_R2;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_01gy_100U;
   where perc_methyl_norm ne .;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(methyl_C / total_C * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if rep = 1 then output methylkit_01gy_100U_R1;
   if rep = 2 then output methylkit_01gy_100U_R2;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;



data methylkit_1gy_0U_R1 methylkit_1gy_0U_R2;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_1gy_0U;
   where perc_methyl_norm ne .;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(methyl_C / total_C * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if rep = 1 then output methylkit_1gy_0U_R1;
   if rep = 2 then output methylkit_1gy_0U_R2;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;

data methylkit_1gy_100U_R1 methylkit_1gy_100U_R2;
   retain chrBase chr stop_pos strand total_C freqC freqT;
   length chrBase $20.;
   length strand $1.;
   format freqC 6.2;
   format freqT 6.2;
   set counts_1gy_100U;
   where perc_methyl_norm ne .;
   chrBase = catx(".",chr,stop_pos);
   strand = "F";
   freqC = round(methyl_C / total_C * 100, 0.01);
   freqT = round(100 - freqC , 0.01);
   if rep = 1 then output methylkit_1gy_100U_R1;
   if rep = 2 then output methylkit_1gy_100U_R2;
   keep chrBase chr stop_pos strand total_C freqC freqT;
   rename stop_pos=base total_c=coverage;
run;











proc export data=methylkit_0gy_0U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_0U_R1.txt" dbms=tab replace; run;
proc export data=methylkit_0gy_0U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_0U_R2.txt" dbms=tab replace; run;
proc export data=methylkit_0gy_100U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_100U_R1.txt" dbms=tab replace; run;
proc export data=methylkit_0gy_100U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_100U_R2.txt" dbms=tab replace; run;

proc export data=methylkit_01gy_0U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_0U_R1.txt" dbms=tab replace; run;
proc export data=methylkit_01gy_0U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_0U_R2.txt" dbms=tab replace; run;
proc export data=methylkit_01gy_100U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_100U_R1.txt" dbms=tab replace; run;
proc export data=methylkit_01gy_100U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_100U_R2.txt" dbms=tab replace; run;

proc export data=methylkit_1gy_0U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_0U_R1.txt" dbms=tab replace; run;
proc export data=methylkit_1gy_0U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_0U_R2.txt" dbms=tab replace; run;
proc export data=methylkit_1gy_100U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_100U_R1.txt" dbms=tab replace; run;
proc export data=methylkit_1gy_100U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_100U_R2.txt" dbms=tab replace; run;



/* methylSig format
The coverage output looks like this (tab-delimited; 1-based genomic coords):
<chromosome>
<start position>
<end position>
<methylation percentage>
<count methylated>
<count unmethylated>


 */


data methylsig_0gy_0U_R1 methylsig_0gy_0U_R2;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_0gy_0U;
   where perc_methyl_norm ne .;
   perc_methyl = round(perc_methyl_norm * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if rep = 1 then output methylsig_0gy_0U_R1;
   if rep = 2 then output methylsig_0gy_0U_R2;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;



data methylsig_0gy_100U_R1 methylsig_0gy_100U_R2;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_0gy_100U;
   where perc_methyl_norm ne .;
   perc_methyl = round(perc_methyl_norm * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if rep = 1 then output methylsig_0gy_100U_R1;
   if rep = 2 then output methylsig_0gy_100U_R2;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;





data methylsig_01gy_0U_R1 methylsig_01gy_0U_R2;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_01gy_0U;
   where perc_methyl_norm ne .;
   perc_methyl = round(perc_methyl_norm * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if rep = 1 then output methylsig_01gy_0U_R1;
   if rep = 2 then output methylsig_01gy_0U_R2;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;



data methylsig_01gy_100U_R1 methylsig_01gy_100U_R2;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_01gy_100U;
   where perc_methyl_norm ne .;
   perc_methyl = round(perc_methyl_norm * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if rep = 1 then output methylsig_01gy_100U_R1;
   if rep = 2 then output methylsig_01gy_100U_R2;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;



data methylsig_1gy_0U_R1 methylsig_1gy_0U_R2;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_1gy_0U;
   where perc_methyl_norm ne .;
   perc_methyl = round(perc_methyl_norm * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if rep = 1 then output methylsig_1gy_0U_R1;
   if rep = 2 then output methylsig_1gy_0U_R2;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;



data methylsig_1gy_100U_R1 methylsig_1gy_100U_R2;
   retain chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
   format perc_methyl 3.0;
   set counts_1gy_100U;
   where perc_methyl_norm ne .;
   perc_methyl = round(perc_methyl_norm * 100, 1);
   unmethyl_C = total_C - methyl_C;
   if rep = 1 then output methylsig_1gy_100U_R1;
   if rep = 2 then output methylsig_1gy_100U_R2;
   keep chr start_pos stop_pos perc_methyl methyl_C unmethyl_C;
run;








proc export data=methylsig_0gy_0U_R1 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_0U_R1.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_0gy_0U_R2 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_0U_R2.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_0gy_100U_R1 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_100U_R1.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_0gy_100U_R2 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_100U_R2.txt"  dbms=tab replace; putnames=no; run;

proc export data=methylsig_01gy_0U_R1 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_0U_R1.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_01gy_0U_R2 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_0U_R2.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_01gy_100U_R1 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_100U_R1.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_01gy_100U_R2 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_100U_R2.txt"  dbms=tab replace; putnames=no; run;

proc export data=methylsig_1gy_0U_R1 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_0U_R1.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_1gy_0U_R2 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_0U_R2.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_1gy_100U_R1 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_100U_R1.txt"  dbms=tab replace; putnames=no; run;
proc export data=methylsig_1gy_100U_R2 outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_100U_R2.txt"  dbms=tab replace; putnames=no; run;





data dss_0gy_0U_R1 dss_0gy_0U_R2;
   retain chr stop_pos total_C methyl_C;
   set counts_0gy_0U;
   if rep = 1 then output dss_0gy_0U_R1;
   if rep = 2 then output dss_0gy_0U_R2;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;
data dss_0gy_100U_R1 dss_0gy_100U_R2;
   retain chr stop_pos total_C methyl_C;
   set counts_0gy_0U;
   if rep = 1 then output dss_0gy_100U_R1;
   if rep = 2 then output dss_0gy_100U_R2;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;

data dss_01gy_0U_R1 dss_01gy_0U_R2;
   retain chr stop_pos total_C methyl_C;
   set counts_01gy_0U;
   if rep = 1 then output dss_01gy_0U_R1;
   if rep = 2 then output dss_01gy_0U_R2;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;

data dss_01gy_100U_R1 dss_01gy_100U_R2;
   retain chr stop_pos total_C methyl_C;
   set counts_01gy_100U;
   if rep = 1 then output dss_01gy_100U_R1;
   if rep = 2 then output dss_01gy_100U_R2;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;

data dss_1gy_0U_R1 dss_1gy_0U_R2;
   retain chr stop_pos total_C methyl_C;
   set counts_1gy_0U;
   if rep = 1 then output dss_1gy_0U_R1;
   if rep = 2 then output dss_1gy_0U_R2;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;
data dss_1gy_100U_R1 dss_1gy_100U_R2;
   retain chr stop_pos total_C methyl_C;
   set counts_1gy_0U;
   if rep = 1 then output dss_1gy_100U_R1;
   if rep = 2 then output dss_1gy_100U_R2;
   keep chr stop_pos total_C methyl_C;
   rename stop_pos = pos total_C = N methyl_C = X;
run;




proc export data=dss_0gy_0U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_0U_R1.txt"  dbms=tab replace;  run;
proc export data=dss_0gy_0U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_0U_R2.txt"  dbms=tab replace;  run;
proc export data=dss_0gy_100U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_100U_R1.txt"  dbms=tab replace; run;
proc export data=dss_0gy_100U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_100U_R2.txt"  dbms=tab replace; run;

proc export data=dss_01gy_0U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_0U_R1.txt"  dbms=tab replace;   run;
proc export data=dss_01gy_0U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_0U_R2.txt"  dbms=tab replace; run;
proc export data=dss_01gy_100U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_100U_R1.txt"  dbms=tab replace;   run;
proc export data=dss_01gy_100U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_100U_R2.txt"  dbms=tab replace;  run;

proc export data=dss_1gy_0U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_0U_R1.txt"  dbms=tab replace;  run;
proc export data=dss_1gy_0U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_0U_R2.txt"  dbms=tab replace; run;
proc export data=dss_1gy_100U_R1   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_100U_R1.txt"  dbms=tab replace;  run;
proc export data=dss_1gy_100U_R2   outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_100U_R2.txt"  dbms=tab replace;   run;



/* Format for metiline:
chr
pos
g1_r1 methyl
g1_r2 methyl
g2_r1 methyl
g2_r2 methyl
*/


data met_0gy_0U;
  set counts_0Gy_0U;
  keep chr stop_pos perc_methyl_norm rep;
  rename stop_pos=pos;
run;

proc sort data=met_0gy_0U;
  by chr pos rep;
proc transpose data=met_0gy_0U out=met_0gy_0U_sbys(rename=(_1=C0gy_0U_R1 _2=C0gy_0U_R2) drop=_NAME_);
  by  chr pos;
  id rep;
  var perc_methyl_norm;
run;


data met_0gy_100U;
  set counts_0Gy_100U;
  keep chr stop_pos perc_methyl_norm rep;
  rename stop_pos=pos;
run;

proc sort data=met_0gy_100U;
  by chr pos rep;
proc transpose data=met_0gy_100U out=met_0gy_100U_sbys(rename=(_1=C0gy_100U_R1 _2=C0gy_100U_R2) drop=_NAME_);
  by  chr pos;
  id rep;
  var perc_methyl_norm;
run;




data met_01gy_0U;
  set counts_01Gy_0U;
  keep chr stop_pos perc_methyl_norm rep;
  rename stop_pos=pos;
run;

proc sort data=met_01gy_0U;
  by chr pos rep;
proc transpose data=met_01gy_0U out=met_01gy_0U_sbys(rename=(_1=C01gy_0U_R1 _2=C01gy_0U_R2) drop=_NAME_);
  by  chr pos;
  id rep;
  var perc_methyl_norm;
run;


data met_01gy_100U;
  set counts_01Gy_100U;
  keep chr stop_pos perc_methyl_norm rep;
  rename stop_pos=pos;
run;

proc sort data=met_01gy_100U;
  by chr pos rep;
proc transpose data=met_01gy_100U out=met_01gy_100U_sbys(rename=(_1=C01gy_100U_R1 _2=C01gy_100U_R2) drop=_NAME_);
  by  chr pos;
  id rep;
  var perc_methyl_norm;
run;






data met_1gy_0U;
  set counts_1Gy_0U;
  keep chr stop_pos perc_methyl_norm rep;
  rename stop_pos=pos;
run;

proc sort data=met_1gy_0U;
  by chr pos rep;
proc transpose data=met_1gy_0U out=met_1gy_0U_sbys(rename=(_1=C1gy_0U_R1 _2=C1gy_0U_R2) drop=_NAME_);
  by  chr pos;
  id rep;
  var perc_methyl_norm;
run;


data met_1gy_100U;
  set counts_1Gy_100U;
  keep chr stop_pos perc_methyl_norm rep;
  rename stop_pos=pos;
run;

proc sort data=met_1gy_100U;
  by chr pos rep;
proc transpose data=met_1gy_100U out=met_1gy_100U_sbys(rename=(_1=C1gy_100U_R1 _2=C1gy_100U_R2) drop=_NAME_);
  by  chr pos;
  id rep;
  var perc_methyl_norm;
run;




proc sort data=met_0gy_0U_sbys;
  by  chr pos;
proc sort data=met_01gy_0U_sbys;
  by  chr pos;
proc sort data=met_1gy_0U_sbys;
  by  chr pos;
proc sort data=met_0gy_100U_sbys;
  by  chr pos;
proc sort data=met_01gy_100U_sbys;
  by  chr pos;
proc sort data=met_1gy_100U_sbys;
  by  chr pos;
run;

data met_0gy_vs_01gy;
  merge met_0gy_0U_sbys (in=in1) met_0gy_100U_sbys (in=in2) met_01gy_0U_sbys (in=in3) met_01gy_100U_sbys (in=in4);
  by  chr pos;
  if in1 and in2 and in3 and in4;
run;


data met_0gy_vs_1gy;
  merge met_0gy_0U_sbys (in=in1) met_0gy_100U_sbys (in=in2) met_1gy_0U_sbys (in=in3) met_1gy_100U_sbys (in=in4);
  by  chr pos;
  if in1 and in2 and in3 and in4;
run;

data met_0gy_vs_01gy_2;
   set met_0gy_vs_01gy;
   C0gy_R1 = C0gy_100U_R1 - C0gy_0U_R1;
   C0gy_R2 = C0gy_100U_R2 - C0gy_0U_R2;
   C01gy_R1 = C01gy_100U_R1 - C01gy_0U_R1;
   C01gy_R2 = C01gy_100U_R2 - C01gy_0U_R2;
   drop  C0gy_100U_R1 C0gy_100U_R2 C0gy_0U_R1 C0gy_0U_R2
         C01gy_100U_R1 C01gy_100U_R2 C01gy_0U_R1 C01gy_0U_R2;
run;


data met_0gy_vs_1gy_2;
   set met_0gy_vs_1gy;
   C0gy_R1 = C0gy_100U_R1 - C0gy_0U_R1;
   C0gy_R2 = C0gy_100U_R2 - C0gy_0U_R2;
   C1gy_R1 = C1gy_100U_R1 - C1gy_0U_R1;
   C1gy_R2 = C1gy_100U_R2 - C1gy_0U_R2;
   drop  C0gy_100U_R1 C0gy_100U_R2 C0gy_0U_R1 C0gy_0U_R2
         C1gy_100U_R1 C1gy_100U_R2 C1gy_0U_R1 C1gy_0U_R2;
run;



proc export data=met_0gy_vs_01gy_2
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_GC.txt"
     dbms=tab replace;
run;


proc export data=met_0gy_vs_1gy_2
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_GC.txt"
     dbms=tab replace;
run;




data met_0gy_vs_01gy_2a;
   set met_0gy_vs_01gy_2;
   if C0gy_R1 < 0 or C0gy_R2 < 0  or  C01gy_R1 < 0 or C01gy_R2 < 0  then delete;
run;


data met_0gy_vs_1gy_2a;
   set met_0gy_vs_1gy_2;
   if C0gy_R1 < 0 or C0gy_R2 < 0  or  C1gy_R1 < 0 or C1gy_R2 < 0  then delete;
run;



proc export data=met_0gy_vs_01gy_2a
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_01gy_GC_noneg.txt"
     dbms=tab replace;
run;


proc export data=met_0gy_vs_1gy_2a
     outfile="/TB14/TB14/sandbox/dtra_sandbox/at_rad_metilene_0gy_1gy_GC_noneg.txt"
     dbms=tab replace;
run;




data sites2keep_gc01;
  set met_0gy_vs_01gy_2a;
  keep chr pos;
  rename pos=base;
run;

data sites2keep_gc1;
  set met_0gy_vs_1gy_2a;
  keep chr pos;
  rename pos=base;
run;

proc sort data=sites2keep_gc01 nodup;
  by chr base;
proc sort data=sites2keep_gc1 nodup;
  by chr base;

proc sort data=methylkit_01gy_0u_r1 nodup;
  by chr base;
proc sort data=methylkit_01gy_0u_r2 nodup;
  by chr base;
proc sort data=methylkit_01gy_100u_r1 nodup;
  by chr base;
proc sort data=methylkit_01gy_100u_r2 nodup;
  by chr base;

proc sort data=methylkit_0gy_0u_r1 nodup;
  by chr base;
proc sort data=methylkit_0gy_0u_r2 nodup;
  by chr base;
proc sort data=methylkit_0gy_100u_r1 nodup;
  by chr base;
proc sort data=methylkit_0gy_100u_r2 nodup;
  by chr base;

proc sort data=methylkit_1gy_0u_r1 nodup;
  by chr base;
proc sort data=methylkit_1gy_0u_r2 nodup;
  by chr base;
proc sort data=methylkit_1gy_100u_r1 nodup;
  by chr base;
proc sort data=methylkit_1gy_100u_r2 nodup;
  by chr base;
run;

%macro exportFiltered(input,sites,output);


data  &input._2;
  merge &input. (in=in1) &sites. (in=in2);
  by chr base;
  if in1 and in2;
run;

proc export data=&input._2 outfile="/TB14/TB14/sandbox/dtra_sandbox/&output..txt" dbms=tab replace; run;

%mend;


%exportFiltered(methylkit_0gy_0u_r1,sites2keep_gc01,at_rad_methylkit_0gy_0U_R1_v01);
%exportFiltered(methylkit_0gy_0u_r2,sites2keep_gc01,at_rad_methylkit_0gy_0U_R2_v01);
%exportFiltered(methylkit_0gy_100u_r1,sites2keep_gc01,at_rad_methylkit_0gy_100U_R1_v01);
%exportFiltered(methylkit_0gy_100u_r2,sites2keep_gc01,at_rad_methylkit_0gy_100U_R2_v01);
%exportFiltered(methylkit_01gy_0u_r1,sites2keep_gc01,at_rad_methylkit_01gy_0U_R1_v01);
%exportFiltered(methylkit_01gy_0u_r2,sites2keep_gc01,at_rad_methylkit_01gy_0U_R2_v01);
%exportFiltered(methylkit_01gy_100u_r1,sites2keep_gc01,at_rad_methylkit_01gy_100U_R1_v01);
%exportFiltered(methylkit_01gy_100u_r2,sites2keep_gc01,at_rad_methylkit_01gy_100U_R2_v01);

%exportFiltered(methylkit_0gy_0u_r1,sites2keep_gc1,at_rad_methylkit_0gy_0U_R1_v1);
%exportFiltered(methylkit_0gy_0u_r2,sites2keep_gc1,at_rad_methylkit_0gy_0U_R2_v1);
%exportFiltered(methylkit_0gy_100u_r1,sites2keep_gc1,at_rad_methylkit_0gy_100U_R1_v1);
%exportFiltered(methylkit_0gy_100u_r2,sites2keep_gc1,at_rad_methylkit_0gy_100U_R2_v1);
%exportFiltered(methylkit_1gy_0u_r1,sites2keep_gc1,at_rad_methylkit_1gy_0U_R1_v1);
%exportFiltered(methylkit_1gy_0u_r2,sites2keep_gc1,at_rad_methylkit_1gy_0U_R2_v1);
%exportFiltered(methylkit_1gy_100u_r1,sites2keep_gc1,at_rad_methylkit_1gy_100U_R1_v1);
%exportFiltered(methylkit_1gy_100u_r2,sites2keep_gc1,at_rad_methylkit_1gy_100U_R2_v1);


