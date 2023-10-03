ods listing; ods html close;
libname rs '!PATCON/arabidopsis/sas_data';
libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';

/* update corresponding up/down indicators, given FDR-corrected P values

   If gene is DD, then we want to flag and indicate if Mock or treat only */

data dtct;
   set rs.arab_flag_gene_on_cpm_gt0;
   drop flag_gene_on_01gy_cpm0 flag_gene_on_1gy_cpm0 flag_gene_on_mock_cpm0 
        flag_gene_on_1hr_cpm0 flag_gene_on_3hr_cpm0 flag_gene_on_24hr_cpm0
        flag_gene_on_72hr_cpm0;
run;

data fdr;
  set rs.fdr_by_gene_log_cpm;
run;

proc sort data=dtct;
  by gene_id;
proc sort data=fdr;
  by gene_id;
run;

data fdr_w_sign;
  merge fdr (in=in1) dtct (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Redo FC flags */

%macro signFDR(fdrsuffix);

data fdr_w_sign_&fdrsuffix.;
  set fdr_w_sign;
  length sign_01gy_v_Mock_1h_&fdrsuffix. $5.;
  length sign_01gy_v_Mock_3h_&fdrsuffix. $5.;
  length sign_01gy_v_Mock_24h_&fdrsuffix. $5.;
  length sign_01gy_v_Mock_72h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_1h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_3h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_24h_&fdrsuffix. $5.;
  length sign_1gy_v_Mock_72h_&fdrsuffix. $5.;

  /* 0.1gy vs Mock */
  if flag_gene_on_01gy_1hr_cpm0=1 and flag_gene_on_Mock_1hr_cpm0=1 then do;
     if flag_01gy_v_Mock_1h_&fdrsuffix.=1 then sign_01gy_v_Mock_1h_&fdrsuffix.=sign_01gy_v_Mock_1h;
     else sign_01gy_v_Mock_1h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_1hr_cpm0=1 and flag_gene_on_Mock_1hr_cpm0 ne 1 then sign_01gy_v_Mock_1h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_1hr_cpm0 ne 1 and flag_gene_on_Mock_1hr_cpm0=1 then sign_01gy_v_Mock_1h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_1h_&fdrsuffix.="";

  if flag_gene_on_01gy_3hr_cpm0=1 and flag_gene_on_Mock_3hr_cpm0=1 then do;
     if flag_01gy_v_Mock_3h_&fdrsuffix.=1 then sign_01gy_v_Mock_3h_&fdrsuffix.=sign_01gy_v_Mock_3h;
     else sign_01gy_v_Mock_3h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_3hr_cpm0=1 and flag_gene_on_Mock_3hr_cpm0 ne 1 then sign_01gy_v_Mock_3h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_3hr_cpm0 ne 1 and gene_gene_on_Mock_3hr_cpm0=1 then sign_01gy_v_Mock_3h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_3h_&fdrsuffix.="";

  if flag_gene_on_01gy_24hr_cpm0=1 and flag_gene_on_Mock_24hr_cpm0=1 then do;
     if flag_01gy_v_Mock_24h_&fdrsuffix.=1 then sign_01gy_v_Mock_24h_&fdrsuffix.=sign_01gy_v_Mock_24h;
     else sign_01gy_v_Mock_24h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_24hr_cpm0=1 and flag_gene_on_Mock_24hr_cpm0 ne 1 then sign_01gy_v_Mock_24h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_24hr_cpm0 ne 1 and flag_gene_on_Mock_24hr_cpm0=1 then sign_01gy_v_Mock_24h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_24h_&fdrsuffix.="";

  if flag_gene_on_01gy_72hr_cpm0=1 and flag_gene_on_Mock_72hr_cpm0=1 then do;
     if flag_01gy_v_Mock_72h_&fdrsuffix.=1 then sign_01gy_v_Mock_72h_&fdrsuffix.=sign_01gy_v_Mock_72h;
     else sign_01gy_v_Mock_72h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_01gy_72hr_cpm0=1 and flag_gene_on_Mock_72hr_cpm0=0 then sign_01gy_v_Mock_72h_&fdrsuffix.="0.1gy";
  else if flag_gene_on_01gy_72hr_cpm0 ne 1 and flag_gene_on_Mock_72hr_cpm0=1 then sign_01gy_v_Mock_72h_&fdrsuffix.="Mock";
  else sign_01gy_v_Mock_72h_&fdrsuffix.="";

  /* 1gy vs Mock */
  if flag_gene_on_1gy_1hr_cpm0=1 and flag_gene_on_Mock_1hr_cpm0=1 then do;
     if flag_1gy_v_Mock_1h_&fdrsuffix.=1 then sign_1gy_v_Mock_1h_&fdrsuffix.=sign_1gy_v_Mock_1h;
     else sign_1gy_v_Mock_1h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_1hr_cpm0=1 and flag_gene_on_Mock_1hr_cpm0 ne 1 then sign_1gy_v_Mock_1h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_1hr_cpm0 ne 1 and flag_gene_on_Mock_1hr_cpm0=1 then sign_1gy_v_Mock_1h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_1h_&fdrsuffix.="";

  if flag_gene_on_1gy_3hr_cpm0=1 and flag_gene_on_Mock_3hr_cpm0=1 then do;
     if flag_1gy_v_Mock_3h_&fdrsuffix.=1 then sign_1gy_v_Mock_3h_&fdrsuffix.=sign_1gy_v_Mock_3h;
     else sign_1gy_v_Mock_3h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_3hr_cpm0=1 and flag_gene_on_Mock_3hr_cpm0 ne 1 then sign_1gy_v_Mock_3h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_3hr_cpm0 ne 1 and flag_gene_on_Mock_3hr_cpm0=1 then sign_1gy_v_Mock_3h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_3h_&fdrsuffix.="";

  if flag_gene_on_1gy_24hr_cpm0=1 and flag_gene_on_Mock_24hr_cpm0=1 then do;
     if flag_1gy_v_Mock_24h_&fdrsuffix.=1 then sign_1gy_v_Mock_24h_&fdrsuffix.=sign_1gy_v_Mock_24h;
     else sign_1gy_v_Mock_24h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_24hr_cpm0=1 and flag_gene_on_Mock_24hr_cpm0 ne 1 then sign_1gy_v_Mock_24h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_24hr_cpm0 ne 1 and flag_gene_on_Mock_24hr_cpm0=1 then sign_1gy_v_Mock_24h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_24h_&fdrsuffix.="";

  if flag_xscript_on_1gy_72hr_cpm0=1 and flag_xscript_on_Mock_72hr_cpm0=1 then do;
     if flag_1gy_v_Mock_72h_&fdrsuffix.=1 then sign_1gy_v_Mock_72h_&fdrsuffix.=sign_1gy_v_Mock_72h;
     else sign_1gy_v_Mock_72h_&fdrsuffix.="N";
     end;
  else if flag_gene_on_1gy_72hr_cpm0=1 and flag_gene_on_Mock_72hr_cpm0 ne 1 then sign_1gy_v_Mock_72h_&fdrsuffix.="1gy";
  else if flag_gene_on_1gy_72hr_cpm0 ne 1 and flag_gene_on_Mock_72hr_cpm0=1 then sign_1gy_v_Mock_72h_&fdrsuffix.="Mock";
  else sign_1gy_v_Mock_72h_&fdrsuffix.="";

  keep gene_id sign_: flag_gene_on_: ;

run;

%mend;
%signFDR(fdr05);
%signFDR(fdr10);
%signFDR(fdr20);

proc sort data=fdr_w_sign_fdr05;
   by gene_id;
proc sort data=fdr_w_sign_fdr10;
   by gene_id;
proc sort data=fdr_w_sign_fdr20;
   by gene_id;
run;

data fdr_w_sign2;
  merge fdr_w_sign_fdr05 fdr_w_sign_fdr10 fdr_w_sign_fdr20;
  by gene_id;
  
if flag_gene_on_01gy_1hr_cpm0=1 and flag_gene_on_mock_1hr_cpm0=1 then flag_gene_01gy_v_Mock_1h_dd=0;
else if flag_gene_on_01gy_1hr_cpm0=1 or flag_gene_on_mock_1hr_cpm0=1 then flag_gene_01gy_v_Mock_1h_dd=1;
else flag_gene_01gy_v_Mock_1h_dd=0;

if flag_gene_on_01gy_3hr_cpm0=1 and flag_gene_on_mock_3hr_cpm0=1 then flag_gene_01gy_v_Mock_3h_dd=0;
else if flag_gene_on_01gy_3hr_cpm0=1 or flag_gene_on_mock_3hr_cpm0=1 then flag_gene_01gy_v_Mock_3h_dd=1;
else flag_gene_01gy_v_Mock_3h_dd=0;

if flag_gene_on_01gy_24hr_cpm0=1 and flag_gene_on_mock_24hr_cpm0=1 then flag_gene_01gy_v_Mock_24h_dd=0;
else if flag_gene_on_01gy_24hr_cpm0=1 or flag_gene_on_mock_24hr_cpm0=1 then flag_gene_01gy_v_Mock_24h_dd=1;
else flag_gene_01gy_v_Mock_24h_dd=0;

if flag_gene_on_01gy_72hr_cpm0=1 and flag_gene_on_mock_72hr_cpm0=1 then flag_gene_01gy_v_Mock_72h_dd=0;
else if flag_gene_on_01gy_72hr_cpm0=1 or flag_gene_on_mock_72hr_cpm0=1 then flag_gene_01gy_v_Mock_72h_dd=1;
else flag_gene_01gy_v_Mock_72h_dd=0;


if flag_gene_on_1gy_1hr_cpm0=1 and flag_gene_on_mock_1hr_cpm0=1 then flag_gene_1gy_v_Mock_1h_dd=0;
else if flag_gene_on_1gy_1hr_cpm0=1 or flag_gene_on_mock_1hr_cpm0=1 then flag_gene_1gy_v_Mock_1h_dd=1;
else flag_gene_1gy_v_Mock_1h_dd=0;

if flag_gene_on_1gy_3hr_cpm0=1 and flag_gene_on_mock_3hr_cpm0=1 then flag_gene_1gy_v_Mock_3h_dd=0;
else if flag_gene_on_1gy_3hr_cpm0=1 or flag_gene_on_mock_3hr_cpm0=1 then flag_gene_1gy_v_Mock_3h_dd=1;
else flag_gene_1gy_v_Mock_3h_dd=0;

if flag_gene_on_1gy_24hr_apn0=1 and flag_gene_on_mock_24hr_apn0=1 then flag_gene_1gy_v_Mock_24h_dd=0;
else if flag_gene_on_1gy_24hr_apn0=1 or flag_gene_on_mock_24hr_apn0=1 then flag_gene_1gy_v_Mock_24h_dd=1;
else flag_gene_1gy_v_Mock_24h_dd=0;

if flag_gene_on_1gy_72hr_cpm0=1 and flag_gene_on_mock_72hr_cpm0=1 then flag_gene_1gy_v_Mock_72h_dd=0;
else if flag_gene_on_1gy_72hr_cpm0=1 or flag_gene_on_mock_72hr_cpm0=1 then flag_gene_1gy_v_Mock_72h_dd=1;
else flag_gene_1gy_v_Mock_72h_dd=0;

run;


data  rs.arab_sign_by_contrast_gene_fdr2;
  set fdr_w_sign2;
run;


proc freq data=fdr_w_sign2;
   tables  sign_01gy_v_Mock_1h_fdr05  sign_01gy_v_Mock_3h_fdr05 sign_01gy_v_Mock_24h_fdr05
           sign_01gy_v_Mock_72h_fdr05 sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_3h_fdr05 
           sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_72h_fdr05 sign_01gy_v_Mock_1h_fdr10
           sign_01gy_v_Mock_3h_fdr10 sign_01gy_v_Mock_24h_fdr10 sign_01gy_v_Mock_72h_fdr10
           sign_1gy_v_Mock_1h_fdr10 sign_1gy_v_Mock_3h_fdr10 sign_1gy_v_Mock_24h_fdr10
           sign_1gy_v_Mock_72h_fdr10 sign_01gy_v_Mock_1h_fdr20 sign_01gy_v_Mock_3h_fdr20
           sign_01gy_v_Mock_24h_fdr20 sign_01gy_v_Mock_72h_fdr20 sign_1gy_v_Mock_1h_fdr20
           sign_1gy_v_Mock_3h_fdr20 sign_1gy_v_Mock_24h_fdr20 sign_1gy_v_Mock_72h_fdr20 
	   flag_gene_01gy_v_Mock_1h_dd flag_gene_01gy_v_Mock_3h_dd 
           flag_gene_01gy_v_Mock_24h_dd flag_gene_01gy_v_Mock_72h_dd
	   flag_gene_1gy_v_Mock_1h_dd flag_gene_1gy_v_Mock_3h_dd 
           flag_gene_1gy_v_Mock_24h_dd flag_gene_1gy_v_Mock_72h_dd

	   flag_gene_01gy_v_Mock_1h_dd*sign_01gy_v_Mock_1h_fdr05 flag_gene_01gy_v_Mock_3h_dd*sign_01gy_v_Mock_3h_fdr05
           flag_gene_01gy_v_Mock_24h_dd*sign_01gy_v_Mock_24h_fdr05 flag_gene_01gy_v_Mock_72h_dd*sign_01gy_v_Mock_72h_fdr05
	   flag_gene_1gy_v_Mock_1h_dd*sign_1gy_v_Mock_1h_fdr05 flag_gene_1gy_v_Mock_3h_dd*sign_1gy_v_Mock_3h_fdr05
           flag_gene_1gy_v_Mock_24h_dd*sign_1gy_v_Mock_24h_fdr05 flag_gene_1gy_v_Mock_72h_dd*sign_1gy_v_Mock_72h_fdr05
           ;
 run;

/*

*/

















