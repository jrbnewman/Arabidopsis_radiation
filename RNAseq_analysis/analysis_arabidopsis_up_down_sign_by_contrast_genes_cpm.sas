ods listing; ods html close;

libname tair '!PATCON/useful_arabidopsis_data/TAIR10/sas_data';
libname rs '!PATCON/arabidopsis/sas_data';

/* For each comparison, flag gene up/down

Going to use the estimates from contrasts to decide what is up/down.

If P<0.05 and estimate is negative, then fusion is "down"
If P<0.05 and estimate is positive, then fusion is "up"
else fusion is not diff

*/
data estimates;
  set rs.arab_gene_cntrs_estim_lcpm;
  keep gene_id label estimate;
run;

data contrasts;
  set rs.arab_gene_cntrs_constr_lcpm;
  keep gene_id label ProbF;
run;

proc sort data=estimates;
   by gene_id label;
proc sort data=contrasts;
   by gene_id label;
run;

data contrast_estim;
  merge contrasts (in=in1) estimates (in=in2);
  by gene_id label;
  if in1 and in2;
run;

/* Flag if up/down */

data up_down;
   set contrast_estim;
   length sign $1.;
   if ProbF = . then sign="";
   else if ProbF < 0.05 then do;
      if estimate < 0 then sign="D";
      else sign = "U";
      end;
   else sign="N";
run;

/* Add in a variable name to transpose on */

data sign_varname;
   set up_down;
   length sign_label $32.;

  if label="0.1gy-Mock: 1h" then sign_label="sign_01gy_v_Mock_1h";
  if label= "0.1gy-Mock: 3h" then sign_label="sign_01gy_v_Mock_3h";
  if label= "0.1gy-Mock: 24h" then sign_label="sign_01gy_v_Mock_24h";
  if label= "0.1gy-Mock: 72h" then sign_label="sign_01gy_v_Mock_72h";

  if label= "1gy-Mock: 1h" then sign_label="sign_1gy_v_Mock_1h";
  if label= "1gy-Mock: 3h" then sign_label="sign_1gy_v_Mock_3h";
  if label= "1gy-Mock: 24h" then sign_label="sign_1gy_v_Mock_24h";
  if label= "1gy-Mock: 72h" then sign_label="sign_1gy_v_Mock_72h";

  if label="1gy-0.1gy: 1h" then sign_label="sign_1gy_v_01gy_1h";
  if label="1gy-0.1gy: 3h" then sign_label="sign_1gy_v_01gy_3h";
  if label="1gy-0.1gy: 24h" then sign_label="sign_1gy_v_01gy_24h";
  if label="1gy-0.1gy: 72h" then sign_label="sign_1gy_v_01gy_72h";

  if label="Mock: 3h-1h" then sign_label="sign_mock_3h_v_1h";
  if label="Mock: 24h-1h" then sign_label="sign_mock_24h_v_1h";
  if label="Mock: 72h-1h" then sign_label="sign_mock_72h_v_1h";
  if label="Mock: 24h-3h" then sign_label="sign_mock_24h_v_3h";
  if label="Mock: 72h-3h" then sign_label="sign_mock_72h_v_3h";
  if label="Mock: 72h-24h" then sign_label="sign_mock_72h_v_24h";

  if label="0.1gy: 3h-1h" then sign_label="sign_01gy_3h_v_1h";
  if label="0.1gy: 24h-1h" then sign_label="sign_01gy_24h_v_1h";
  if label="0.1gy: 72h-1h" then sign_label="sign_01gy_72h_v_1h";
  if label="0.1gy: 24h-3h" then sign_label="sign_01gy_24h_v_3h";
  if label="0.1gy: 72h-3h" then sign_label="sign_01gy_72h_v_3h";
  if label="0.1gy: 72h-24h" then sign_label="sign_01gy_72h_v_24h";

  if label="1gy: 3h-1h" then sign_label="sign_1gy_3h_v_1h";
  if label="1gy: 24h-1h" then sign_label="sign_1gy_24h_v_1h";
  if label="1gy: 72h-1h" then sign_label="sign_1gy_72h_v_1h";
  if label="1gy: 24h-3h" then sign_label="sign_1gy_24h_v_3h";
  if label="1gy: 72h-3h" then sign_label="sign_1gy_72h_v_3h";
  if label="1gy: 72h-24h" then sign_label="sign_1gy_72h_v_24h";

  if label="0.1gy 3h-1h = Mock 3h-1h" then sign_label="sign_01gy_v_Mock_3h_sub_1h";
  if label="0.1gy 24h-1h = Mock 24h-1h" then sign_label="sign_01gy_v_Mock_24h_sub_1h";
  if label="0.1gy 72h-1h = Mock 72h-1h" then sign_label="sign_01gy_v_Mock_72h_sub_1h";
  if label="1gy 3h-1h = Mock 3h-1h" then sign_label="sign_1gy_v_Mock_3h_sub_1h";
  if label="1gy 24h-1h = Mock 24h-1h" then sign_label="sign_1gy_v_Mock_24h_sub_1h";
  if label="1gy 72h-1h = Mock 72h-1h" then sign_label="sign_1gy_v_Mock_72h_sub_1h";

keep gene_id sign sign_label;
run;

proc sort data= sign_varname;
   by gene_id sign_label;
proc transpose data=sign_varname out=sign_sbys;
   by gene_id;
   id sign_label;
   var sign;
run;


data rs.arab_sign_by_contrast_gene_cpm;
   set sign_sbys;
   drop _NAME_;
run;


proc freq data=sign_sbys noprint;
   tables sign_01gy_v_Mock_3h_sub_1h*sign_01gy_v_Mock_24h_sub_1h*sign_01gy_v_Mock_72h_sub_1h
            / out=sign_check;
run;

proc print data=sign_check;
run;


proc freq data=sign_sbys noprint;
   tables sign_1gy_v_Mock_3h_sub_1h*sign_1gy_v_Mock_24h_sub_1h*sign_1gy_v_Mock_72h_sub_1h
            / out=sign_check;
run;

proc print data=sign_check;
run;


