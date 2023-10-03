/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;

/* compare TE counting methods */

data bwa_multi;
set arabrna.counts_by_te_bwa_multi;
cpm = region_depth * 10000000 / mapped_reads;
length TE_id $50.;
length TE_locus $50.;
TE_id=compress(scan(element_id,1,"/"));
TE_locus=compress(scan(element_id,2,"/"));
keep sample_id TE_id  TE_locus sample_id apn cpm;
rename apn=bwa_apn_multi cpm=bwa_cpm_multi;
run;

data bwa_uniq;
set arabrna.counts_by_te_bwa_uniq;
cpm = region_depth * 10000000 / mapped_reads;
length TE_id $50.;
length TE_locus $50.;
TE_id=compress(scan(element_id,1,"/"));
TE_locus=compress(scan(element_id,2,"/"));
keep sample_id TE_id  TE_locus sample_id apn cpm;rename apn=bwa_apn_uniq cpm=bwa_cpm_uniq;
run;

data local_multi;
   set arabrna.cpm_by_te_local_multi;
   where flag_TE=1;
   length TE_id $50.;
   length TE_locus $50.;
   TE_id=compress(scan(gene_id,2,":"));
   TE_locus=compress(scan(gene_id,1,":"));
   keep sample_id TE_id TE_locus cpm;
   rename cpm=local_cpm_multi;
run;


data local_uniq;
   set arabrna.cpm_by_te_local_uniq;
   where flag_TE=1;
   length TE_id $50.;
   length TE_locus $50.;
   TE_id=compress(scan(gene_id,2,":"));
   TE_locus=compress(scan(gene_id,1,":"));
   keep sample_id TE_id TE_locus cpm;
   rename cpm=local_cpm_uniq;
run;
data xs_multi;
   set arabrna.cpm_by_te_xscript_multi;
   where flag_TE=1;
   length TE_id $50.;
   length TE_locus $50.;
   TE_id=compress(scan(gene_id,1,"/"));
   TE_locus=compress(scan(scan(gene_id,1,":"),2,"/"));
   keep sample_id TE_id TE_locus cpm;
   rename cpm=xs_cpm_multi;
run;


data xs_uniq;
   set arabrna.cpm_by_te_xscript_uniq;
   where flag_TE=1;
   length TE_id $50.;
   length TE_locus $50.;
   TE_id=compress(scan(gene_id,1,"/"));
   TE_locus=compress(scan(scan(gene_id,1,":"),2,"/"));
   keep sample_id TE_id TE_locus cpm;
   rename cpm=xs_cpm_uniq;
run;



proc sort data=bwa_multi; by sample_id TE_id TE_locus;
proc sort data=bwa_uniq; by sample_id TE_id TE_locus;
proc sort data=local_multi; by sample_id TE_id TE_locus;
proc sort data=local_uniq; by sample_id TE_id TE_locus;
proc sort data=xs_multi; by sample_id TE_id TE_locus;
proc sort data=xs_uniq; by sample_id TE_id TE_locus;
run;

data te_locus_compare;
  merge bwa_multi bwa_uniq local_multi local_uniq xs_multi xs_uniq;
  by sample_id TE_id TE_locus;
run;

proc corr data=te_locus_compare pearson;
   by sample_id;
   var bwa_apn_multi bwa_cpm_multi local_cpm_multi xs_cpm_multi
       bwa_apn_uniq bwa_cpm_uniq local_cpm_uniq xs_cpm_uniq;
run;

/*  CORRELATIONS: (one sample)


                                                             Simple Statistics

                     Variable                  N          Mean       Std Dev           Sum       Minimum       Maximum

                     bwa_apn_multi         31189       0.77125      19.37575         24054             0          2943
                     bwa_cpm_multi         31189      29.70763     512.85256        926551             0         60676
                     local_cpm_multi       31189       0.06007       1.17401          1873             0      72.68250
                     xs_cpm_multi          31189       0.06007       1.17400          1873             0      72.68139
                     bwa_apn_uniq          31189       3.87353      96.42249        120812             0         14614
                     bwa_cpm_uniq          31189     149.78863          2574       4671758             0        304620
                     local_cpm_uniq        31189       0.04933       1.09144          1539             0      70.11011
                     xs_cpm_uniq           31189       0.04933       1.09141          1538             0      70.10813


                                                Pearson Correlation Coefficients, N = 31189
                                                         Prob > |r| under H0: Rho=0

                           bwa_apn_      bwa_cpm_         local_       xs_cpm_      bwa_apn_      bwa_cpm_        local_       xs_cpm_
                              multi         multi      cpm_multi         multi          uniq          uniq      cpm_uniq          uniq

      bwa_apn_multi         1.00000       0.84128        0.03969       0.03969       0.99967       0.84238       0.03307       0.03307
                                           <.0001         <.0001        <.0001        <.0001        <.0001        <.0001        <.0001

      bwa_cpm_multi         0.84128       1.00000        0.34504       0.34504       0.84030       0.99933       0.32455       0.32455
                             <.0001                       <.0001        <.0001        <.0001        <.0001        <.0001        <.0001

      local_cpm_multi       0.03969       0.34504        1.00000       1.00000       0.03991       0.34689       0.93666       0.93666
                             <.0001        <.0001                       <.0001        <.0001        <.0001        <.0001        <.0001

      xs_cpm_multi          0.03969       0.34504        1.00000       1.00000       0.03991       0.34689       0.93666       0.93666
                             <.0001        <.0001         <.0001                      <.0001        <.0001        <.0001        <.0001

      bwa_apn_uniq          0.99967       0.84030        0.03991       0.03991       1.00000       0.84198       0.03334       0.03334
                             <.0001        <.0001         <.0001        <.0001                      <.0001        <.0001        <.0001

       bwa_cpm_uniq          0.84238       0.99933        0.34689       0.34689       0.84198       1.00000       0.32658       0.32658
                              <.0001        <.0001         <.0001        <.0001        <.0001                      <.0001        <.0001

       local_cpm_uniq        0.03307       0.32455        0.93666       0.93666       0.03334       0.32658       1.00000       1.00000
                              <.0001        <.0001         <.0001        <.0001        <.0001        <.0001                      <.0001

       xs_cpm_uniq           0.03307       0.32455        0.93666       0.93666       0.03334       0.32658       1.00000       1.00000
                              <.0001        <.0001         <.0001        <.0001        <.0001        <.0001        <.0001
*/

/*sum and merge with gene-level counts */

proc sort data=te_locus_compare;
   by sample_id TE_id;
proc means data=te_locus_compare noprint;
  by sample_id TE_id;
   var bwa_apn_multi bwa_cpm_multi local_cpm_multi xs_cpm_multi
       bwa_apn_uniq bwa_cpm_uniq local_cpm_uniq xs_cpm_uniq;
  output out=te_locus_compare_summed sum=;
run;


data te_multi;
   set arabrna.cpm_by_te_multi;
   where flag_TE=1;
   length TE_id $50.;
   TE_id=compress(scan(gene_id,1,":"));
   keep sample_id TE_id cpm;
   rename cpm=te_cpm_multi;
run;


data te_uniq;
   set arabrna.cpm_by_te_uniq;
   where flag_TE=1;
   length TE_id $50.;
   TE_id=compress(scan(gene_id,1,":"));
   keep sample_id TE_id cpm;
   rename cpm=te_cpm_uniq;
run;

proc sort data=te_locus_compare_summed;
  by sample_id te_id;
proc sort data=te_multi;
  by sample_id te_id;
proc sort data=te_uniq;
  by sample_id te_id;
run;

data te_compare;
  merge te_locus_compare_summed te_multi te_uniq;
  by samplE_id te_id;
run;



proc corr data=te_compare pearson;
   by sample_id;
   var bwa_apn_multi bwa_cpm_multi local_cpm_multi xs_cpm_multi
       bwa_apn_uniq bwa_cpm_uniq local_cpm_uniq xs_cpm_uniq te_cpm_multi te_cpm_uniq;
run;


/*

                                          Simple Statistics

  Variable                  N          Mean       Std Dev           Sum       Minimum       Maximum

  bwa_apn_multi           320      75.17017     255.84732         24054             0          2945
  bwa_cpm_multi           320          2895          8837        926551             0         97212
  local_cpm_multi         320       5.85445      23.81535          1873             0     374.76641
  xs_cpm_multi            320       5.85437      23.81498          1873             0     374.76070
  bwa_apn_uniq            320     377.53621          1277        120812             0         14627
  bwa_cpm_uniq            320         14599         44536       4671758             0        490021
  local_cpm_uniq          320       4.80782      23.12263          1539             0     371.99014
  xs_cpm_uniq             320       4.80769      23.12197          1538             0     371.97962
  te_cpm_multi            320       5.92610      23.84133          1896             0     374.92630
  te_cpm_uniq             320       4.81991      23.12563          1542             0     372.00741


                                                            Pearson Correlation Coefficients, N = 320
                                                                    Prob > |r| under H0: Rho=0

                        bwa_apn_      bwa_cpm_         local_       xs_cpm_      bwa_apn_      bwa_cpm_        local_       xs_cpm_       te_cpm_       te_cpm_
                           multi         multi      cpm_multi         multi          uniq          uniq      cpm_uniq          uniq         multi          uniq

   bwa_apn_multi         1.00000       0.83012        0.07197       0.07197       0.99983       0.82980       0.05479       0.05479       0.07277       0.05495
                                        <.0001         0.1991        0.1991        <.0001        <.0001        0.3285        0.3285        0.1942        0.3272

   bwa_cpm_multi         0.83012       1.00000        0.48337       0.48337       0.83094       0.99985       0.46085       0.46085       0.48463       0.46111
                          <.0001                       <.0001        <.0001        <.0001        <.0001        <.0001        <.0001        <.0001        <.0001

   local_cpm_multi       0.07197       0.48337        1.00000       1.00000       0.07156       0.48471       0.98265       0.98265       0.99998       0.98273
                          0.1991        <.0001                       <.0001        0.2017        <.0001        <.0001        <.0001        <.0001        <.0001

   xs_cpm_multi          0.07197       0.48337        1.00000       1.00000       0.07156       0.48471       0.98265       0.98265       0.99998       0.98273
                          0.1991        <.0001         <.0001                      0.2017        <.0001        <.0001        <.0001        <.0001        <.0001

   bwa_apn_uniq          0.99983       0.83094        0.07156       0.07156       1.00000       0.83080       0.05438       0.05438       0.07237       0.05453
                          <.0001        <.0001         0.2017        0.2017                      <.0001        0.3322        0.3322        0.1966        0.3308

   bwa_cpm_uniq          0.82980       0.99985        0.48471       0.48471       0.83080       1.00000       0.46215       0.46215       0.48597       0.46241
                          <.0001        <.0001         <.0001        <.0001        <.0001                      <.0001        <.0001        <.0001        <.0001

   local_cpm_uniq        0.05479       0.46085        0.98265       0.98265       0.05438       0.46215       1.00000       1.00000       0.98228       1.00000
                          0.3285        <.0001         <.0001        <.0001        0.3322        <.0001                      <.0001        <.0001        <.0001

   xs_cpm_uniq           0.05479       0.46085        0.98265       0.98265       0.05438       0.46215       1.00000       1.00000       0.98228       1.00000
                          0.3285        <.0001         <.0001        <.0001        0.3322        <.0001        <.0001                      <.0001        <.0001

   te_cpm_multi          0.07277       0.48463        0.99998       0.99998       0.07237       0.48597       0.98228       0.98228       1.00000       0.98236
                          0.1942        <.0001         <.0001        <.0001        0.1966        <.0001        <.0001        <.0001                      <.0001

   te_cpm_uniq           0.05495       0.46111        0.98273       0.98273       0.05453       0.46241       1.00000       1.00000       0.98236       1.00000
                          0.3272        <.0001         <.0001        <.0001        0.3308        <.0001        <.0001        <.0001        <.0001



Use TElocal counts!
*/

