ods listing; ods html close;

libname rs "!PATCON/arabidopsis/sas_data";

 /* Import coverage counts for fusions and sum tech reps */

    data WORK.FUSION_COUNTS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!PATCON/arabidopsis/alignment_output/coverage_counts_fusions_arab.csv' delimiter = ',' MISSOVER DSD
lrecl=32767 firstobs=2 ;
       informat sample_id $78. ;
       informat fusion_id $10. ;
       informat mapped_reads best32. ;
       informat read_length best32. ;
       informat region_length best32. ;
       informat region_depth best32. ;
       informat reads_in_region best32. ;
       informat apn best32. ;
       informat rpkm best32. ;
       informat mean best32. ;
       informat std best32. ;
       informat cv best32. ;
       format sample_id $78. ;
       format fusion_id $10. ;
       format mapped_reads best12. ;
       format read_length best12. ;
       format region_length best12. ;
       format region_depth best12. ;
       format reads_in_region best12. ;
       format apn best12. ;
       format rpkm best12. ;
       format mean best12. ;
       format std best12. ;
       format cv best12. ;
     input
                 sample_id $
                 fusion_id $
                 mapped_reads
                 read_length
                 region_length
                 region_depth
                 reads_in_region
                 apn
                 rpkm
                 mean
                 std
                 cv
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

/* Merge in design file */

data design;
   set rs.arab_design_file;
   keep sample_id sample_number treatment time replicate;
run;

proc sort data=design nodup;
  by sample_id;
proc sort data=fusion_counts;
  by sample_id;
run;

data counts_w_key;
  merge design (in=in1) fusion_counts (in=in2);
  by sample_id;
  if in1 and in2;
run;

proc sort data=counts_w_key;
   by sample_number treatment time replicate fusion_id;
proc means data=counts_w_key noprint;
   by sample_number treatment time replicate  fusion_id;
   var apn;
   output out=summed_counts_w_key sum=;
run;

data rs.arab_fusion_counts_w_key;
  set summed_counts_w_key;
  drop _TYPE_ _FREQ_;
run;


