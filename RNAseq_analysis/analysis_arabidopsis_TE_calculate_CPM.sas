ods listing; ods html close;
libname arabRNA "!HOME/concannon/DTRA/arabidopsis/sas_data";

/* Normalize gene read counts using "copies per million" CPM: reads aligning to gene / 1,000,000
   [reads aligning to gene] * 10^6 / [total number of aligned reads in one sample] */

   %macro calcCPM(dataIn,dataOut);


	   data counts_by_gene;
		     set arabRNA.&dataIn.;
	     run;

	     /* sum expected_count by sample */

	     proc sort data=counts_by_gene;
		       by sample_id;
	       proc means data=counts_by_gene noprint;
		         by sample_id;
			   var count;
			     output out=total_counts_per_sample sum=total_counts;
		     run;

		     proc sort data=total_counts_per_sample;
			       by sample_id;
		       proc sort data=counts_by_gene;
			         by sample_id;
			 run;

			 data counts_by_gene2;
				   merge counts_by_gene (in=in1) total_counts_per_sample (in=in2);
				     by sample_id;
				       if in1 and in2;
			       run;

			       /* Calculate CPM from expected_counts */

			       data calc_cpm;
				         set counts_by_gene2;
					   cpm = (count * 1000000)/total_counts;
				   run;


				   data arabRNA.&dataOut.;
					     set calc_cpm;
					       drop _TYPE_ _FREQ_;
				       run;

   %mend;


   %calcCPM(counts_by_te_uniq,cpm_by_te_uniq);
   %calcCPM(counts_by_te_multi,cpm_by_te_multi);
   %calcCPM(counts_by_te_xscript_uniq,cpm_by_te_xscript_uniq);
   %calcCPM(counts_by_te_xscript_multi,cpm_by_te_xscript_multi);
   %calcCPM(counts_by_te_local_uniq,cpm_by_te_local_uniq);
   %calcCPM(counts_by_te_local_multi,cpm_by_te_local_multi);

