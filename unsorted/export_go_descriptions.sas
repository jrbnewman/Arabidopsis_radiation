
libname tair "!PATCON/useful_arabidopsis_data/TAIR10/sas_data";

data go_terms;
  set tair.gene2go;
    length go_id $1200.;
      go_id=tranwrd(GO_number_biopro_cat, "GO", "go");
        keep go_id go_biological_process_cat;
	run;

	data go_terms2;
	  length go_id2 $12.;
	    length go_description $255.;
	      set go_terms;
	        do i=1 by 1 while(scan(go_id,i,"|") ^= "");
		      go_id2=scan(go_id,i,"|");
		            go_description=scan(go_biological_process_cat,i,"|");
			          output; end;
				    keep go_id2 go_description;
				      rename go_id2=go_id;
				      run;

				      proc sort data=go_terms2 nodup;
				         by go_id go_description;
					 run;

					 proc freq data=go_terms2 noprint;
					   tables go_id /out=check;
					   proc sort data=check;
					     by descending count;
					     run;


					     proc export data=go_terms2
					          outfile="!PATCON/DTRA/arabidopsis_go_terms.csv"
						       dbms=csv replace;
						       run;

