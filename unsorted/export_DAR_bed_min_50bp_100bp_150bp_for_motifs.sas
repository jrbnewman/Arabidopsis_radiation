/* Counts for arabidopsis paper */

libname arabRNA '!HOME/concannon/DTRA/arabidopsis/sas_data';
libname arabMAP '/TB14/TB14/sandbox/dtra_sandbox/arab_wgbs';
libname wgbsA '!HOME/concannon/DTRA/arabidopsis_wgbs/sas_data';

ods html close;
ods listing;

proc datasets kill lib=work  noprint;
run;
quit;


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER_annotation.txt"
	out=dar_annot dbms=tab replace;
	guessingrows=max;
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


data results_by_dar;
	set wgbsA.results_by_dar_5sites_50bp;
run;

proc freq data=results_by_dar;
	tables comparison*flag_fdr05;
run;


proc sort data=dar_annot2;
	by comparison site_type chr  region_num;
proc sort data=results_by_dar;
	by comparison site_type chr  region_num;
run;


data dar_w_annot;
	merge dar_annot2 (in=in1) results_by_dar (in=in2);
	by comparison site_type chr region_num;
	if in1 and in2;
run;



data dar_w_annot2;
	set dar_w_annot;
	if mean_methyl_diff < 0 then flag_direction=-1;
	else if  mean_methyl_diff > 0 then flag_direction=1;
	else mean_methyl_diff = 0;
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

proc sort data=gc_01_up nodup; by _all_;
proc sort data=gc_01_dn nodup; by _all_;
proc sort data=gc_1_up nodup; by _all_;
proc sort data=gc_1_dn nodup; by _all_;
run;

proc freq data=gc_01_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_01_dn;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_dn;  tables flag_DAC10_ge2*flag_fdr05; run;


%macro exportBED(inputData);

	data format_for_bedfile;
		retain chr region_start2 region_stop2 region_id;
		set &inputData.;
		length region_id $100.;
		region_start2 = region_start - 6;
		region_stop2 = region_stop + 6;
		region_id=compress(catx("_",chr,region_start,region_stop));
		keep chr region_start2 region_stop2 region_id;
	run;

	proc sort data=format_for_bedfile nodup;
		by chr region_start2 region_stop2 region_id;
	run;


	proc export data=format_for_bedfile
		outfile="!HOME/concannon/DTRA/at_rad_motif_analysis/input_data/min5_50bp_&inputData..bed"
		dbms=tab replace;
		putnames=no;
	run;

%mend;


%exportBED(gc_01_up);
%exportBED(gc_1_up);
%exportBED(gc_01_dn);
%exportBED(gc_1_dn);







proc datasets kill lib=work  noprint;
run;
quit;


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER_annotation.txt"
	out=dar_annot dbms=tab replace;
	guessingrows=max;
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


data results_by_dar;
	set wgbsA.results_by_dar_5sites_100bp;
run;

proc freq data=results_by_dar;
	tables comparison*flag_fdr05;
run;


proc sort data=dar_annot2;
	by comparison site_type chr  region_num;
proc sort data=results_by_dar;
	by comparison site_type chr  region_num;
run;


data dar_w_annot;
	merge dar_annot2 (in=in1) results_by_dar (in=in2);
	by comparison site_type chr region_num;
	if in1 and in2;
run;



data dar_w_annot2;
	set dar_w_annot;
	if mean_methyl_diff < 0 then flag_direction=-1;
	else if  mean_methyl_diff > 0 then flag_direction=1;
	else mean_methyl_diff = 0;
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

proc sort data=gc_01_up nodup; by _all_;
proc sort data=gc_01_dn nodup; by _all_;
proc sort data=gc_1_up nodup; by _all_;
proc sort data=gc_1_dn nodup; by _all_;
run;

proc freq data=gc_01_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_01_dn;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_dn;  tables flag_DAC10_ge2*flag_fdr05; run;


%macro exportBED(inputData);

	data format_for_bedfile;
		retain chr region_start2 region_stop2 region_id;
		set &inputData.;
		length region_id $100.;
		region_start2 = region_start - 6;
		region_stop2 = region_stop + 6;
		region_id=compress(catx("_",chr,region_start,region_stop));
		keep chr region_start2 region_stop2 region_id;
	run;

	proc sort data=format_for_bedfile nodup;
		by chr region_start2 region_stop2 region_id;
	run;


	proc export data=format_for_bedfile
		outfile="!HOME/concannon/DTRA/at_rad_motif_analysis/input_data/min5_100bp_&inputData..bed"
		dbms=tab replace;
		putnames=no;
	run;

%mend;


%exportBED(gc_01_up);
%exportBED(gc_1_up);
%exportBED(gc_01_dn);
%exportBED(gc_1_dn);



proc datasets kill lib=work  noprint;
run;
quit;


proc import datafile="/TB14/TB14/sandbox/dtra_sandbox/DARs_min_5_sites_for_HOMER_annotation.txt"
	out=dar_annot dbms=tab replace;
	guessingrows=max;
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


data results_by_dar;
	set wgbsA.results_by_dar_5sites_150bp;
run;

proc freq data=results_by_dar;
	tables comparison*flag_fdr05;
run;


proc sort data=dar_annot2;
	by comparison site_type chr  region_num;
proc sort data=results_by_dar;
	by comparison site_type chr  region_num;
run;


data dar_w_annot;
	merge dar_annot2 (in=in1) results_by_dar (in=in2);
	by comparison site_type chr region_num;
	if in1 and in2;
run;



data dar_w_annot2;
	set dar_w_annot;
	if mean_methyl_diff < 0 then flag_direction=-1;
	else if  mean_methyl_diff > 0 then flag_direction=1;
	else mean_methyl_diff = 0;
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

proc sort data=gc_01_up nodup; by _all_;
proc sort data=gc_01_dn nodup; by _all_;
proc sort data=gc_1_up nodup; by _all_;
proc sort data=gc_1_dn nodup; by _all_;
run;

proc freq data=gc_01_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_01_dn;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_up;  tables flag_DAC10_ge2*flag_fdr05; run;
proc freq data=gc_1_dn;  tables flag_DAC10_ge2*flag_fdr05; run;


%macro exportBED(inputData);

	data format_for_bedfile;
		retain chr region_start2 region_stop2 region_id;
		set &inputData.;
		length region_id $100.;
		region_start2 = region_start - 6;
		region_stop2 = region_stop + 6;
		region_id=compress(catx("_",chr,region_start,region_stop));
		keep chr region_start2 region_stop2 region_id;
	run;

	proc sort data=format_for_bedfile nodup;
		by chr region_start2 region_stop2 region_id;
	run;


	proc export data=format_for_bedfile
		outfile="!HOME/concannon/DTRA/at_rad_motif_analysis/input_data/min5_150bp_&inputData..bed"
		dbms=tab replace;
		putnames=no;
	run;

%mend;


%exportBED(gc_01_up);
%exportBED(gc_1_up);
%exportBED(gc_01_dn);
%exportBED(gc_1_dn);




