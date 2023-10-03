ods listing;
ods html close;

libname cold '!PATCON/DTRA/arabidopsis_wgbs_cold/sas_data';
libname coldloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';

proc datasets lib=work kill noprint;
run;
quit;

/* Import HOMER annotations and compare base tair10 to ensemble GTF and whether or not strand makes a difference -- use superwindows for this */


%macro importHOMER(dataIN,dataOUT);

    data &dataOUT.   ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!PATCON/DTRA/arabidopsis_wgbs_cold/analysis_output/HOMER_annotation/&dataIN." delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat regionID $50.;
       informat Chr $3. ;
       informat Start best32. ;
       informat End best32. ;
       informat Strand $1. ;
       informat Peak_Score best32. ;
       informat Focus_Ratio_Region_Size $23. ;
       informat Annotation $100. ;
       informat Detailed_Annotation $100. ;
       informat Distance_to_TSS best32. ;
       informat Nearest_PromoterID $18. ;
       informat Entrez_ID best32. ;
       informat Nearest_Unigene $15. ;
       informat Nearest_Refseq $14. ;
       informat Nearest_Ensembl $15. ;
       informat Gene_Name $15. ;
       informat Gene_Alias $500. ;
       informat Gene_Description $200. ;
       informat Gene_Type $14. ;
        format regionID $50. ;
        format Chr $3. ;
        format Start best12. ;
        format End best12. ;
        format Strand $1. ;
        format Peak_Score best12. ;
        format Focus_Ratio_Region_Size $23. ;
        format Annotation $100. ;
        format Detailed_Annotation $100. ;
        format Distance_to_TSS best12. ;
        format Nearest_PromoterID $18. ;
        format Entrez_ID best32. ;
        format Nearest_Unigene $15. ;
        format Nearest_Refseq $14. ;
        format Nearest_Ensembl $15. ;
        format Gene_Name $15. ;
        format Gene_Alias $500. ;
        format Gene_Description $200. ;
        format Gene_Type $14. ;
     input
                 regionID $
                 Chr $
                 Start
                 End
                 Strand $
                 Peak_Score
                 Focus_Ratio_Region_Size $
                 Annotation $
                 Detailed_Annotation $
                 Distance_to_TSS
                 Nearest_PromoterID $
                 Entrez_ID
                 Nearest_Unigene $
                 Nearest_Refseq $
                 Nearest_Ensembl $
                 Gene_Name $
                 Gene_Alias $
                 Gene_Description $
                 Gene_Type $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


%mend;

%importHOMER(all_windows_for_HOMER_annotated_02jrbn.txt,windows_annot);
%importHOMER(CG_sites_for_HOMER_annotated2.txt,CG_site_annot);
%importHOMER(CHG_sites_for_HOMER_annotated2.txt,CHG_site_annot);
%importHOMER(GC_sites_for_HOMER_annotated2.txt,GC_site_annot);
%importHOMER(CHH_sites_for_HOMER_annotated2_6.txt,CHH_site_annot_6);
%importHOMER(CHH_sites_for_HOMER_annotated2_7.txt,CHH_site_annot_7);
%importHOMER(CHH_sites_for_HOMER_annotated2_4.txt,CHH_site_annot_4);
%importHOMER(CHH_sites_for_HOMER_annotated2_2.txt,CHH_site_annot_2);
%importHOMER(CHH_sites_for_HOMER_annotated2_3.txt,CHH_site_annot_3);
%importHOMER(CHH_sites_for_HOMER_annotated2_1.txt,CHH_site_annot_1);
%importHOMER(CHH_sites_for_HOMER_annotated2_5.txt,CHH_site_annot_5);


/* Make permenant */

data cold.homer_annotations_windows_new;
   set windows_annot;
   length geneID $20.;
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"exon") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(upcase(scan(Nearest_PromoterID,1,".")));
   rename regionID=windowID;
run;

data cold.homer_annotations_sites;
   set CG_site_annot CHG_site_annot GC_site_annot
       CHH_site_annot_1 CHH_site_annot_2 CHH_site_annot_3 CHH_site_annot_4
       CHH_site_annot_5 CHH_site_annot_6 CHH_site_annot_7;
   length geneID $20.;
   if count(annotation,"promoter") >= 1 then flag_promoter=1; else flag_promoter=0;
   if count(annotation,"exon") >= 1 or count(annotation,"exon") >= 1 then flag_genebody=1; else flag_genebody=0;
   geneID=compress(upcase(scan(Nearest_PromoterID,1,".")));
   rename regionID=siteID;
run;

