libname rs '!PATCON/arabidopsis/sas_data';
ods listing; ods html close;

/* ABA transporters:
			0.1gy	1gy
AT1G15520	UNNN	NNNN	*
AT1G71960	NNNN	NNNN
AT1G69850	UNUU	UUNN	***
AT5G52050	UNNN	UNNN	**

ABA homeostasis:
AT1G07240	UNNN	UNNN	**
AT1G42550	DNNN	DNNN	**

ABA catabolism:
AT3G19270	NNNN	DNNN	*
AT4G19230	NNNN	NNNN
AT2G29090	NNNN	NNNN
AT3G21780	NNNN	NNNN
AT5G45340	NNNN	NNNN

ABA glucosyltransferase:
AT4G34138	NNNN	NUNN
AT1G07240	UNNN	UNNN	**
AT1G05560	NNNN	NNNN
AT2G31750	NNNN	NUNN	*
AT3G21780	NNNN	NNNN
AT4G34131	NNNN	NNNN
AT1G05530	TNT-	TC-T	***

*/

/* check if genes are significant */

data sign_check;
  set rs.arab_sign_by_cntrst_gene_cpm_fdr;
  if gene_id in ('AT1G15520','AT1G71960','AT1G69850','AT5G52050',
                 'AT1G07240','AT1G42550','AT3G19270','AT4G19230',
                 'AT2G29090','AT3G21780','AT5G45340','AT4G34138',
                 'AT1G07240','AT1G05560','AT2G31750','AT2G23250',
                 'AT2G23260','AT3G21780','AT4G34131','AT1G05530');
  keep gene_id sign_01gy_v_Mock_1h_fdr05 sign_01gy_v_Mock_3h_fdr05
       sign_01gy_v_Mock_24h_fdr05 sign_01gy_v_Mock_72h_fdr05
       sign_1gy_v_Mock_1h_fdr05 sign_1gy_v_Mock_3h_fdr05
       sign_1gy_v_Mock_24h_fdr05 sign_1gy_v_Mock_72h_fdr05;
run;

proc print data=sign_check;
run;


/* What is the fold change? */

data fold_change;
  set rs.arab_gene_mean_cpm_by_trt_time;
  if gene_id in ('AT1G15520','AT1G69850','AT5G52050','AT1G07240',
                 'AT1G42550','AT3G19270','AT4G34138','AT1G07240',
                 'AT2G31750','AT1G05530');
  fc_01gy_1h=log(mean_cpm_01gy_1h/mean_cpm_mock_1h);
  fc_01gy_3h=log(mean_cpm_01gy_3h/mean_cpm_mock_3h);
  fc_01gy_24h=log(mean_cpm_01gy_24h/mean_cpm_mock_24h);
  fc_01gy_72h=log(mean_cpm_01gy_72h/mean_cpm_mock_72h);

  fc_1gy_1h=log(mean_cpm_1gy_1h/mean_cpm_mock_1h);
  fc_1gy_3h=log(mean_cpm_1gy_3h/mean_cpm_mock_3h);
  fc_1gy_24h=log(mean_cpm_1gy_24h/mean_cpm_mock_24h);
  fc_1gy_72h=log(mean_cpm_1gy_72h/mean_cpm_mock_72h);

  keep gene_id fc_01gy_1h fc_01gy_3h fc_01gy_24h fc_01gy_72h
       fc_1gy_1h fc_1gy_3h fc_1gy_24h fc_1gy_72h;
run;

/* nothin screaming interesting really, mostly pretty obvious... */

/*
ABA transporters
AT1G15520	UNNN	NNNN	*	upreg response to lead, lead resistance; ethylene, jasmonic acid, ozone, salicyclic acaid, terpenoid; highly expressed in mature flower and sensescent leaf
AT1G69850	UNUU	UUNN	***	nitrate uptake, response to nematodes, calcium binding; moderate expression in mature leaf and root
AT5G52050	UNNN	UNNN	**	reg response water deprivation, ABA response

ABA homeostasis:
AT1G07240	UNNN	UNNN	**	drought recovery, +reg seed germination, 
AT1G42550	DNNN	DNNN	**	response to osmotic stress, reg seed germination

ABA catabolism:		
AT3G19270	NNNN	DNNN	*	oxid-red process, sterol metabolism, 
				
ABA glucosyltransferase:	
AT4G34138	NNNN	NUNN	*	metabolic process
AT1G07240	UNNN	UNNN	**	drought recovery, +reg seed germination,
AT2G31750	NNNN	NUNN	*	auxin metabolic process,
AT1G05530	TNT-	TC-T	***	metabolic process

*/f
