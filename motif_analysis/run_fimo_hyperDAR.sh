PROJ=/mnt/c/Users/Yamanto/Desktop/dtra_arabidopsis_manuscript/at_rad_motif_analysis
MOTIF1=/mnt/e/meme/motif_databases/ARABD/ArabidopsisPBM_20140210.meme
MOTIF2=/mnt/e/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme


for IN1 in up_dar_01_72 up_dar_1_72 
do
	        echo ${IN1}

		        FASTAIN=$PROJ/dreme_sequences/${IN1}.fasta.masked
			        STREAMOUT=$PROJ/streme_analysis/${IN1}/streme.txt

				        FIMOOUT1=$PROJ/fimo_output_PBM/${IN1}
					        mkdir -p $FIMOOUT1

						        FIMOOUT2=$PROJ/fimo_output_DAP/${IN1}
							        mkdir -p $FIMOOUT2

								        TOMOUT1=$PROJ/tomtom_output_PBM/${IN1}
									        mkdir -p $TOMOUT1

										        TOMOUT2=$PROJ/tomtom_output_DAP/${IN1}
											        mkdir -p $TOMOUT2

												        echo "${IN1} FIMO"
												fimo --oc $FIMOOUT1 $MOTIF1 $FASTAIN
												#fimo --oc $FIMOOUT2 $MOTIF2 $FASTAIN
												echo "${IN1} Tomtom"
												tomtom --oc $TOMOUT1 $STREMEOUT $MOTIF1
												#tomtom --oc $TOMOUT2 $STREMEOUT $MOTIF2

											done
