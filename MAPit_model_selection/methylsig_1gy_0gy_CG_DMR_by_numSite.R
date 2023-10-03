library(bsseq)
library(methylSig)
library(genomation)

files=c("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r1_CG.txt",
	                "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r2_CG.txt",
			                        "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r1_CG.txt",
			                        "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_1gy_r2_CG.txt")

bs = bsseq::read.bismark(files=files,
			                          colData = data.frame(row.names = c("ctrl1","ctrl2","test1","test2"), type=c("control","control","case","case")),
						                           rmZeroCov = FALSE,
						                           strandCollapse = FALSE)


# 2 sites

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_100cGy_minSites_2.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = tile_by_regions(bs = bs, gr = myRanges)



bsseq::pData(myRegionData)

binTest = diff_binomial(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'))

write.table(binTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_2.binomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

betaTest = diff_methylsig(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'),
			                            disp_groups = c('case' = TRUE, 'control' = TRUE),
						                              local_window_size = 0,
						                              t_approx = TRUE,
									                                n_cores = 1)

write.table(betaTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_2.betabinomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)




# 2 sites

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_100cGy_minSites_3.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = tile_by_regions(bs = bs, gr = myRanges)



bsseq::pData(myRegionData)

binTest = diff_binomial(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'))

write.table(binTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_3.binomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

betaTest = diff_methylsig(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'),
			                            disp_groups = c('case' = TRUE, 'control' = TRUE),
						                              local_window_size = 0,
						                              t_approx = TRUE,
									                                n_cores = 1)

write.table(betaTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_3.betabinomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)




# 5 sites

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_100cGy_minSites_5.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = tile_by_regions(bs = bs, gr = myRanges)



bsseq::pData(myRegionData)

binTest = diff_binomial(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'))

write.table(binTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CG_minSites_5.binomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

betaTest = diff_methylsig(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'),
			                            disp_groups = c('case' = TRUE, 'control' = TRUE),
						                              local_window_size = 0,
						                              t_approx = TRUE,
									                                n_cores = 1)

write.table(betaTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_5.betabinomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)




# 10 sites

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CG_100cGy_minSites_10.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = tile_by_regions(bs = bs, gr = myRanges)



bsseq::pData(myRegionData)

binTest = diff_binomial(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'))

write.table(binTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CG_minSites_10.binomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

betaTest = diff_methylsig(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'),
			                            disp_groups = c('case' = TRUE, 'control' = TRUE),
						                              local_window_size = 0,
						                              t_approx = TRUE,
									                                n_cores = 1)

write.table(betaTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_1cGy_CG_minSites_10.betabinomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)
