library(bsseq)
library(methylSig)
library(genomation)

files=c("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r1_CHG.txt",
	        "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r2_CHG.txt",
		        "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r1_CHG.txt",
		        "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r2_CHG.txt")

bs = bsseq::read.bismark(files=files,
			 colData = data.frame(row.names = c("ctrl1","ctrl2","test1","test2"), type=c("control","control","case","case")),
			 rmZeroCov = FALSE,
			 strandCollapse = FALSE)

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_10cGy.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = tile_by_regions(bs = bs, gr = myRanges)



bsseq::pData(myRegionData)

binTest = diff_binomial(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'))

write.table(binTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHG.binomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

betaTest = diff_methylsig(bs=myRegionData, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'),
			  disp_groups = c('case' = TRUE, 'control' = TRUE),
			  local_window_size = 0,
			  t_approx = TRUE,
			  n_cores = 1)

write.table(betaTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHG.betabinomial.DMR.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

