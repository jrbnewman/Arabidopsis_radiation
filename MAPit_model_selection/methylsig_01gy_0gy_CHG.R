
library(bsseq)
library(methylSig)


files=c("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r1_CHG.txt",
	"/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0gy_r2_CHG.txt",
	"/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r1_CHG.txt",
	"/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_01gy_r2_CHG.txt")

bs = bsseq::read.bismark(files=files,
			 colData = data.frame(row.names = c("ctrl1","ctrl2","test1","test2"), type=c("control","control","case","case")),
			 rmZeroCov = FALSE,
			 strandCollapse = FALSE)

bsseq::pData(bs)

binTest = diff_binomial(bs=bs, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'))
betaTest = diff_methylsig(bs=bs, group_column='type', comparison_groups = c('case' = 'case', 'control' = 'control'),
			  disp_groups = c('case' = TRUE, 'control' = TRUE),
			  local_window_size = 0,
			  t_approx = TRUE,
			  n_cores = 1)


write.table(binTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHG.binomial.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(betaTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylsig_0cGy_01cGy_CHG.betabinomial.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)



