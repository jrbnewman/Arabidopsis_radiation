library(DSS)
library(bsseq)

dat1.1 <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r1_CHH.txt", header=TRUE)
dat1.2 <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_r2_CHH.txt", header=TRUE)
dat2.1 <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r1_CHH.txt", header=TRUE)
dat2.2 <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_1gy_r2_CHH.txt", header=TRUE)
BSobj = makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2), c("C1","C2", "T1", "T2") )

dmlTest = DMLtest(BSobj, group1=c("T1","T2"), group2=c("C1","C2"))

write.table(dmlTest, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0cGy_1cGy_CHH.results.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)




