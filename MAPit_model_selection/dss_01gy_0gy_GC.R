
library(DSS)
library(bsseq)

treatment = c("ctrl", "ctrl", "ctrl", "ctrl", "test", "test", "test", "test")
units = c("0U", "100U", "0U", "100U", "0U", "100U", "0U", "100U")
design = data.frame(treatment,units)
design

dat1.1A <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_0U_R1.txt", header=TRUE)
dat1.1B <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_100U_R1.txt", header=TRUE)
dat1.2A <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_0U_R2.txt", header=TRUE)
dat1.2B <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0gy_100U_R2.txt", header=TRUE)
dat2.1A <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_0U_R1.txt", header=TRUE)
dat2.1B <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_100U_R1.txt", header=TRUE)
dat2.2A <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_0U_R2.txt", header=TRUE)
dat2.2B <- read.table("/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_01gy_100U_R2.txt", header=TRUE)

BSobj = makeBSseqData( list(dat1.1A, dat1.1B, dat1.2A, dat1.2B, dat2.1A, dat2.1B, dat2.2A, dat2.2B), c("C1A","C1B","C2A","C2B", "T1A","T1B", "T2A","T2B") )

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~treatment+units+treatment:units)

DMLtest.treatment = DMLtest.multiFactor(DMLfit, coef="treatmenttest")
DMLtest.int = DMLtest.multiFactor(DMLfit, coef="treatmenttest:units100U")

write.table(DMLtest.treatment, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0cGy_01cGy_GC.treatment.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(DMLtest.int, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_dss_0cGy_01cGy_GC.trt_by_units.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)



