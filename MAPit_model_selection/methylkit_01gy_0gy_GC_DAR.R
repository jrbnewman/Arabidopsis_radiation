### running methylKit

library(methylKit)
library(genomation)
file.list=list("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_0U_R1_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_0U_R2_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_100U_R1_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_100U_R2_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_0U_R1_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_0U_R2_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_100U_R1_v01.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_100U_R2_v01.txt")


covariates = data.frame(units=c(0,0,1,1,0,0,1,1), treat_units=c(0,0,0,0,0,0,1,1))

# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
	       sample.id=list("ctrl1_0U","ctrl2_0U","ctrl1_100","ctrl2_100","test1_0U","test2_0U","test1_100U","test2_100U"),
	       assembly="tair10",
	       treatment=c(0,0,0,0,1,1,1,1),
	       context="GC")
meth=unite(myobj, destrand=FALSE)

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DAR_GC_10cGy.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = regionCounts(myobj,myRanges)	

meth=unite(myRegionData, destrand=FALSE)

myDiff=calculateDiffMeth(meth, covariates=covariates)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0cGy_01cGy_GC.DAR.results2.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

myDiffOverDispers <- calculateDiffMeth(meth, covariates=covariates, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0cGy_01cGy_GC.DAR.results2.overdispersion.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)



