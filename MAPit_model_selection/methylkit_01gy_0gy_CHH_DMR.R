library(methylKit)
library(genomation)

file.list=list("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r1_CHH.txt",
	                      "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r2_CHH.txt",
			                     "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r1_CHH.txt",
			                     "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_01gy_r2_CHH.txt")


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
	       sample.id=list("ctrl1","ctrl2","test1","test2"),
	       assembly="tair10",
	       treatment=c(0,0,1,1),
	       context="CHH")

myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CHH_10cGy.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = regionCounts(myobj,myRanges)	


meth=unite(myRegionData, destrand=FALSE)

myDiff=calculateDiffMeth(meth)

myDiffOverDispers <- calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_01gy_CHH.DMR.results.txt",
	                append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_01gy_CHH.DMR.results.overdispersion.txt",
	                append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)
