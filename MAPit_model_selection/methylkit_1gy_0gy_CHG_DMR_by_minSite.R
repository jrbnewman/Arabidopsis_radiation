library(methylKit)
library(genomation)

file.list=list("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r1_CHG.txt",
	                                     "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r2_CHG.txt",
					                                                  "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r1_CHG.txt",
					                                                  "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r2_CHG.txt")


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
	                      sample.id=list("ctrl1","ctrl2","test1","test2"),
			                     assembly="tair10",
			                     treatment=c(0,0,1,1),
					                    context="CHG")


#min sites 2
myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_100cGy_minSites_2.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = regionCounts(myobj,myRanges)


meth=unite(myRegionData, destrand=FALSE)

myDiff=calculateDiffMeth(meth)

myDiffOverDispers <- calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_2.DMR.results.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_2.DMR.results.overdispersion.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

#min sites 3
myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_100cGy_minSites_3.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = regionCounts(myobj,myRanges)


meth=unite(myRegionData, destrand=FALSE)

myDiff=calculateDiffMeth(meth)

myDiffOverDispers <- calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_3.DMR.results.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_3.DMR.results.overdispersion.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)


#min sites 5
myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_100cGy_minSites_5.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = regionCounts(myobj,myRanges)


meth=unite(myRegionData, destrand=FALSE)

myDiff=calculateDiffMeth(meth)

myDiffOverDispers <- calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_5.DMR.results.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_5.DMR.results.overdispersion.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)



#min sites 10
myBed = "/TB14/TB14/sandbox/dtra_sandbox/DMR_CHG_100cGy_minSites_10.bed"
myRanges = readBed(myBed, track.line=FALSE, remove.unusual=FALSE)
myRegionData = regionCounts(myobj,myRanges)


meth=unite(myRegionData, destrand=FALSE)

myDiff=calculateDiffMeth(meth)

myDiffOverDispers <- calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_10.DMR.results.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_1gy_CHG_minSites_10.DMR.results.overdispersion.txt",
	                            append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)



