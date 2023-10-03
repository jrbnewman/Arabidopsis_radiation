
library(methylKit)
file.list=list("/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r1_CHH.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0gy_r2_CHH.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r1_CHH.txt",
	       "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_1gy_r2_CHH.txt")


# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
	       sample.id=list("ctrl1","ctrl2","test1","test2"),
	       assembly="tair10",
	       treatment=c(0,0,1,1),
	       context="CHH")
meth=unite(myobj, destrand=FALSE)


myDiff=calculateDiffMeth(meth)
# get hyper methylated bases
#myDiff25p.hyper=getMethylDiff(myDiff,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
#myDiff25p.hypo=getMethylDiff(myDiff,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
#myDiff25p=getMethylDiff(myDiff,difference=10,qvalue=0.01)

myDiffOverDispers <- calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores=1)

write.table(myDiff, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0cGy_1cGy_CHH.results.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)

write.table(myDiffOverDispers, "/TB14/TB14/sandbox/dtra_sandbox/at_rad_methylkit_0cGy_1cGy_CHH.results.overdispersion.txt", append = FALSE, sep = "\t", dec = ".",   row.names = FALSE, col.names = TRUE)




