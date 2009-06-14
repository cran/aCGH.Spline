`WriteGFF` <-
function(filename, x) {
	print(paste("Writing results to - ", filename, "...", sep=""))
	sorted <- x[order(x[,4], x[,5]),]
	GFF = cbind(paste("chr", sorted[sorted[,8]==0,4], sep=""), "DNAcopy", filename, sorted[sorted[,8]==0,5], sorted[sorted[,8]==0,6], sorted[sorted[,8]==0,1], ".", ".", "; color 000000")
	write.table(GFF, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
}

