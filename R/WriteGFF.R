`WriteGFF` <-
function(filename, x) {
	print(paste("Writing results to - ", filename, "...", sep=""))
	sorted <- x[order(x[,3], x[,4]),]
	GFF = cbind(paste("chr", sorted[sorted[,6]==0,3], sep=""), "DNAcopy", filename, sorted[sorted[,6]==0,4], sorted[sorted[,6]==0,5], log2(sorted[sorted[,6]==0,1] / sorted[sorted[,6]==0,2]), ".", ".", "; color 000000")
	write.table(GFF, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
}

