`WriteQC` <-
function(filename, batch, x) {
	
	Xmed = median(log2(x[x[,3] == 23,1] / x[x[,3] == 23,2]), na.rm=TRUE)
	dL = dLRs(log2(x[,1] / x[,2]))
	l = log2(x[x[,3]<23,1] / x[x[,3]<23,2])
	al = abs(l)
	rp = quantile(al, 0.68, na.rm=TRUE)
	qc = cbind(rp, dL, Xmed)
	
	if (batch==TRUE) { 
	qc = cbind(filename, qc)
	return(qc) 
	}
	
	if (batch==FALSE) { 
	colnames(qc) = c("68thp", "dLRs", "X_median")
	write.table(qc, file=paste(filename, "_QCvalues.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
	}	
}

