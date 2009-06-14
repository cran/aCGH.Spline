`WriteQC` <-
function(filename, batch, x) {
	
	Xmed = median(x[x[,4] == 23,1], na.rm=TRUE)
	dL = dLRs(x[order(x[,4], x[,5]),1])
	l = x[x[,4]<23,1]
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

