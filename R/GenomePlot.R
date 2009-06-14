`GenomePlot` <-
function(x, ylim=c(-3,3), axes=FALSE, pch=46, col="black") {
	x <- x[order(x[,4], x[,5]),]
	uni= unique(x[,4])
	uni = uni[complete.cases(uni)]
	si = ceiling(length(uni) / 4)
	par(mfrow=c(si,4), mar=c(1,1,1,1)) 
	
	for (i in 1:length(uni)) {
	pin = uni[i]
	plot(x[x[,4]==pin,4], log2(x[x[,4]==pin,2]/x[x[,4]==pin,3]),pch=pch, col=col, ylim=ylim, xlim=c(min(x[x[,4]==pin,5], na.rm=TRUE),max(x[x[,4]==pin,5], na.rm=TRUE)), axes=axes, main=paste("chr", pin, sep=""))
	abline(h=0, col="red")
	}
}

