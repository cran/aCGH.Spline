`ProPlot` <-
function(x, pch=46, col="black", ylim=c(-3,3), xlim=c(0,length(x[,1])), xlab="Index", ylab="log2 ratio", main="Profile Plot") {
	sortR <-x[order(x[,3], x[,4]),1]
	sortG <- x[order(x[,3], x[,4]),2]  
	plot(log2(sortR / sortG), ylim=ylim, xlim=xlim, pch=pch, col=col, xlab=xlab, ylab=ylab, main=main)
	abline(h=0, col="red")
	abline(h=1, col="blue")
	abline(h= -1, col="blue")	
}

