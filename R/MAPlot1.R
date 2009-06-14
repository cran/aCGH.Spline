`MAPlot1` <-
function(x, pch=46, Retcol="black", Excol="red", ylim=c(-5,5), xlim=c(0,20), xlab="log2(mean(cy5, cy3))", ylab="log2(cy5/cy3)", main="MA-plot") {
	intes = cbind(x[,1], x[,2], x[,7])
	plot(log2(apply(intes[intes[,3]==0,1:2], 1, mean)), log2(intes[intes[,3]==0,1] / intes[intes[,3]==0,2]), ylim=ylim, xlim=xlim, pch=pch, xlab=xlab, ylab=ylab, main=main, col=Retcol)
	matplot(log2(apply(intes[intes[,3]==1,1:2], 1, mean)), log2(intes[intes[,3]==1,1] / intes[intes[,3]==1,2]), pch=pch, add=TRUE, col=Excol)
	abline(h=0, col="red")
	legend("topleft", col=c(Retcol, Excol), pch=pch, c("Retained", "Excluded"))
}

