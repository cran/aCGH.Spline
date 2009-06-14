`MAPlot2` <-
function(x, raw, ntyp="percentile", p=0.68, fact=4.5, segN=FALSE, sn=0.75, pch=46, Fitcol="green", Intcol="red", ylim=c(-5,5), xlim=c(0,20), xlab="log2(mean(cy5, cy3))", ylab="log2(cy5/cy3)", main="MA-plot") {
	
	raw = raw[order(raw[,3], raw[,4]),]
	cleandata <- na.omit(raw)
	r1 = log2(cleandata[,1] / cleandata[,2])
	t = 0;

	if (segN == TRUE) {
		r2 = r1[order(cleandata[,3], cleandata[,4])]
		r3 = segN(r2, sn)
		t = f.Noise(r3, fact, p, ntyp)
	}
	
	if (segN == FALSE) {
	t = f.Noise(r1, fact, p, ntyp)
	}
    
	a = abs(log2(x[,1] / x[,2]))
	intes = cbind(x[,1],x[,2], a)
	plot(log2(apply(intes[intes[,3] < t,1:2], 1, mean)), log2(intes[intes[,3] < t,1] / intes[intes[,3] < t,2]), ylim=ylim, xlim=xlim, pch=pch, col=Fitcol, xlab=xlab, ylab=ylab, main=main)
	matplot(log2(apply(intes[intes[,3] > t,1:2], 1, mean)), log2(intes[intes[,3] > t,1] / intes[intes[,3] > t,2]), pch=pch, add=TRUE, col=Intcol)
	abline(h=0, col="red")
	legend("topleft", col=c("green", "red"), pch=pch, c("Fitted", "Interpolated"))
 }

