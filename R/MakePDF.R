`MakePDF` <-
function(filename, x, raw, ntyp="percentile", p=0.68, fact=4.5, segN=FALSE, sn=0.75) {
	print(paste("Writing results to - ", filename, "...", sep=""))
	pdf(file=filename)
	ProPlot(x)
	MAPlot1(x)
	MAPlot2(x, raw, ntyp, p, fact, segN, sn)
	GenomePlot(x)
	dev.off()	
}

