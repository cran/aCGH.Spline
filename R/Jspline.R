`Jspline` <-
function(x, offset=5, knots=1000, ntyp="percentile", p=0.68, fact=4.5, robust=TRUE, segN=FALSE, sn=0.75) {
	
	# Set java class path and heap size
	SetJ()
	gc(reset=TRUE)
	
	# Make sure no negative values
	x[x[,2]<1,2] = 1
	x[x[,3]<1,3] = 1
	
	if (offset < 1) {
	print("Offset cannot be less than one - settings offset to one!")
	offset=1
	}
	
	print("Calculating noise...")
	
	cleandata <- na.omit(x)
	r1 = log2(cleandata[,2] / cleandata[,3])
	r1 = r1[order(cleandata[,4], cleandata[,5])]
	t = 0;

	if (segN == TRUE) {
		r2 = r1 - segN(r1)
		t = f.Noise(r2, fact, p, ntyp)
		rm(r2)
	}
	
	if (segN == FALSE) {
		t = f.Noise(r1, fact, p, ntyp)
	}
	
	print("Performing Spline fitting and interpolation...")

	spVals = c(offset, knots)
	.jcall("Jspline", "[D", "SetValues", spVals)

	if (robust==FALSE) {
	ra <- cleandata[,c(2, 3)]
	Qiindex <- cleandata[,-c(1, 2, 3)]
	arg <- 0
	rm(cleandata, r1)
	gc(reset=TRUE)
	}
	
	if (robust==TRUE) {
	r1 = abs(r1 - median(r1, na.rm=TRUE))
	rxy <- cleandata[r1>t,c(2, 3)]
	ra <- cleandata[r1<t,c(2, 3)]
	Qind = cleandata[,-c(1, 2, 3)]
	Qiindex <- rbind(Qind[r1<t,], Qind[r1>t,])				 	 		
	arg <- length(rxy[,1])	
	fxy <- rxy
	rm(cleandata, r1, Qind)
	gc(reset=TRUE)
	}

	# run Jspline - Jspline.java
	fitR = .jcall("Jspline","[D", "RunSpline",ra[,1],ra[,2])
	fitG = .jcall("Jspline","[D", "RunSpline",ra[,2],ra[,1])
	fit = cbind(fitR, fitG)
	rm(fitR, fitG)
	gc(reset=TRUE)
	
	if (arg > 0) {	
	print(paste("Threshold = ", t))
	print(paste("Data points = ", length(rxy[,1])))
	
	# run interpolation/extrapolation
	fxy[,1] = .jcall("Jspline","[D", "INTERPOL2",rxy[,1],ra[,1],fit[,1])
	fxy[,2] = .jcall("Jspline","[D", "INTERPOL2",rxy[,2],ra[,2],fit[,2])
	
	fxx1 = c(fit[,1], fxy[,1])
	fxx2 = c(fit[,2], fxy[,2])
	fxx = cbind(fxx1, fxx2)
	tt = cbind(fxx, Qiindex)
	rm(ra, rxy, fxx1, fxx2, fxx, fxy)
	gc(reset=TRUE)
	}

	if (arg < 1) {		
	tt = cbind(fit, Qiindex)
	rm(ra, fit, Qiindex)
	gc(reset=TRUE)
	}

	
	# exclude the top and bottom 0.001 percentile (paired intensity values)
	p = mean(quantile(tt[,1], probs=0.999, na.rm=TRUE), quantile(tt[,2], probs=0.999, na.rm=TRUE))
	pp = mean(quantile(tt[,1], probs=0.001, na.rm=TRUE), quantile(tt[,2], probs=0.001, na.rm=TRUE))
	
	tt[tt[,1] >p & tt[,2]>p,1:2] = NA
	tt[tt[,1] <pp & tt[,2]<pp,1:2] = NA
	tt[tt[,1] >p & tt[,2]>p,7] = 1
	tt[tt[,1] <pp & tt[,2]<pp,7] = 1

	# Recombine data
	newresult <- matrix(ncol=length(x[1,]), nrow=length(x[,1]))
	colnames(newresult) <- colnames(x)
	tt = cbind(log2(tt[,1]/tt[,2]), tt)
	newresult <- rbind(tt, x[!complete.cases(x),])
	rm(tt)
	newresult[,1] <- log2(newresult[,2]/newresult[,3])
	final <- newresult[order(newresult[,7]),]
	rm(newresult)
	gc(reset=TRUE)

return(final)
}

