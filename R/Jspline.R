`Jspline` <-
function(x, offset=5, knots=1000, ntyp="percentile", p=0.68, fact=4.5, robust=TRUE, segN=FALSE, sn=0.75, index1=1, index2=2) {
	
	# Set java class path and heap size
	SetJ()
	
	# Make sure no negative values
	x[x[,1]<1,1] = 1
	x[x[,2]<1,2] = 1
	
	if (offset < 1) {
	print("Offset cannot be less than one - settings offset to one!")
	offset=1
	}
	
	print("Calculating noise...")
	
	cleandata <- na.omit(x)
	r1 = log2(cleandata[,1] / cleandata[,2])
	r2 = r1[order(cleandata[,3], cleandata[,4])]
	t = 0;

	if (segN == TRUE) {
		r3 = segN(r2, sn)
		t = f.Noise(r3, fact, p, ntyp)
	}
	
	if (segN == FALSE) {
	t = f.Noise(r2, fact, p, ntyp)
	}
	
	print("Performing Spline fitting and interpolation...")

	spVals = c(offset, knots)
	.jcall("Jspline", "[D", "SetValues", spVals)

	if (robust==FALSE) {
	ra <- cleandata[,c(index1, index2)]
	Qiindex <- cleandata[,-c(index1, index2)]
	arg <- 0
	}
	
	if (robust==TRUE) {
	r1 = abs(r1 - median(r1, na.rm=TRUE))
	rxy <- cleandata[r1>t,c(index1, index2)]
	ra <- cleandata[r1<t,c(index1, index2)]
	Qind = cleandata[,-c(index1, index2)]
	Qiindex <- rbind(Qind[r1<t,], Qind[r1>t,])				 	 		
	arg <- length(rxy[,1])	
	fxy <- rxy
	}
	
	# run Jspline - Jspline.java
	fitR = .jcall("Jspline","[D", "RunSpline",ra[,1],ra[,2])
	fitG = .jcall("Jspline","[D", "RunSpline",ra[,2],ra[,1])
	fit = cbind(fitR, fitG)
	
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
	}
	
	if (arg < 1) {		
	tt = cbind(fit, Qiindex)	
	}
	
	# exclude the top and bottom 0.001 percentile (paired intensity values)
	p = mean(quantile(tt[,1], probs=0.999, na.rm=TRUE), quantile(tt[,2], probs=0.999, na.rm=TRUE))
	pp = mean(quantile(tt[,1], probs=0.001, na.rm=TRUE), quantile(tt[,2], probs=0.001, na.rm=TRUE))
	
	tt[tt[,1] >p & tt[,2]>p,1:2] = NA
	tt[tt[,1] <pp & tt[,2]<pp,1:2] = NA
	tt[tt[,1] >p & tt[,2]>p,7] = 1
	tt[tt[,1] <pp & tt[,2]<pp,7] = 1

	# Recombine data
	colnames(tt) <- colnames(x)
	newresult <- matrix(ncol=length(x[1,]), nrow=length(x[,1]))
	colnames(newresult) <- colnames(x)
	newresult <- rbind(tt, x[!complete.cases(x),]) 
	#final <- newresult[order(newresult[,3], newresult[,4]),]
	final <- newresult[order(newresult[,6]),]

return(final)
}

