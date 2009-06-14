`f.Noise` <-
function(r1, fact=4.5, p=0.68, typ="percentile") {
		
		noise = 0
		r1<-(r1-median(r1, na.rm=TRUE))		
		
		if (typ == "percentile") { 
		r1 = abs(r1)
 		noise <- quantile(r1, probs=p, na.rm=TRUE)
		}		
		
		if (typ == "derivative") {
		noise <- dLRs(r1)
		}		
		
		if (typ == "combined") { 	
		noise1 <- dLRs(r1)
		r1 = abs(r1)
 		noise2 <- quantile(r1, probs=p, na.rm=TRUE)
		noise = (noise1 + noise2) / 2
		}
		
		t <- noise * fact
		
return(as.numeric(t))
}

