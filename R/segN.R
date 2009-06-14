`segN` <-
function(r, t) {
d = abs(diff(r))
n = quantile(d, probs=t, na.rm=TRUE)
c =which(d>n, arr.ind=TRUE)
rr = r	

	if (c[1] != 1) {
	rr[1:c[1]] = r[1:c[1]] - median(r[1:c[1]], na.rm=TRUE)
	}
	
		for (x in 1:length(c)-1) {
		i = x+1
		rr[c[x:i]] = r[c[x:i]] - median(r[c[x:i]], na.rm=TRUE)
		}
	
	if (c[length(c)] != length(r)) {		
	rr[c[length(c)]:length(r)] = r[c[length(c)]:length(r)] - median(r[c[length(c)]:length(r)], na.rm=TRUE)	
	}

return (rr)
}

