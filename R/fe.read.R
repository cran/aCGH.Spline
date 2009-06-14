`fe.read` <-
function(file, raw=TRUE) {
	
	SetJ()
	options(warn=-1)	# set to remove warnings refering to NA values (which are set to NA deliberatly)
	raw = as.character(raw)
	size = .jcall("Jspline", "I","InAll", file, raw)
	dat = matrix(ncol=7, nrow=size)
	dat[,1] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "red"))
	dat[,2] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "green"))
	dat[,3] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "chr"))
	dat[,4] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "start"))
	dat[,5] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "stop"))
	dat[,6] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "index"))
	dat[,7] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "flag"))
	
	.jcall("Jspline", "V", "FreeMem")
	
return(dat)
}

