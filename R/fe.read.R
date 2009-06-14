`fe.read` <-
function(file, Raw=TRUE) {
	
	SetJ()
    gc(reset=TRUE)
	options(warn=-1)	# set to remove warnings refering to NA values (which are set to NA deliberatly)
	raw = as.character(Raw)
	size = .jcall("Jspline", "I","InAll", file, raw)
	dat = matrix(ncol=8, nrow=size)
	dat[,1] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "log"))
	dat[,2] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "red"))
	dat[,3] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "green"))
	dat[,4] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "chr"))
	dat[,5] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "start"))
	dat[,6] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "stop"))
	dat[,7] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "index"))
	dat[,8] = as.numeric(.jcall("Jspline", "[S", "ReturnHelper", "flag"))
	
	.jcall("Jspline", "V", "FreeMem")
	
return(dat)
}

