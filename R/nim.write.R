`nim.write` <-
function(file, x) {

	SetJ()
	gc(reset=TRUE)
	.jcall("Jspline", "V", "SetOutput", file)
	x = x[order(x[,7]),]	
    moi = .jcall("Jspline", "[D", "OutAllNimbelgen", x[,1], x[,2], x[,3], x[,8])

}

