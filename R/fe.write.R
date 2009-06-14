`fe.write` <-
function(file, x) {
	SetJ()
	.jcall("Jspline", "V", "SetOutput", file)
	x = x[order(x[,6]),]	
    moi = .jcall("Jspline", "[D", "OutAll", x[,1], x[,2], x[,7])
}

