`fe.write` <-
function(file, x) {
	SetJ()
	.jcall("Jspline", "V", "SetOutput", file)	
    moi = .jcall("Jspline", "[D", "OutAll", x[,1], x[,2], x[,7])
}

