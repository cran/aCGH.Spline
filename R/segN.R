`segN` <-
function(ddd) {
	SetJ()
return (as.numeric(.jcall("Jspline", "[D", "runSegN", ddd)))
}

