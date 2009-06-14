`fe.extract` <-
function(file, list) {

head = scan(file, what="", skip=9, nlines=1)
ind_list = vector(length=length(list))

	for (x in 1:length(list)) {
	ind_list[x] = grep(list[x], head)
	}

lines = scan(file, what="", skip=10, sep="\n")
iLines = grep("chr", lines)
nLines = lines[iLines]
newLines =  strsplit(nLines, "\t", extended = TRUE, fixed = FALSE, perl = FALSE)

format = matrix(ncol=length(list), nrow=length(nLines))
 
	for (x in 1:length(nLines)) {
		for (i in 1:length(list)) {
		t = newLines[x]
		tt = unlist(t)
		format[x,i] = tt[ind_list[i]]
		}
	}
	
return (data.frame(format))
}

