`nim.read` <-
function(file) {

data <- read.table(file, header=TRUE)

	format = matrix(nrow=length(data[,1]), ncol=7)
	red = 2^data$EXP_SPATIAL
	green = 2^data$REF_SPATIAL
	check = red + green

	loc = as.vector(data$CHROMOSOME)
	flag = vector(length=length(data[,1]))
	flag = 0;

		for (x in 1:length(loc)) {
		format[x,3] = substr(loc[x],4, nchar(loc[x]))
		}

	format[,1] = red
	format[,2] = green
	format[,4] = data$POSITION
	format[,5] = data$POSITION
	format[,6] = data$INDEX
	format[,7] = flag

	format[format[,3]=="X",3] = 23
	format[format[,3]=="Y",3] = 24

	format[check < 100,7] = 1
	format[red > 60000 | green > 60000 | red < 0 | green < 0, 1:2] = NA
	format[red > 60000 | green > 60000 | red < 0 | green < 0, 7] = 1

	format2 = matrix(ncol=length(format[1,]), nrow=length(format[,1]))
		
		for (x in 1:length(format[1,])) {
		format2[,x] = as.numeric(format[,x])
		}

return(format2)
}

