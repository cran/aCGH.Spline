`nim.read` <-
function(file, Raw=TRUE) {
	
	gc(reset=TRUE)
 	data <- read.table(file, header = TRUE)
 	options(warn=-1)	# set to remove warnings refering to NA values (which are set to NA deliberatly)
    format = matrix(nrow = length(data[, 1]), ncol = 8)
    red = 0
    green = 0
    if (Raw == TRUE) {
        red = 2^data$EXP_SPATIAL
        green = 2^data$REF_SPATIAL
    }
    if (Raw == FALSE) {
        red = 2^data$EXP_NORM
        green = 2^data$REF_NORM
    }
    check = red + green
    loc = data$CHROMOSOME
    ind = vector(length=length(loc))
    ind[grep("chr",loc)] = TRUE
    loc[ind==FALSE] = NA
	loc = as.character(loc) 
	loc[loc == "chrX"] = "chr23"
    loc[loc == "chrY"] = "chr24"
	
    flag = vector(length = length(data[, 1]))
    flag = 0
    format[, 1] = as.numeric(log2(red/green))
    format[, 2] = as.numeric(red)
    format[, 3] = as.numeric(green)
    format[, 4] = as.numeric(substr(loc, 4, 1000))
    format[, 5] = as.numeric(data$POSITION)
    format[, 6] = as.numeric(data$POSITION)
    format[, 7] = as.numeric(data$INDEX)
    format[, 8] = as.numeric(flag)
    rm(data, loc, ind)
    gc()
    format[check < 100, 8] = 1
    rm(check)
    gc()
	format[is.na(format[,4]), 1:3] = NA
	format[is.na(format[,4]), 8] = 1
    format[red > 60000 | green > 60000 | red < 0 | green < 0, 1:3] = NA
    format[red > 60000 | green > 60000 | red < 0 | green < 0, 8] = 1
	rm(red, green)
	gc()
    return(format)
}

