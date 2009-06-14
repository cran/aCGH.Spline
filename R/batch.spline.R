`batch.spline` <-
function	(dir,
							format="FE", #c("NIM")
							raw = TRUE, 
							robust=TRUE, 
							offset =5,
							knots = 1000,
							ntyp="percentile", # c("derivative", "combined")
							p=0.68, 
							fact=4.5,
							segN = FALSE,
							sn = 0.75,
							writeFE = FALSE,
							QC=FALSE, 
							GFF=FALSE, 
							PDF=FALSE) {
	
	
	slash = substr(dir,nchar(dir),nchar(dir))
	if (slash != "/") {
	dir = paste(dir, "/", sep="")
	}	
	
	names <- dir(dir)

	index <- grep(".*txt$", names)

	if (QC==TRUE) { qc = matrix(ncol=4, nrow=length(index)) }

	for (x in 1 : length(index)) 
	{
	
	pin <- index[x]
	file <- paste(dir, names[pin], sep="")
	
	if (format == "FE") { print("Using fe.read (Agilent feture extraction files)"); print(paste("Reading ", file, "....", sep="")); data <- fe.read(file, Raw=raw); } 
	if (format == "NIM") { print("Using nim.read (Nimblegen seg_MNT.txt files)"); print(paste("Reading ", file, "....", sep="")); data <- nim.read(file, Raw=raw); }
	
	final <- Jspline(data, offset, knots, ntyp, p, fact, robust, segN, sn)	
    justname = substr(file, 0, nchar(file) -4)
    tname = paste(justname, "_processed.temp", sep="")
	pname = paste(justname, ".pdf", sep="")
	gname = paste(justname, ".gff", sep="")
	print(paste("Writing results to - ", tname, "....", sep=""))
	
	write.table(final, file=tname, row.names=FALSE, sep="\t", quote=FALSE, col.names=FALSE)
		
		if (writeFE == TRUE) {
		fe.write(file, final)
		}
		
		if (PDF==TRUE) {
		MakePDF(pname, final, data, ntyp, p, fact, segN, sn)
		} 
		
		if (GFF==TRUE) {
		WriteGFF(gname, final)
		}
		
		if (QC==TRUE) {
		qc[x,] = WriteQC(file, batch=TRUE, final)
		}
		
	print("")	
	
	}
	
	if (QC==TRUE) { 
	print("Writing QC report")
	colnames(qc) = c("Name", "rp68", "dLRs", "Xmedian")
	write.table(qc, file="QC_values.txt", row.names=FALSE, sep="\t", quote=FALSE) 
	}
	
	print("Finished!!!")
}

