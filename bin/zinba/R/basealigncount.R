basealigncount=function(inputfile,outputfile,twoBitFile,extension=NULL,filetype="bowtie", binary=0){
		library(R.utils)
    if(!file.exists(inputfile)){
    	stop(paste("Input file not found,",inputfile,sep=" "))
    }else{
			inputfile=getAbsolutePath(inputfile)
		}
    if(!file.exists(twoBitFile)){
        stop(paste("twoBit file not found,",twoBitFile,sep=" "))
    }
    if(filetype != "bowtie" && filetype != "tagAlign" && filetype != "bed"){
        stop(paste("Incorrect filetype:",filetype,"[bowtie|tagAlign|bed]",sep=" "))
    }
    if(is.null(extension)){
        stop(paste("Need to specify an extension length",,sep=" "))
    }
		#if(is.null(getParent(outputfile))){
			#
		#}else if(!file.exists(getParent(outputfile))){
		#	stop("Directory where output is being sent doesn't exist")
		#}else{
		#	outputfile=getAbsolutePath(outputfile)
		#}


    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(outputfile), as.character(twoBitFile),as.integer(extension),as.character(filetype),as.integer(binary), PACKAGE="zinba")
		if(binary==1){
			library(R.utils)
			cat("Because binary option is selected, appending '.bin.gz' to output path.\n")
			newpath=paste(outputfile, ".bin.gz", sep="")
			tar(tarfile=newpath, files=dir(path=getParent(outputfile), pattern=paste(outputfile,"_", sep="")), compression="gzip")
			cat("Compressed binary basecount file can be found at", newpath, "\n")
		}
}
