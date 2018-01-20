getwindowcounts=function(seq,twoBit,winSize=500,offset=0,filetype="bowtie",extension, Nthresh=.1){	
	if(!file.exists(seq)){
		stop(paste("Seq file not found,",seq,sep=" "))
	}
	if(!file.exists(twoBit)){
		stop(paste("twoBit file not found,",twoBit,sep=" "))
	}
	if(filetype != "bowtie" && filetype != "tagAlign" && filetype != "bed"){
		stop(paste("Incorrect filetype:",filetype,"[bowtie|tagAlign|bed]",sep=" "))
	}
	if(is.null(extension)){
		stop(paste("extension was not specified"))
	}
	cReturn <- .C("getWindowCounts",as.character(seq),as.character(twoBit),as.integer(winSize),as.integer(offset),as.character(filetype),as.integer(extension), as.double(Nthresh),PACKAGE="zinba")
	gc()
}
