buildcustomwindowdata=function(seq,input="none",twoBit,winSize=500,offset=0,filelist,filetype="bowtie",extension, outdir="default"){

	if(!file.exists(seq)){
		stop(paste("Seq file not found,",seq,sep=" "))
	}
	if(input != "none" && !file.exists(input)){
		stop(paste("Input file not found,",input,sep=" "))
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
	if(is.null(filelist)){
		stop(paste("filelist was not specified"))
	}
	if(!file.exists(outdir) && outdir!="default"){
		stop(paste("buildwindowdata: pecified output directory to place built datafiles doesnt exist"))
	}
	cReturn <- .C("buildCustomWindows",as.character(seq),as.character(input),as.character(twoBit),as.integer(winSize),as.integer(offset),as.character(filetype),as.character(filelist),as.integer(extension), as.character(outdir),PACKAGE="zinba")
	gc()
}
