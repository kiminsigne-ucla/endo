buildwindowdata2=function(
  seq,
  input="none",
  twoBit,
  winSize=500,
  offset=0,
  filelist,
  filetype="bowtie",
  extension, 
  cnvWinSize=100000,
  cnvOffset=0,
  outdir="default"
){

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
		stop(paste("buildwindowdata: specified output directory to place built datafiles doesnt exist"))
		if(!isDirectory(outdir)) stop("specified output directory path is not a directory")
	}
        if(substr(x=outdir, nchar(outdir), nchar(outdir)) !="/"){
		outdir = paste(outdir, "/", sep="")
        }
	
  cReturn <- .C(
    "buildWindowsfast",
    as.character(seq),
    as.character(input),
    #as.character(align),
    as.character(twoBit),
    as.integer(winSize),
    as.integer(offset),
    as.integer(cnvWinSize),
    as.integer(cnvOffset),
    as.character(filetype),
    as.character(filelist),
    as.integer(extension), 
    as.character(outdir),
    #as.integer(binary),
    PACKAGE="zinba"
  )
    

}
