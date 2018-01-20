coord.sbpc=function(inputfile,coordfile,outputfile,twobitfile){
    if(!file.exists(inputfile)){
    	stop(paste("Input file not found,",inputfile,sep=" "))
    }
    if(!file.exists(coordfile)){
    	stop(paste("Coord file not found,",coordfile,sep=" "))
    }
    if(!file.exists(twobitfile)){
        stop(paste("twoBit file not found,",twobitfile,sep=" "))
    }
    if(is.null(outputfile)){
    	stop(paste("Need to specify an outputfile,",outputfile,sep=" "))
    }
    cReturn <- .C("coordSeqCount",as.character(inputfile),as.character(coordfile),as.character(outputfile),as.character(twobitfile),PACKAGE="zinba")
}
