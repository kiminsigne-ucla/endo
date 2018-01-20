alignAdjust=function(inputfile, outDir="",twoBitFile,athresh=1,adjustsize=0){
    if(!file.exists(inputfile)){
	stop(paste("Input file not found,",inputfile,sep=" "))
    }
    if(!file.exists(twoBitFile)){
	stop(paste("twoBit file not found,",twoBitFile,sep=" "))
    }
    .C("alignAdjust",as.character(inputfile),as.character(outDir),as.character(twoBitFile),as.integer(athresh),as.integer(adjustsize),PACKAGE="zinba")
    gc()
}
