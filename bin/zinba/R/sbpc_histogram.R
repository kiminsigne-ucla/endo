sbpc.histogram=function(inputfile,output,twoBitFile,max=500,min=0,numBin=500){
    if(!file.exists(inputfile)){
	stop(paste("Input file not found,",inputfile,sep=" "))
    }
    if(!file.exists(twoBitFile)){
	stop(paste("twoBit file not found,",twoBitFile,sep=" "))
    }
    .C("sbpchistogram",as.character(inputfile),as.character(output),as.character(twoBitFile),as.integer(max),as.integer(min),as.integer(numBin),PACKAGE="zinba")
    if(file.exists(output)){
	hdata <- read.table(output,header=T)
	jout <- paste(output,".jpeg",sep="")
	jpeg(jout,width=960,height=960)
	plot(hdata$BIN_STOP,hdata$COUNT,type="l",lwd=2)
	dev.off()
    }
}
