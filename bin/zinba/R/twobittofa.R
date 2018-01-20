twobittofa=function(chrm,start=NULL,end=NULL,twoBitFile,gcSeq=NULL, outdir=" "){
	if(outdir!=" " & !file.exists(outdir)){
		stop("Specified directory 'outdir' doesnt exist or is not correct")
	}else{
		if(strsplit(outdir,"")[[1]][length(strsplit(outdir,"")[[1]])]!='/') outdir=paste(outdir, '/', sep='')
	}
	if(!file.exists(twoBitFile)){
		stop("Specified .2bit file doesnt exist or path is not correct")
	}

	if(chrm=='all'){		
		twobitinfo(twoBitFile, "twobit.temp")
		info=read.table("twobit.temp")	
		for(i in 1:dim(info)[1]){
			print(paste("printing", info[i,1]))
			.C("twoBitToFa",as.character(info[i,1]),as.integer(1),as.integer(info[i,2]),as.character(twoBitFile),as.character(paste(outdir, info[i,1], ".fa", sep="")),PACKAGE="zinba")
		}
		cat("\ntwobittofa complete\n")
		unlink("twobit.temp")
	}else{
		.C("twoBitToFa",as.character(chrm),as.integer(start),as.integer(end),as.character(twoBitFile),as.character(gcSeq),PACKAGE="zinba")
		cat("\ntwobittofa complete\n")
	}
}
	
