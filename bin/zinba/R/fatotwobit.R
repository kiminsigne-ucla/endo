fatotwobit=function(faFile=" ",fadir = " ",outFile=NULL){
    	if(fadir!=" " & !file.exists(fadir)){
		stop("Specified directory 'fadir' doesnt exist or is not correct")
	}else{
		if(strsplit(fadir,"")[[1]][length(strsplit(fadir,"")[[1]])]!='/' & fadir!=" ") fadir=paste(fadir, '/', sep='')
	}
   
	if(fadir!=" " & faFile==" "){
		fafiles=paste(paste(fadir, dir(fadir, pattern="\\.fa"), sep=""), collapse="@")
			cat("processing...\n")
			.C("faToTwoBit",as.character(fafiles),as.character(outFile),PACKAGE="zinba")		
		cat("Conversion of fasta files in specified directory to .2bit complete\n")  
	}else if(faFile!=" " & fadir==" "){
		cat("processing...\n")
		.C("faToTwoBit",as.character(faFile),as.character(outFile),PACKAGE="zinba")
		cat("Conversion of specified fasta file to .2bit complete\n")
	}else{
		stop("Either need to specify a directory containing .fa files to be converted or specify a specific faFile")
	}

	
}

