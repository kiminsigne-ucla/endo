basecountimport=function(inputfile,winlist,threshold=.01,method='pscl',printFullOut=0,outputfile,twobitfile,chromosome='all', winGap=0, FDR=FALSE){
		library(R.utils)
    if(!file.exists(inputfile)){
    	stop(paste("Input file,",inputfile,", not found.  Check to see if you mean to specify the compressed version of bascountfile instead ( filepath + .bin.gz)",sep=" "))
    }
    if(file.info(inputfile)$size==0){
    	stop(paste("Input file ",inputfile, "has size 0, chegeck for disk space issues",sep=" "))
    }
    if(!file.exists(twobitfile)){
        stop(paste("twoBit file not found,",twobitfile,sep=" "))
    }
    if(!file.exists(winlist)){
        stop(paste("Window list file,",winlist,", not found.  You may have run out of disk space or critical errors have occurred during mixture regression analysis.",sep=" "))
    }else if(file.info(winlist)$size == 0){
        stop(paste("Window list file,",winlist,", has 0 size, you may have run out of disk space or error occurred during mixture regression",sep=" "))			
		}else if(countLines(winlist)==0){
        stop(paste("0 lines found in window list file,",winlist,", you may have run out of disk space or error occurred during mixture regression",sep=" "))			
		}else{
			#check to see if files listied in winlist file exist
			files=as.matrix(read.table(winlist))
			exist=unlist(lapply(files,file.exists))
			if(any(exist==FALSE)){
				cat("Warning: Some of the files listed in the winlist file do not exist, deleting these files from winlist\n")
				if(length(which(exist==F)) == countLines(winlist)){
					stop(paste("None of the .wins files in ",winlist,
							" exist, either you have run out of disk space or serious error has occurred during model selection"))
				}else{
					cat("The following files do not exist:\n")
					print(files[which(exist==F)])
					cat("Removing files from winlist:\n")
					write.table(files[exist==T], winlist, col.names=F, row.names=F, quote=F)
				}

			}
			
			#only existing files are left, or if all did not exist, program has stopped here
			haslines=unlist(lapply(files,countLines))			
			hassize = unlist(lapply(files,function(x) file.info(x)$size))	

			if(any(hassize==0)){
					cat("Warning: Some of the files listed in the winlist file have 0 size, deleting these files from winlist\n")
				if(length(which(hassize==0)) == countLines(winlist)){
					stop(paste("None of the .wins files in ",winlist,
							" have non-zero size, either you have run out of disk space or serious error has occurred during model selection"))
				}else{
					cat("The following files have 0 size:\n")
					print(files[which(hassize==0)])
					cat("Removing files from winlist:\n")
					write.table(files[hassize==0], winlist, col.names=F, row.names=F, quote=F)
				}
			}else if(any(haslines==0)){
					cat("Warning: Some of the files listed in the winlist file have 0 lines, deleting these files from winlist\n")
				if(length(which(haslines==0)) == countLines(winlist)){
					stop(paste("None of the .wins files in ",winlist,
							" have more than 0 lines, either you have run out of disk space or serious error has occurred during model selection"))
				}else{
					cat("The following files have 0 size:\n")
					print(files[which(haslines==0)])
					cat("Removing files from winlist:\n")
					write.table(files[haslines==0], winlist, col.names=F, row.names=F, quote=F)
				}
			}
		}

    if(is.null(outputfile)){
    	stop(paste("Need to specify an outputfile,",outputfile,sep=" "))
    }

		if(length(grep(".bin.gz", inputfile))>0){
			cat("Using compressed binary version of basecount file\n")
			cat("Uncompressing file", inputfile,"\n")
			library(R.utils)
			files=untar(tarfile=inputfile, exdir=getParent(inputfile), list=T)
			write.table(files, getAbsolutePath(paste(getParent(inputfile),files, sep="")),col.names=F, row.names=F, quote=F)
			untar(tarfile=inputfile, exdir=getParent(inputfile), compressed="gzip")
			binary=1
			inputfile2=gsub("\\.gz$", "", inputfile)
		}else{
			cat("Assuming uncompressed version of basecount file (file in path does not have .bin.gz extension used for compressed files\n")
			binary=0
			inputfile2=inputfile
		}

    cReturn <- .C("getSeqCountProfile",as.character(inputfile2),as.character(winlist),
				as.double(threshold),as.character(method),as.integer(printFullOut),
				as.character(outputfile),as.character(twobitfile),as.character(chromosome), 
				as.integer(winGap), as.integer(FDR^2),as.integer(binary),PACKAGE="zinba")
	
		if(binary==1) unlink(gsub("\\.gz$", "", inputfile))
}
