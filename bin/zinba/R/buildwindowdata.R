buildwindowdata=function(seq,input="none",align,twoBit,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,filelist,filetype="bowtie",extension, outdir="default", numProc=1){
	if(!file.exists(align)){
		stop(paste("Specified Alignability directory not found,",align,sep=" "))
	}else{
		#if alignability path does not end in /, then put it in
		if( strsplit(align,"")[[1]][length(strsplit(align,"")[[1]])]!='/') align=paste(align, '/', sep='')
		if(!isDirectory(align)) stop("specified alignability directory path is not a directory")
		if(length(dir(align))==0) stop("There is nothing in your specified alignability directory, run generateAlignability again")
	}	
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
	
	require(R.utils)
	#require(multicore)
	#check if gzipped binary wig files
	#bwiggzfiles = file.access(dir(align,pattern="\\.bwig.gz$",full.names=T))
	#bwigfiles = file.access(dir(align,pattern="\\.bwig$",full.names=T))
	wigfiles = file.access(dir(align,pattern="\\.wig$",full.names=T))
	allfiles = dir(align,pattern="\\.*wig*",full.names=T)

	if(length(allfiles) != sum(file.access(allfiles)==0)){
			cat("The following alignability files could not be accessed\n")
			print(allfiles[file.access(allfiles)!=0])
			stop("It is likely you do not have permission to access these files\n")
	} 

	#if( (length(bwiggzfiles) > 0 | length(bwigfiles) > 0) & length(wigfiles) == 0 ){
		#all though either be bwig.gz or bwig, otherwise .wig
		#here, unzip those that are .gz
	#	cat("Found compressed binary alignability files in", align, ", uncompressing\n")
	#	if(length(bwiggzfiles) > 0) mc=mclapply(names(bwiggzfiles)[bwiggzfiles == 0], function(x) gunzip(x, remove=F, overwrite=T), mc.cores = numProc)
	#	cat("Uncompressiing complete\n")
		#check for disk space issues
	#  bwigfiles=file.access(dir(align,pattern="\\.bwig$",full.names=T))
	#	bwig_sizes=file.info(names(bwigfiles))$size
	#	if(any(bwig_sizes==0)){
	#		cat("The following alignability files have a size of zero\n")
	#		print(names(bwigfiles)[which(bwig_sizes==0)])
	#		stop("It is likely you have run of out disk space\n")
	#	}

		#check for unzipping issues
	#	if(length(bwigfiles) != length(bwiggzfiles)){
	#		cat("Not all files could be unzipped, listed below:\n")
	#		#print(setdiff(allfiles, names(bwigfiles)))
	#		stop("Check to see whether you are used all of your disk space")
	#	}
		
		#set buildwindows binary flag
	#	binary=1
	#}else{
		#then all we have are the regular uncompressed wig files
		#check for disk space issues
		wig_sizes=file.info(names(wigfiles))$size
		if(any(wig_sizes==0)){
			cat("The following alignability files have a size of zero\n")
			print(names(bwigfiles)[which(bwig_sizes==0)])
			stop("It is likely you have run of out disk space\nPlease free some disk space and regenerate your alignability directory")
		}

		#set buildwindows binary flag
		binary=0
		
	cReturn <- .C("buildWindows",as.character(seq),as.character(input),as.character(align),as.character(twoBit),as.integer(winSize),as.integer(offset),as.integer(cnvWinSize),as.integer(cnvOffset),as.character(filetype),as.character(filelist),as.integer(extension), as.character(outdir), as.integer(binary),PACKAGE="zinba")

	#remove all bwig files if they exist, and no wigfiles exist
	#if( length(bwigfiles) > 0 & length(wigfiles) == 0 ){
		#If there are any bwig files and no wig files, then all bwig files should be gzipped back
	#	cat("Removing uncompressed binary alignability files in", align, "\n\n")
	#	mc=mclapply(names(bwigfiles)[bwigfiles == 0], unlink, mc.cores = numProc)
	#	bwiggzfiles = file.access(dir(align,pattern="\\.bwig.gz$",full.names=T))
		
	#}

	gc()
}
