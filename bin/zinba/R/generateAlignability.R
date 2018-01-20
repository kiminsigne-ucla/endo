process=function(maptargz,mapdir, chr, outdir, twoBitFile, athresh, extension,maplist,i){
	untar(maptargz,compressed=T, file=chr[i], exdir=getParent(mapdir))
	maplist2=paste(maplist,"_", i, sep="")
	write.table(getAbsolutePath(chr[i]), maplist2, quote = F, row.names = F, col.names = F)
	b=.C("binAlignAdjust", as.character(maplist2), as.character(outdir),as.character( twoBitFile) , as.integer(athresh), as.integer(extension/2),package="zinba")
	unlink(maplist2)
	unlink(getAbsolutePath(chr[i]))
	adjmapchri=paste(outdir,gsub("b.out",".bwig",basename(chr[i])), sep="")
	gzip(adjmapchri, overwrite=T, remove=TRUE)
	cat("processed", getAbsolutePath(adjmapchri), "\n")
	return(getAbsolutePath(adjmapchri))
}


generateAlignability=function(mapdir, outdir="", athresh=1, extension=0, 
	twoBitFile, binary=0, cleanup=T,savespace=T,numProc=1,maptargz=""){
	library(R.utils)
	if(file.exists(outdir)){
		outdir=paste(getAbsolutePath(outdir),"/",sep="")
		cat("\n\nOverwriting existing directory at", outdir,"\n")
		unlink(outdir, recursive=T) 
		suppressWarnings(dir.create(outdir))
	}else{
		outdir=paste(getAbsolutePath(outdir),"/",sep="")
		cat("Creating mappability directory at", outdir,"\n")
		dir.create(outdir)
	} 	 	

	if(!file.exists(twoBitFile)){
		stop("Specified twoBit file doesnt exist or is not correct")
	}else{
		twoBitFile=getAbsolutePath(twoBitFile)
	} 
	
  cat("\nStarting generateAlignability\n")
	if(binary==0){
		if(!file.exists(mapdir))
			{stop("Specified directory 'mapdir' doesnt exist or is not correct")
		}else{
			#if alignability path does not end in /, then put it in
			if( strsplit(mapdir,"")[[1]][length(strsplit(mapdir,"")[[1]])]!='/') mapdir=paste(mapdir, '/', sep='')
		} 
		currdir=getwd()
		setwd(mapdir)
		write.table(dir("." ,pattern="\\.out"), "map.list", quote=F, row.names=F, col.names=F)
		convertmappability(inputfile='map.list', outputfile='temp.wig')
		setwd(currdir)
		alignAdjust(inputfile=paste(mapdir,'temp.wig',sep=''), 		
			outDir=outdir,twoBitFile=twoBitFile,athresh=athresh,adjustsize=round(extension/2))
		unlink('temp.wig')
	}else{
		cat("Binary mode: Reading file names within",maptargz,"\n")
		chrtar=untar(maptargz, list=T)
		chr = chrtar[grep("chr", chrtar)]
		mapdir=getAbsolutePath(chrtar[1])
		maplist=paste(mapdir, "/map.list", sep="")
		currdir=getwd()

		if(savespace==F){
			if(file.exists(mapdir)) cat("Overwriting existing mappability directory",mapdir,"\n")
			unlink(mapdir, recursive=T)
			cat("Unpacking compressed mappability file",maptargz,"\n")
			untar(maptargz,compressed=T)
			write.table(dir(mapdir, full.names=T, pattern = "\\.out"), maplist, quote = F, row.names = F, 				col.names = F)
			cat("Creating adjusted mappability files\n")
			b=.C("binAlignAdjust", as.character(maplist), as.character(outdir),as.character
			( twoBitFile) , as.integer(athresh), as.integer(extension/2),package="zinba")
			adjmapchr=paste(outdir,dir(outdir, pattern = "\\.bwig"), sep="")
			cat("Compressing adjusted mappability files\n")
			for(i in 1:length(adjmapchr)){
				gzip(adjmapchr[i], overwrite=T, remove=TRUE)
			}
		}else{
			suppressPackageStartupMessages(require(multicore))
			suppressPackageStartupMessages(require(doMC))
			suppressPackageStartupMessages(require(foreach)) 
			registerDoMC(numProc)
			getDoParWorkers()
			if(file.exists(mapdir)) cat("\n\nOverwriting existing mappability directory",mapdir,"\n")
			unlink(mapdir, recursive=T)
			dir.create(mapdir)
			cat("Starting processing\n")
			mapfiles <- foreach(i=1:length(chr),.combine='rbind',.inorder=FALSE,.errorhandling="remove" ) %dopar%{
				process(maptargz,mapdir, chr, outdir, twoBitFile, athresh, extension,maplist, i)
			}
			cat("\nFinished Processing and compressing files\n")
			print(dir(outdir))
		}

		if(cleanup==T){ 
			unlink(mapdir, recursive=T)
		}
	cat("\n---------BINARY  ALIGN ADJUST COMPLETED SUCCESSFULLY ----------------")
	}

}
