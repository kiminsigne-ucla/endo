run.zinba=function(filelist=NULL,formula=NULL,formulaE=NULL,formulaZ=NULL,
		outfile=NULL,seq=NULL,align=NULL,input="none",twoBit=NULL,
		winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,
		basecountfile=NULL,threshold=0.01,peakconfidence=.8,tol=10^-5,
		numProc=1,buildwin=1, winGap=0,pWinSize=200,pquant=1,refinepeaks=1,
		printFullOut=0,method="mixture",initmethod="count",diff=0,
		filetype="bowtie",extension, cleanup=FALSE, selectmodel=FALSE, 
		selectchr=NULL, selecttype="dirty", selectcovs=NULL, FDR=FALSE, 
		interaction=TRUE, model = "FMR"
		){

	#since peakconfidence is now deprecated, set as 1-threshold	
	peakconfidence=1-threshold	

	parameters=list(filelist=filelist,formula=formula,formulaE=formulaE,
								formulaZ=formulaZ,outfile=outfile,seq=seq,align=align,
								input=input,twoBit=twoBit,winSize=winSize,offset=offset,
								cnvWinSize=cnvWinSize,cnvOffset=cnvOffset,basecountfile=basecountfile,
								threshold=threshold,peakconfidence=peakconfidence,tol=tol,
                numProc=numProc,buildwin=buildwin, winGap=winGap,pWinSize=pWinSize,
								pquant=pquant,refinepeaks=refinepeaks,printFullOut=printFullOut,
								method=method,initmethod=initmethod,diff=diff,filetype=filetype,
								extension=extension, cleanup=cleanup, selectmodel=selectmodel, 
								selectchr=selectchr, selecttype=selecttype, selectcovs=selectcovs, FDR=FDR,
								interaction=interaction)

	print("parameters in effect:")
	print(parameters)

  rmc <- require(doParallel)
  rdmc <- require(doMC)
  rfor <- require(foreach)
  if(rmc == FALSE){
    stop(paste("doParallel package not available, required for ZINBA"))
  }
  
	if(rdmc == FALSE){
	  stop(paste("doMC package not available, required for ZINBA"))
	}
  
	if(rfor == FALSE){
	  stop(paste("foreach package not available, required for ZINBA"))
	}
  
	require(R.utils)
	require(quantreg)
	
	time.start <- Sys.time()
	
	
        #####################################################################################################
	if(is.null(outfile)) stop("output prefix must be specified")
	if(is.null(align)) stop("alignability directory must be specified")
	if(is.null(seq)) stop("path to experimental mapped reads 'seq' must be specified")
	if(is.null(twoBit)) stop("path to .2bit file 'twoBit' must be specified")
	if(is.null(basecountfile) & refinepeaks == 1) stop("path to read overlap file 'basecountfile' must be specified")


	#create subdirectory to hold intermediate files to be used later
	outfile_subdir=paste(outfile,"_files/", sep="")
	if(!dir.create(outfile_subdir, showWarnings=FALSE )& buildwin==1){
		#if cannot create directory, delete existing directory and try to create again
		cat(paste("\nOverwriting previously existing directory ",outfile_subdir, "\n",sep=""))
		unlink(outfile_subdir, recursive=T)
		if(!dir.create(outfile_subdir,showWarnings=F)){
			#if fails again, print error and set output directory to same as main files at "outpath"
			cat("\nCould not create subdirectory to hold intermediate files, placing intermediate files with main output files\n")
			outfile_subpath=outfile
		}else{
			#if retry successful, set path prefix to where built 
			#window analysis files will be place in subdirectory			
			slashindex=which(substring(outfile,1:nchar(outfile),1:nchar(outfile))=="/")
			if(length(slashindex>0)){
				outfile_subpath=paste(outfile,"_files","/"
					,substr(outfile,slashindex[length(slashindex)]+1,nchar(outfile)),sep="")
			}else{
				outfile_subpath=paste(outfile,"_files","/",outfile,sep="")	
			}
		}
	}else{
		#if initial directoy creation sucessful, set path prefix to where 
		#built window analysis files will be place in subdirectory			
		slashindex=which(substring(outfile,1:nchar(outfile),1:nchar(outfile))=="/")
		if(length(slashindex>0)){
			outfile_subpath=paste(outfile,"_files","/",substr(outfile,slashindex[length(slashindex)]+1,
				nchar(outfile)),sep="")
		}else{
			outfile_subpath=paste(outfile,"_files","/",outfile,sep="")	
		}
	}
        
	#####################################################################################################
	#buildwindows
	if(buildwin==1){
	    if(is.null(filelist)) filelist=paste(outfile_subpath,".list",sep="")	
	    cat(paste("\n--------BEGIN BUILDING WINDOW DATA--------",as.character(Sys.time()),"\n"))
      buildwindowdata(seq=seq,align=align,input=input,twoBit=twoBit,winSize=winSize,offset=offset,
				cnvWinSize=cnvWinSize,cnvOffset=cnvOffset,filelist=filelist,filetype=filetype,
				extension=extension, outdir=outfile_subdir, numProc=numProc)
	}

	if(refinepeaks==1 && is.null(basecountfile)){
		stop(paste("Basecount file must be specified if refinepeaks = 1, currently",basecountfile,sep=" "))
	}else if (is.null(filelist)){
		stop(paste("Need list of files ",filelist,sep=" "))
	}else if (!file.exists(filelist)){
		stop(paste("File list ",filelist,
		" does not exist, check whether your reads have been formatted properly" ,sep=" "))
	}else if (file.info(filelist)$size == 0){
		stop(paste("File list ",filelist,
		" has 0 size, check whether you have run out of disk space" ,sep=" "))
	}else if (countLines(filelist)==0 ){
		stop(paste("File list ",filelist,
		" has 0 lines, check whether you have run out of disk space or if your reads are formatted properly" ,sep=" "))
	}else if(method != 'pscl' && method != 'mixture'){
		stop(paste("Method should be either pscl or mixture, currently",method))
	}else{
		#set prefixes of outputfiles
		#Before reading file list, check if all built files exist		
		checkfiles=as.character(unlist(read.table(filelist, sep=";")))
		file.info2=function(x){ return(file.info(x)$size) }
			
		if(any(unlist(lapply(checkfiles, file.exists))==FALSE)){
			cat("The following files from buildwindows do no exist\n")
			print(checkfiles[which(unlist(lapply(checkfiles, file.exists))==FALSE)])
			stop("Some files did not print, check the format of your reads, whether you have run out of disk space, or if you have permission to write to the _files/ directory")
		}else if(any(unlist(lapply(checkfiles, file.info2))==0)){
			#if all exist, check if all built files have non-zero size (no disk space issues)
			cat("The following files from buildwindows have 0 size\n")
			print(checkfiles[which(unlist(lapply(checkfiles, file.info2))==0)])
			stop("Check whether you have run out of disk space, or if you have permission to write to the _files/ directory")
		}else if(any(unlist(lapply(checkfiles, countLines))==0)){
			#if all exist, check if all built files have non-zero size (no disk space issues)
			cat("The following files from buildwindows have 0 lines\n")
			print(checkfiles[which(unlist(lapply(checkfiles, countLines))==0)])
			stop("Check whether you have run out of disk space, or if you have permission to write to the _files/ directory")
		}

		params=scan(filelist,what=character(0),quiet=T)
		winlist=paste(outfile_subpath,".winlist",sep="")
		peakout=paste(outfile,".peaks",sep="")
		broadpeakout=paste(outfile,".peaks.broad.bed",sep="")
		bpout=paste(outfile_subpath,".bpcount",sep="")

		if(selectmodel==FALSE){
			if(is.null(formulaE)){
				formulaE=exp_count~1
			}
			if(is.null(formulaZ)){
				formulaZ=formula
			}
			if(is.null(formula)){
				stop("No background formula specified, so must perform model selection procedure")
			}
			if(!inherits(formula, "formula")){
				stop("Check your background component formula, not entered as a formula object")
			}else if(!inherits(formulaE, "formula")){
				stop("Check your enrichment component formula, not entered as a formula object")
			}else if(!inherits(formulaZ, "formula")){ 
				stop("Check your zero-inflated component formula, not entered as a formula object")
			}
	
		}else if(selectmodel==TRUE){
			#start optional model selection
			if(is.null(selectchr)){
				stop("need to specify which chromosome to apply model selection to")
			}else if(length(selectchr)>1){
				stop("only one chromsome can be specified for selectchr, typically a smaller chromsome in larger genomes")
			}

			if(is.null(selectcovs)) stop("need to specify which covariates to use in model selection")
			supported=c("gcPerc", "align_perc", "input_count", "exp_cnvwin_log")
			if( sum(selectcovs %in% supported)!=length(selectcovs)){
				stop(paste("Covariate ",selectcovs[which(selectcovs %in% supported == FALSE)], 
				" not found in list of supported covariates"))
			}
			cat(paste("--------STARTING MODEL SELECTION--------",as.character(Sys.time()),"\n\n")) 
			data=NULL
			splitlist=strsplit(unlist(strsplit(params,";")),"_")
			for(i in 1:length(splitlist)){
				if(any(splitlist[[i]]==selectchr)){
					data=unlist(strsplit(params,";"))[i]
					break
				}
			}	
			if(is.null(data)){
				data=unlist(strsplit(params,";"))[i]
				cat(paste("\nSpecified chromosome not found or specified, using", data,"\n"))
			}
			model=covariateselect(file=data, selection=selecttype,
					loc=paste(outfile_subpath,".model",sep=""),
					covs=selectcovs, numProc=numProc, interaction=interaction)
			formula=model[[1]]
			formulaE=model[[2]]
			formulaZ=model[[3]]
			cat("Background formula is:\n\t");print(formula);
			cat("Enrichment formula is:\n\t");print(formulaE);
			cat("Zero-inflated formula is:\n\t");print(formulaZ);
			cat(paste("--------MODEL SELECTION COMPLETE--------",as.character(Sys.time()),"\n\n")) 
		}	
		 
		#begin mixture regression (parallelized)
	  if(rmc == TRUE && rdmc == TRUE && rfor == TRUE){
			cat(paste("--------GETTING ENRICHED WINDOWS--------",as.character(Sys.time()),"\n\n")) 		
	    #registerDoMC(numProc)
	  	#mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
	    #getDoParWorkers()
			cl <- makeCluster(numProc)
			registerDoParallel(cl)
	    winfiles <- foreach(i=1:length(params),.combine='rbind',.inorder=FALSE,
				.errorhandling="remove") %dopar%
			{
				zinba::getsigwindows(file=params[i],formula=formula,formulaE=formulaE,
					formulaZ=formulaZ,threshold=threshold,winout=outfile_subpath,
					peakconfidence=peakconfidence,tol=tol,method=method,printFullOut=printFullOut,
					initmethod=initmethod, FDR=FDR,  model = model
				)
	    }
	    write.table(winfiles,winlist,quote=F,row.names=F,col.names=F)
		  cat(paste("--------WINDOW ANALYSIS COMPLETE--------",as.character(Sys.time()),"\n\n"))		
	  }else{
		  #if parallelization fails due to lack of packages, resort to non-parallelized version
		  cat(paste("--------GETTING ENRICHED WINDOWS--------",as.character(Sys.time()),"\n\n")) 	
	    winfiles = rep("", length(params))
	    for(i in 1:length(params)){
	       winfiles[i] <- zinba::getsigwindows(file=params[i],formula=formula,formulaE=formulaE,
	                               formulaZ=formulaZ,threshold=threshold,winout=outfile_subpath,
	                               peakconfidence=peakconfidence,tol=tol,method=method,printFullOut=printFullOut,
	                               initmethod=initmethod, FDR=FDR, model = model
	              ) 
				
	    }
			write.table(winfiles,winlist,quote=F,row.names=F,col.names=F)
		  cat(paste("--------WINDOW ANALYSIS COMPLETE--------",as.character(Sys.time()),"\n\n"))		
		}
	    
		if(refinepeaks==1){
			#merge windows and refine regions		
			cat(paste("--------MERGE WINDOWS AND REFINE PEAKS (no parallelization)--------",
						as.character(Sys.time()),"\n"))
			getrefinedpeaks(winlist=winlist,basecountfile=basecountfile,bpout=bpout,peakout=peakout,
					twoBit=twoBit,winSize=winSize,pWinSize=pWinSize,pquant=pquant,printFullOut=printFullOut,
					peakconfidence=peakconfidence,threshold=threshold,method=method, 
					winGap=winGap, extension=extension, FDR=FDR
			)
	  }else{
			#merge windows only
			if(method=="mixture") threshold = peakconfidence
			cat(paste("--------MERGE WINDOWS --------",as.character(Sys.time()),"\n"))
			collapsewindows(winlist=winlist,printFullOut=printFullOut,thresholds=threshold,method=method, 
				winGap=winGap, FDR=FDR, output=broadpeakout)
		}
		
	}
        #####################################################################################################
	cat(paste("\n--------ZINBA COMPLETE--------",as.character(Sys.time()),"\n\n"))
	
	if(cleanup==TRUE) unlink(outfile_subdir, recursive=T)
	time.end <- Sys.time()
	print(difftime(time.end,time.start))

}
