zinba=function(outfile=NULL,seq=NULL,align=NULL,input="none",twoBit=NULL,basecountfile=NULL,threshold=0.05,numProc=1,
refinepeaks=1,printFullOut=0,filetype="bowtie",extension=NULL,broad=FALSE, mode="peaks", interaction=TRUE, FDR=TRUE,
genome="hg"
){
	if(is.null(outfile)){stop("output prefix must be specified")}
	if(is.null(seq)){stop("path to mapped experimental reads must be specified")}
	if(!file.exists(seq)){stop("sequencing file does not exist at the specified path, check path for errors")}
	if(!file.exists(align)){stop("mappability directory does not exist at the specified path, check path for errors")}
	if(!file.exists(input) & input!="none"){stop("sequencing file does not exist at the specified path, check path for errors")}
	if(!file.exists(twoBit)){stop(".2bit file does not exist at the specified path, check path for errors")}
	if(is.null(basecountfile) & refinepeaks==1){stop("basecount file needs to be specified for peak refinement")}
	if(is.null(extension)){stop("sequencing read extension length (average length of fragment library) needs to be specified")}

	#base parameters
	buildwin=1
	winSize=250	
	offset=125
	cnvWinSize=100000
	cnvOffset=2500
	pquant=1
	initmethod="count"
	printFullOut=1
	diff=0
	pWinSize=200
	method="mixture"
	selectmodel=TRUE
	selectchr="chr22"
	selecttype="dirty"
	winGap=0
	#peakconfidence=1-threshold	  this is now done in run.zinba

	#options	
	if(genome=="hg" || genome=="mm"){
		if(broad==TRUE) winGap=5000
		if(input=="none" & broad==FALSE){ 
			selectcovs=c("gcPerc", "align_perc", "exp_cnvwin_log")
		}else if(input=="none" & broad==TRUE){
			selectcovs=c("gcPerc", "align_perc")
		}else if(input!="none"){
			selectcovs=c("input_count")
			selecttype="complete"
		}
		if(genome=="mm") selectchr="chr19"

	}else if(genome=="othersmall" || genome=="ce" ||genome=="dm"){
		winSize=200	
		offset=40
		if(broad==TRUE) winGap=2000
		if(input=="none" & broad==FALSE){ 
			selectcovs=c("gcPerc", "exp_cnvwin_log", "align_perc")
		}else if(input=="none" & broad==TRUE){
			selectcovs=c("gcPerc","align_perc")
		}else if(input!="none"){
			selectcovs=c("input_count")
			selecttype="complete"
		}
		if(genome=="dm")	selectchr="chr2L"
		else if(genome=="ce")	selectchr="chrI"
	}else{
		stop("'genome' should be either 'hg', 'mm', 'dm', 'ce', 'othersmall'")
	}

	


	if(mode=="peaks"){
		run.zinba(
			align=align,
			numProc=numProc,
			seq=seq,
			input=input,
			basecountfile=basecountfile,
			filetype=filetype,
			offset=offset,
			buildwin=1,
			outfile=outfile,
			threshold=threshold,
			twoBit=twoBit,
			cnvOffset=cnvOffset,
			pquant=1,
			winGap=winGap,
			cnvWinSize=cnvWinSize,	
			initmethod="count",
			printFullOut=printFullOut,
			winSize=winSize,
			diff=0,
			pWinSize=200,	
			extension=extension,
			method="mixture",
			refinepeaks=refinepeaks,
			selectmodel=TRUE,
			selectchr=selectchr,
			selecttype=selecttype,
			selectcovs=selectcovs,
			FDR=FDR,
			interaction=interaction #,
			#peakconfidence=peakconfidence  this is now deprecated
		)
	}else if(mode=="CNV"){
		run.zinba(
			align=align,
			numProc=numProc,
			seq=seq,
			input=input,
			basecountfile=basecountfile,
			filetype=filetype,
			offset=2500,
			buildwin=1,
			outfile=outfile,
			threshold=0.05,
			twoBit=twoBit,
			cnvOffset=2500,
			pquant=1,
			winGap=10000,
			cnvWinSize=100000,	
			initmethod="count",
			printFullOut=1,
			winSize=10000,
			diff=0,
			pWinSize=200,	
			extension=extension,
			method="mixture",
			refinepeaks=refinepeaks,
			selectmodel=FALSE,
			formula=exp_count~1,
			formulaE=exp_count~1,
			formulaZ=exp_count~1,
			FDR=FALSE,
			#peakconfidence=0.95,
			interaction=interaction
		)
	}
} 
