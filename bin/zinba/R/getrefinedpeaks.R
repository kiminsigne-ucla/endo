getrefinedpeaks=function(winlist,basecountfile,bpout,peakout,twoBit,winSize,pWinSize=200,printFullOut=0,pquant=0.75,threshold=.01,peakconfidence=0.8,method='pscl',minscore=0, winGap=0, extension=200, FDR=FALSE){
    if(!file.exists(winlist)){stop("Specified winlist doesnt exist or is not correct")}
    if(!file.exists(basecountfile)){stop("Specified basecountfile doesnt exist or is not correct")}
    if(!file.exists(twoBit)){stop("Specified twoBit file doesnt exist or is not correct")}
    if(method=='mixture' & FDR==FALSE){
        threshold = peakconfidence  #need to reconcile this with latest change
    }
    cat(paste("\nThreshold is ",threshold,"\n",sep=""))
    basecountimport(inputfile=basecountfile,winlist=winlist,threshold=threshold,
				method=method,printFullOut=printFullOut,outputfile=bpout,
				twobitfile=twoBit, winGap=winGap, FDR=FDR)
    peakbound(bpprofile=bpout,output=peakout,pwinSize=pWinSize,winSize=winSize,
				quantile=pquant,minscore=minscore, extension=extension)
    unlink(bpout)
}
