peakbound=function(bpprofile,output,pwinSize=200, winSize,quantile=.75,minscore=0, extension=200){
	    if(!file.exists(bpprofile)){stop("Error, file pertaining to SBP for merged regions does not exist.  Check .winlist file for errors.")}
	    localMaximum=function(x) {
		    #from MassSpecWavelet Package in BioConductor, written by Pan Du and Simon Lin, modified by Naim Rashid	
		    len <- length(x)
		    rNum <- ceiling(len/pwinSize)
		    y <- matrix(c(x, rep(x[len], rNum * pwinSize - len)), nrow = pwinSize)
		    y.maxInd <- apply(y, 2, which.max)
		    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > 
		        x[pwinSize]))
		    localMax <- rep(0, len)
		    localMax[(selInd - 1) * pwinSize + y.maxInd[selInd]] <- 1
		    shift <- floor(pwinSize/2)
		    rNum <- ceiling((len + shift)/pwinSize)
		    y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * pwinSize - 
		        len - shift)), nrow = pwinSize)
		    y.maxInd <- apply(y, 2, which.max)
		    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > 
		        x[pwinSize]))
		    localMax[(selInd - 1) * pwinSize + y.maxInd[selInd] - shift] <- 1
		    maxInd <- which(localMax > 0)
		    selInd <- which(diff(maxInd) < pwinSize)
	    	if (length(selInd) > 0) {
		        selMaxInd1 <- maxInd[selInd]
		        selMaxInd2 <- maxInd[selInd + 1]
		        temp <- x[selMaxInd1] - x[selMaxInd2]
	        	localMax[selMaxInd1[temp <= 0]] <- 0
		        localMax[selMaxInd2[temp > 0]] <- 0
		    }
	    	#delete maxes lower than prob threshold and in flanking
	    	localMax[which(x<quantile(x,quantile))]=0
		localMax[which(x<minscore)]=0
		localMax[c(1:250,(len-250+1):len)]=0
	    	#ensure global max is always tagged
	    	localMax[which.max(x[251:(len-250)])+250]=1
		med=median(x[251:(len-250)])
		return(c(length(which(localMax>0)),med,which(localMax>0)))
	} 
      peakbound2 <- function(f, bpprofile, output) {
       f.check <- function(x) {
         x <- f(x)
         if(!is.numeric(x)) stop("Need a numeric result")
         as.integer(x)
       }
       .Call("peakboundc", body(f.check), as.character(bpprofile), as.character(output),new.env(), PACKAGE='zinba')
      }
      peakbound2(localMaximum, bpprofile, output)
}
