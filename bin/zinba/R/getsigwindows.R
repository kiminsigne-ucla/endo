logsumexp=function(v){
	if(any(is.infinite(v))){
		stop("infinite value in v\n")
	}
	if(length(v)==1){ return(v[1]) }
	sv  = sort(v, decreasing=TRUE)
	res = sum(exp(sv[-1] - sv[1]))
	lse = sv[1] + log(1+res)
	lse
}

loglikfun=function(parms){
	mu1=parms$start$count1
	mu2=parms$start$count2
	prob0 =parms$start$zero
	prop1=parms$prop1
	prop2=parms$prop2
	theta1=parms$start$theta1
	theta2=parms$start$theta2
	Y=parms$Y
	loglik0=log(prob0+(1-prob0)*prop1*dnbinom(0, size = theta1, mu = mu1)+(1-prob0)*prop2*dnbinom(0, size = theta2, mu = mu2))
	loglik1=log((1-prob0)*prop1*dnbinom(Y, size = theta1, mu = mu1)+(1-prob0)*prop2*dnbinom(Y, size = theta2, mu = mu2))
	NAs=which(loglik1==-Inf)
	if(length(NAs>0)){
		loglik1[NAs]=apply(cbind(log((1-prob0[NAs])*prop1)+dnbinom(Y[NAs], size = theta1, mu = mu1[NAs],log=TRUE),log((1-prob0[NAs])*prop2)+dnbinom(Y[NAs], size = theta2, mu = mu2[NAs], log=TRUE)), 1, logsumexp)
	}
	loglik=sum(loglik0*(Y <= 0) + loglik1*(Y > 0))
	loglik
}			

getsigwindows=function(file,formula,formulaE,formulaZ,winout,
	threshold=.01,peakconfidence=.8, printFullOut=0,tol=10^-5,
	method="mixture",initmethod="count", diff=0,modelselect=FALSE, trace=0, 
	FDR=FALSE)
{

	time.start <- Sys.time()
	library(MASS)
	options(scipen=999)

	if(!inherits(formula, "formula"))  
		stop("Check your background component formula, not entered as a formula object")
	if(!inherits(formulaE, "formula")) 
		stop("Check your enrichment component formula, not entered as a formula object")
	if(!inherits(formulaZ, "formula")) 
		stop("Check your zero-inflated component formula, not entered as a formula object")
	
	winfile=NULL
	printflag = 0
	if(method=='pscl'){
	#not actively maintained
	  suppressPackageStartupMessages(library(quantreg))
		suppressPackageStartupMessages(require(qvalue))
		files = unlist(strsplit(file,";"))
			for(i in 1:length(files)){
				data=read.table(files[i], header=TRUE)
				chrm = data$chromosome[1]
				mf <- model.frame(formula=formula, data=data)
				X <- model.matrix(attr(mf, "terms"), data=mf)
				if(i == 1){
					a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE)
				}else{
						a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE,start=param)
				}
				q25=quantile(data$exp_count, 0.25)
				leverage=hat(X, intercept=FALSE)
				fdrlevel=threshold
				standardized=residuals.zeroinfl(a)/sqrt(1-leverage)
				pval=1-pnorm(as.matrix(standardized))
				fdr=qvalue(pval)
				numpeaks=length(which(fdr[[3]]<fdrlevel))
			if(printFullOut == 1){
				data=cbind(data, ((data$exp_count>q25)^2), 
					formatC(fdr[[3]],format="f",digits=16), standardized)
				colnames(data)[c(dim(data)[2]-2, dim(data)[2]-1, 
					dim(data)[2])]=c('q25','qvalue', 'residual')
			}else{
				data=cbind(data[1:3],((data$exp_count>q25)^2), 
					formatC(fdr[[3]],format="f",digits=16), standardized)
				colnames(data)=c('chromosome','start','stop','q25',
					'qvalue', 'residual')
			}
			param=list(count=a$coefficients$count, zero=a$coefficients$zero, theta=a$theta)
	
		### PRINT SIGNIFICANT WINDOWS
			print(paste("\nFor ",files[i],", found ",as.character(numpeaks),
					" significant wins",sep=''))
			winfile = paste(winout,"_",chrm,".wins",sep="")
			if(printflag==1){
				write.table(data,winfile,quote=F,sep="\t",row.names=F,col.names=F,
					append=T)
			}else{
				write.table(data,winfile,quote=F,sep="\t",row.names=F)
				printflag=1
			}
		}
		time.end <- Sys.time()
		print(difftime(time.end,time.start))
		return(winfile)
	
	}else if(method=='mixture'){
		files = unlist(strsplit(file,";"))
		fnum=1
		retry=0	#flag for retrial of offset
		while(fnum <= length(files)){
			fail=0 #failure flag for GLM, reset to zero for each offset
			#print(files[fnum])
			data=read.table(files[fnum], header=TRUE)
			chrm = data$chromosome[1]
			q25=quantile(data$exp_count, 0) #not used, set to 0
			mf <- model.frame(formula=formula, data=data)
			mfE <- model.frame(formula=formulaE, data=data)
			mfZ <- model.frame(formula=formulaZ, data=data)
			X <- model.matrix(attr(mf, "terms"), data=mf)
			XE <- model.matrix(attr(mfE, "terms"), data=mfE)
			XZ <- model.matrix(attr(mfZ, "terms"), data=mfZ)
			
			Y <- model.response(mf)
			n <- length(Y)
			kx <- NCOL(X)
	
			linkstr <- 'logit'
			linkobj <- make.link(linkstr)
			linkinv <- linkobj$linkinv
	
			#starting params for ze0 component
		   model_zero <-.C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), 
				M=as.integer(ncol(as.data.frame(XZ[,-c(1)]))), as.double(Y <= 0), as.double(rep(1,n)), as.double(rep(0,n)), 
				as.double(unlist(as.data.frame(XZ[,-c(1)]))),  as.integer(rep(1,n)),init=as.integer(1), 
				rank=integer(1), double(n*ncol(as.data.frame(XZ[,-c(1)]))), 
				fitted=as.double((rep(1,n) * (Y <= 0) + 0.5)/(rep(1,n) + 1)), double(n), 
				double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')
	
			#INITIALIZATION OF THE MIXTURE MODEL
			if(initmethod=='quantile'){
			  require(quantreg)
			 	if(fnum == 1){
					startprop=startenrichment(c(.15, .001), data, formula, formulaE, formulaZ, initmethod)
				}
				data2=data
				if(sum(colnames(data)=='input_count')==1){data2$input_count=exp(data2$input_count)-1}
				if(sum(colnames(data)=='exp_cnvwin_log')==1){data2$exp_cnvwin_log=exp(data2$exp_cnvwin_log)-1}
				prop2=max(c(startprop, 1000/n))
				prop1=1-prop2
				t=rq(formula, tau=1-prop2, data=data2, method='pfn')
				priorCOUNTweight=rep(10^-10, length(Y))	  
				priorCOUNTweight[as.double(which(t$residuals>quantile(t$residuals,1-prop2)))]=1-10^-10
				rm(data2)
  	  }else if(initmethod=='count'){
				if(fnum == 1){
					startprop=startenrichment(c(.15, .001), data, formula, formulaE,formulaZ,initmethod)
				}
				prop2=max(c(startprop, 1000/n))
				prop1=1-prop2
				n1  = round(length(Y) * (1 - prop2))
				priorCOUNTweight=rep(1-10^-10, length(Y))
				odY = order(Y)
				priorCOUNTweight[odY[1:n1]]=10^-10
			}else if(initmethod=='pscl'){
			  require(quantreg)
				mf2 <- model.frame(formula=formula, data=data)
				X2 <- model.matrix(attr(mf2, "terms"), data=mf2)
				if(fnum == 1){
					a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE)
				}else{
					a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE,start=param)
				}
				leverage=hat(X2, intercept=FALSE)
				fdrlevel=threshold
				standardized=residuals(a)/sqrt(1-leverage)
				pval=1-pnorm(as.matrix(standardized))
				fdr=qvalue(pval)
				peaks=which(fdr[[3]]<fdrlevel)
				prop2=length(peaks)/length(Y)
				prop1=1-prop2
				param=list(count=a$coefficients$count, zero=a$coefficients$zero, theta=a$theta)
				priorCOUNTweight=rep(10^-10, length(Y))	  
				priorCOUNTweight[peaks]=1-10^-10
				rm(X2)
				rm(mf2)
				rm(standardized)
				rm(leverage)
				rm(pval)
				rm(fdr)
				rm(a)
				gc()
	  	}
			

			model_count1 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), 
				M=as.integer(ncol(as.data.frame(X[,-c(1)]))), as.double(Y), as.double(1-priorCOUNTweight), 
				as.double(rep(0,length(Y))), as.double(unlist(as.data.frame(X[,-c(1)]))),  
				as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), 
				double(length(Y)*ncol(as.data.frame(X[,-c(1)]))), fitted=as.double(Y+(Y==0)/6), 
				double(length(Y)), double(length(Y)),scale=double(1), df_resid=integer(1), 
				theta=as.double(-1), package='zinba')
			model_count2 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), 
				M=as.integer(ncol(as.data.frame(XE[,-c(1)]))), as.double(Y), as.double(priorCOUNTweight), 
				as.double(rep(0,length(Y))), as.double(unlist(as.data.frame(XE[,-c(1)]))),  
				as.integer(rep(1,length(Y))),init=as.integer(1), 
				rank=integer(1),double(length(Y)*ncol(as.data.frame(XE[,-c(1)]))), fitted=as.double(Y+(Y==0)/6), 
				double(length(Y)), double(length(Y)),scale=double(1), df_resid=integer(1), 
				theta=as.double(-1), package='zinba')
				
			#starting prior probs
			mui1  <- model_count1$fitted
			mui2  <- model_count2$fitted
			probi0 <- model_zero$fitted
			#start mean vector
			start <- list(count1 = model_count1$fitted, count2 = model_count2$fitted,
				zero = model_zero$fitted, zerocoef=model_zero$coefficients)
			start$theta1 <- 1
			start$theta2 <- 1
						
			probi1  <- (1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)/
				(probi0*(Y <= 0)+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1) 
				+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi2  <- (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2)/
				(probi0*(Y <= 0)+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1) 
				+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi0=probi0/(probi0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)
				+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi0[ Y > 0]=0
				
			NAs=which(probi1=='NaN'| probi2=='NaN')		
			if(length(NAs>0)){
					probi1[NAs]=0
					probi2[NAs]=1
			}	
			ll_new <- loglikfun(list(start=start, prop1=prop1, prop2=prop2, Y=Y))
			ll_old <- 2 * ll_new	  
	
			if(!require("MASS")) {
				ll_old <- ll_new
				warning("EM estimation of starting values not available")
			}
			ll=matrix(0, 1, 10000)
			ll[1]=ll_new
			i=2
				
			while(abs((ll_old - ll_new)/ll_old) > tol) {
				ll_old <- ll[max(1, i-10)]
				prop1=sum(probi1)/(sum(probi1)+sum(probi2))
				prop2=sum(probi2)/(sum(probi1)+sum(probi2))
				if(trace==1){			 
					print(c(ll_new, prop2))
				}
	
				if(prop1<.5){
					if(modelselect==F){
						print(paste("The estimated proportion of enrichment for  ", 
							files[fnum], 
							" has exceeded 0.5, suggesting difficulty in estimating enrichment.  Switching to more conservative model"))
					}
					break
				}
				#updated values for parameters of component means
					 
				model_zero <- .C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), 
					M=as.integer(ncol(as.data.frame(XZ[,-c(1)]))),as.double(probi0), as.double(rep(1,n)), 
					as.double(rep(0,n)), as.double(unlist(as.data.frame(XZ[,-c(1)]))),  
					as.integer(rep(1,n)),init=as.integer(1), rank=integer(1), 
					double(n*ncol(as.data.frame(XZ[,-c(1)]))), fitted=as.double((rep(1,n)*probi0 + 0.5)/(rep(1,n) + 1)), 
					double(n), double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1),
					package='zinba')	
				model_count1 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), 
					M=as.integer(ncol(as.data.frame(X[,-c(1)]))), as.double(Y), as.double(probi1), 
					as.double(rep(0,length(Y))), as.double(unlist(as.data.frame(X[,-c(1)]))), 
					as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), 
					double(length(Y)*ncol(as.data.frame(X[,-c(1)]))), fitted=as.double(start$count1), double(length(Y)), 
					double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta1), 
					package='zinba')  
					model_count2 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), 
					M=as.integer(ncol(as.data.frame(XE[,-c(1)]))), as.double(Y), as.double(probi2), 
					as.double(rep(0,length(Y))), as.double(unlist(as.data.frame(XE[,-c(1)]))),  
					as.integer(rep(1,length(Y))),init=as.integer(1), 
					rank=integer(1), double(length(Y)*ncol(as.data.frame(XE[,-c(1)]))), fitted=as.double(start$count2), 
					double(length(Y)), double(length(Y)),scale=double(1), df_resid=integer(1), 
					theta=as.double(start$theta2), package='zinba')  
	
				#check for invalid fitted values, continue to next file if foun
				start <- list(count1 = model_count1$fitted, count2 = model_count2$fitted,zero = model_zero$fitted, zerocoef=model_zero$coefficients)
				start$theta1 <- model_count1$theta
				start$theta2 <- model_count2$theta
	
				mui1  <- model_count1$fitted
				mui2  <- model_count2$fitted
				probi0 <- model_zero$fitted
	
				if(any(!is.finite(mui1)) | any(!is.finite(mui2)) | any(!is.finite(probi0))){
					fail=1
					print(paste("none finite values obtained, likely failure of GLMs in", files[fnum]))
					break
				}
	
			probi1  <- (1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)/
					(probi0*(Y <= 0)+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)
					+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
				probi2  <- (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2)/
					(probi0*(Y <= 0)+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)
					+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
				probi0=probi0/(probi0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)
					+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
				probi0[ Y > 0]=0
				#extremely large counts can lead to small denominators in calculation of probi1, probi2
				NAs=which(probi1=='NaN'| probi2=='NaN')		
				if(length(NAs>0)){
					probi1[NAs]=10^-10
					probi2[NAs]=1
				}
				probi1[probi1<10^-10]=10^-10
				probi2[probi2<10^-10]=10^-10

				ll_new <- loglikfun(list(start=start, prop1=prop1, prop2=prop2, Y=Y))
				ll[i]=ll_new
				i=i+1 
				if(i>300){break}
				cat(".")
			}
			
			if(modelselect==TRUE){
				#if model selection is occuring, return all objects needed
				#note that we need to state that failure has occurred when we update ZINBA for model selection
				return(list(start=start,prop1=prop1, prop2=prop2, probi1=probi1, 
					probi2=probi2, probi0=probi0,ll=ll_new, logdimdata=log(dim(data)[1]), 
					fail=(prop1<.5)^2))
			}else if(fail==1 &  modelselect==FALSE & fnum<length(files)){
				#if GLM function fails, skip to next offset, fail is reset at start
				fnum=fnum+1
				next
			}else if(prop1<.5 & modelselect==FALSE & fnum<length(files)){
			#retry offset, do not go to next offset until # tries exceeds 2	
				#model did not converge, switch to more convervative ICL model 
				#if model selection file is available, otherwise switch to more
				#conservative intercept mode.  If all fails, then return current
				#and state that model selection needs to be performed with this chrm
				model=paste(winout,".model",sep="")
				if(file.exists(model) & retry==0){
					print("using ICL existing model file")
					final=read.table(model, header=T, sep="\t")
					bestICL=which.min(final$ICL)
					formula=as.formula(paste("exp_count~",final$formula[bestICL]))
					formulaE=as.formula(paste("exp_count~",final$formulaE[bestICL]))
					formulaZ=as.formula(paste("exp_count~",final$formulaZ[bestICL]))
					#first pass is made at ICL if available		
					retry=1
					next
				}else if(retry<=1){	
					print("defaulting enrichment to intercept")
					formulaE=exp_count~1
					retry=retry+1	
					next
				}else{
					#if number of tries exceeded, go to next offset
					#since not on last offset
					fnum=fnum+1
					#reset retry		
					retry=0
					next
				}
			}else if(fnum==length(files) & printflag==0 & (fail==1 | prop1 <.5)){
				#no files have printed successfully,and the last offset 
				#also has failed, return error
				stop("no offsets converged successfully")
			}else{
			  #success, print peaks
			  p=1-probi2
			  p2=rep(0,length(p))
			  p2[order(p)]=cumsum(p[order(p)])/(1:length(p))
		   if(FDR){
					numpeaks=length(which(p2<threshold))
		   }else{
		 	   	numpeaks=length(which(probi2>peakconfidence))
		   }
		
		   if(printFullOut == 1){
				data=cbind(data,((data$exp_count>q25)^2),formatC(probi2,format="f",digits=16), formatC(p2,format="f",digits=16)) #all non-zero set to 1
				colnames(data)[c(dim(data)[2]-2 ,dim(data)[2]-1,dim(data)[2])]=c('q25','peakprob','qvalue')
		   }else{
		   	data=cbind(data[1:3],((data$exp_count>q25)^2),formatC(probi2,format="f",digits=16), formatC(p2,format="f",digits=16))
				colnames(data)=c('chromosome','start','stop','q25','peakprob','qvalue') #all non-zero set to 1
		   }
		   if(diff==0){
				#if differential expression (comparing two samples) is not occuring
				data$q25[Y-mui1<0]=0
		   }
		   line = paste("\nProcessed ",files[fnum],"\n\t","Selected number of peaks: ",as.character(numpeaks),"\n\t",as.character(Sys.time()),"\t",sep='')
				### PRINT WINDOWS
				cat(line)
				winfile = paste(winout,"_",chrm,".wins",sep="")
				if(printflag==1){
					write.table(data,winfile,quote=F,sep="\t",row.names=F,col.names=F,append=T)
				}else{
					write.table(data,winfile,quote=F,sep="\t",row.names=F)
					printflag=1
				}
				fnum=fnum+1
				rm(data); rm(Y); rm(X); rm(XE); rm(XZ); rm(probi0); rm(probi1); rm(probi2); 
				rm(mui1); rm(mui2); rm(start); rm(prop1); rm(p);rm(p2); gc();
			}
		}
	}
	time.end <- Sys.time()
	cat("\n")
	print(difftime(time.end,time.start))
	return(winfile)
#final
}
