startenrichment=function(range, data, formula,formulaE, formulaZ, initmethod, trace=0){
	library(zinba)
	suppressPackageStartupMessages(library(quantreg))
            mf <- model.frame(formula=formula, data=data)
	    mfE <- model.frame(formula=formulaE, data=data)
	    mfZ <- model.frame(formula=formulaZ, data=data)
            X <- model.matrix(attr(mf, "terms"), data=mf)
	    XE <- model.matrix(attr(mfE, "terms"), data=mfE)
	    XZ <- model.matrix(attr(mfZ, "terms"), data=mfZ)
	    XNB=as.data.frame(X[,-c(1)])
	    XNBE=as.data.frame(XE[,-c(1)])
	    XNBZ=as.data.frame(XZ[,-c(1)])
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
                    loglik0=log(prob0+(1-prob0)*prop1*dnbinom(0, size = theta1, mu = mu1)+(1-prob0)*prop2*dnbinom(0, size = theta2, mu = mu2))
                    loglik1=log((1-prob0)*prop1*dnbinom(Y, size = theta1, mu = mu1)+(1-prob0)*prop2*dnbinom(Y, size = theta2, mu = mu2))
                    NAs=which(loglik1==-Inf)
                    if(length(NAs>0)){
                            loglik1[NAs]=apply(cbind(log((1-prob0[NAs])*prop1)+dnbinom(Y[NAs], size = theta1, mu = mu1[NAs],log=TRUE),log((1-prob0[NAs])*prop2)+dnbinom(Y[NAs], size = theta2, mu = mu2[NAs], log=TRUE)), 1, logsumexp)
                    }
                    loglik=sum(loglik0*Y0+loglik1*Y1)
                    loglik
            }
            Y <- model.response(mf)
            n <- length(Y)
            kx <- NCOL(X)
            Y0 <- Y <= 0
            Y1 <- Y > 0
            linkstr <- 'logit'
            linkobj <- make.link(linkstr)
            linkinv <- linkobj$linkinv


	probs=seq(range[1], range[2],,5)
	result=rep(0, length(probs))
	a=Sys.time()


	for(k in 1:length(probs)){
           model_zero <-.C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), M=as.integer(ncol(XNBZ)), y=as.double(Y0), prior=as.double(rep(1,n)), offset=as.double
(rep(0,n)), X=as.double(unlist(XNBZ)),  stratum=as.integer(rep(1,n)),init=as.integer(1), rank=integer(1), Xb=double(n*ncol(XNBZ)), fitted=as.double((rep(1,n) * Y0 + 0.5)/(rep
(1,n) + 1)), resid=double(n), weights=double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')


		#starting params for count componenets
	   if(initmethod=='quantile'){
	        prop2=probs[k]
	  	prop1=1-prop2
		data2=data
	        if(sum(colnames(data)=='input_count')==1){data2$input_count=exp(data2$input_count)-1}
		if(sum(colnames(data)=='exp_cnvwin_log')==1){data2$exp_cnvwin_log=exp(data2$exp_cnvwin_log)-1}
		prop2=startprop
            	prop1=1-prop2
		t=rq(formula, tau=1-prop2, data=data2, method='pfn')
		priorCOUNTweight=rep(10^-10, length(Y))      
		priorCOUNTweight[as.double(which(t$residuals>quantile(t$residuals,1-prop2)))]=1-10^-10
		rm(data2)
  	   }else if(initmethod=='count'){
	        prop2=probs[k]
   	  	prop1=1-prop2
	        n1  = round(length(Y) * (1 - prop2))
		priorCOUNTweight=rep(1-10^-10, length(Y))
                odY = order(Y)
		priorCOUNTweight[odY[1:n1]]=10^-10
  	  }

	    model_count1 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(1-priorCOUNTweight), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(Y+(Y==0)/6), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')
            model_count2 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), M=as.integer(ncol(XNBE)), y=as.double(Y), prior=as.double(priorCOUNTweight), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNBE)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNBE)), fitted=as.double(Y+(Y==0)/6), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')
            
            #starting prior probs
            mui1  <- model_count1$fitted
            mui2  <- model_count2$fitted
            probi0 <- model_zero$fitted

            #start mean vector
            start <- list(count1 = model_count1$fitted, count2 = model_count2$fitted,zero = model_zero$fitted, zerocoef=model_zero$coefficients)
            start$theta1 <- 1
            start$theta2 <- 1
                        
            probi1  <- (1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)/(probi0*Y0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
            probi2  <- (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2)/(probi0*Y0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
            probi0=probi0/(probi0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
            probi0[Y1]=0
            
            NAs=which(probi1=='NaN'| probi2=='NaN')			
            if(length(NAs>0)){
                    probi1[NAs]=0
                    probi2[NAs]=1
            }	
            ll_new <- loglikfun(list(start=start, prop1=prop1, prop2=prop2))
            ll_old <- 2 * ll_new      

            if(!require("MASS")) {
                    ll_old <- ll_new
                    warning("EM estimation of starting values not available")
            }
            ll=matrix(0, 1, 10000)
            ll[1]=ll_new
            i=2

		while(i<3) {
		#print(ll_new)
                ll_old <- ll[max(1, i-10)]
                prop1=sum(probi1)/(sum(probi1)+sum(probi2))
                prop2=sum(probi2)/(sum(probi1)+sum(probi2))
		if(prop1<.5){
			break
		}
                #updated values for parameters of component means
                 
                model_zero <- .C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), M=as.integer(ncol(XNBZ)), y=as.double(probi0), prior=as.double(rep(1,n)), 
offset=as.double(rep(0,n)), X=as.double(unlist(XNBZ)),  stratum=as.integer(rep(1,n)),init=as.integer(1), rank=integer(1), Xb=double(n*ncol(XNBZ)), fitted=as.double((rep(1,n)*
probi0 + 0.5)/(rep(1,n) + 1)), resid=double(n), weights=double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')    
                model_count1 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(probi1), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(start$count1), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta1), package='zinba')  
                model_count2 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), M=as.integer(ncol(XNBE)), y=as.double(Y), prior=as.double(probi2), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNBE)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNBE)), fitted=as.double(start$count2), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta2), package='zinba')  

                start <- list(count1 = model_count1$fitted, count2 = model_count2$fitted,zero = model_zero$fitted, zerocoef=model_zero$coefficients)
                start$theta1 <- model_count1$theta
                start$theta2 <- model_count2$theta

                mui1  <- model_count1$fitted
                mui2  <- model_count2$fitted
                probi0 <- model_zero$fitted

		if(any(!is.finite(mui1)) | any(!is.finite(mui2)) | any(!is.finite(probi0))){
                        fail=1
                       	break
                }
	        probi1  <- (1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)/(probi0*Y0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
            	probi2  <- (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2)/(probi0*Y0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
                probi0=probi0/(probi0+(1-probi0)*prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ (1-probi0)*prop2*dnbinom(Y, size = start$theta2, mu = mui2))
                probi0[Y1]=0
		NAs=which(probi1=='NaN'| probi2=='NaN')			
                if(length(NAs>0)){
                        probi1[NAs]=0
                        probi2[NAs]=1
                }

                ll_new <- loglikfun(list(start=start, prop1=prop1, prop2=prop2))
                ll[i]=ll_new
                i=i+1 
		cat(".")
		}
		result[k]=ll_new
	}
	if(trace==1) print(result)
	return(probs[which.max(result)])
	rm(data); rm(Y); rm(X); rm(XNB); rm(XZ); rm(XNBZ); rm(XE);rm(XNBE);rm(probi0); rm(probi1); rm(probi2); rm(mui1); rm(mui2); rm(start); rm(prop1); rm(prop2);gc();
}
