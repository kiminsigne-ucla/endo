
covariateselect=function(file=NULL, covs=NULL, selection="dirty",formulaZ=NULL,numProc=1,loc=NULL, start=1, append=FALSE, interaction=TRUE){

if(!is.null(file)){
	if(!file.exists(file)) cat("covselect.R: file at path", file, "does not exist")
}

if(is.null(covs)){
	cat("covselect.R: covariates must be specified for model selection, see zinba website for details")
}




###################################local functions#####################################3
#From MuMIn R Package
`getAllTerms.default` <-
function(x, ...) getAllTerms(as.formula(formula(x)), ...)


`getAllTerms.terms` <-
function(x, offset = TRUE, ...) {
	if (!is.null(attr(x, "offset"))){
		offs <- sapply((attr(x, "variables")[-1])[attr(x, "offset")], deparse)
	} else {
		offs <- NULL
	}
	ret <- attr(x, "term.labels")

	# Get term names, with higher order term components arranged alphabetically
	if (length(ret) > 0) {
		factors <- attr(x, "factors")
		factors1 <- rownames(factors)
		ret <- apply(factors > 0, 2, function(i) paste(sort(factors1[i]), collapse=":"))
	}

	# Leave out random terms (lmer type)
	#ran <- attr(x, "variables")[-1][-c(attr(x, "offset"), attr(x, "response"))]
	ran <- attr(x, "variables")[-1]
	ran <- as.character(ran[sapply(ran,
		function(x) length(x) == 3 && x[[1]] == as.name("|"))])
	ifx <- !(ret %in% ran)

	ret <- ret[ifx]

	# finally, sort by order and then alphabetically
	ret <- unname(ret[order(attr(x, "order")[ifx], ret)])


	if (!is.null(offs[1])) {
		if (offset)
			ret <- c(ret, offs)
		attr(ret, "offset") <- offs
	}
	attr(ret, "intercept") <- attr(x, "intercept")
	if (length(ran) > 0) {
		attr(ret, "random.terms") <- ran
		attr(ret, "random") <- reformulate(c(".", paste("(", ran, ")",
			sep = "")), response = ".")
	}

	return(ret)
}



`getAllTerms.formula` <-
function(x, ...) getAllTerms.terms(terms(x), ...)


`getAllTerms.lme` <-
function(x, ...) {
	ret <- getAllTerms(terms(x))
	attr(ret, "random") <- . ~ .

	# Code from nlme:::print.reStruct, modified slightly
	reStruct <- x$modelStruct$reStruct
	nobj <- length(reStruct)
	if (is.null(namx <- names(reStruct)))
		names(reStruct) <- nobj:1
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	reStruct[] <- rev(reStruct)
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	attr(ret, "random.terms") <- paste(lapply(lapply(reStruct, attr, "formula"),
		"[[", 2), "|",
		rev(apply(aux, 1, function(z) paste(z[z != ""], collapse = " %in% "))))

	return(ret)
}

`getAllTerms.glmer` <- # For backwards compatibility
`getAllTerms.lmer` <-  # with older versions of lme4
`getAllTerms.mer` <-
function(x, ...) getAllTerms(formula(x), ...)

`getAllTerms` <-
function(x, ...) UseMethod("getAllTerms")

`formulaAllowed` <-
function(frm, except=NULL) {
	if(isTRUE(except)) return(TRUE)
	factors <- attr(terms(frm), "factors")
	if(length(factors) == 0) return(TRUE)
	if(is.character(except))
		factors <- factors[!(rownames(factors) %in% except), ]
	return(all(factors < 2))
}

covanal=function(file, formula, formulaE, formulaZ,i,loc, size){
	#k is sum of the number of covariates within each component + the rest is fixed regardless of the model (intercepts always used in each component + two dispersion parameters)
	k=sum(attr(terms(formula),"order"))+sum(attr(terms(formulaE),"order"))+sum(attr(terms(formulaZ),"order"))+ (3+2)
	mix0=getsigwindows(file=file,formula=formula,formulaE=formulaE,formulaZ=formulaZ,method='mixture',initmethod="count",modelselect=TRUE)	
	write.table(t(c(i,as.character(formula)[3], as.character(formulaE)[3], as.character(formulaZ)[3],-2*mix0$ll+k*mix0$logdimdata, 2*k-2*mix0$ll,-2*mix0$ll+sum(-2*mix0$probi0*log(mix0$probi0 + (mix0$probi0==0))-2*mix0$probi1*log(mix0$probi1 + (mix0$probi1==0))-2*mix0$probi2*log(mix0$probi2 + (mix0$probi2==0)))+k*mix0$logdimdata,k,mix0$ll, mix0$fail)),file=loc,quote=F, append=T, row.names=F,col.names=F, sep="\t")
	rm(mix0)		
	gc()
	cat(paste("\nmodel ",i," completed out of ",size,"\t"))
	return(1)
}


########################################################################################
#Modified from dredge.R in MuMIn package
require(multicore)
require(doMC)
require(foreach)
library(MASS)
library(R.utils)

supported=c("gcPerc", "align_perc", "input_count", "exp_cnvwin_log")
if( sum(covs %in% supported)!=length(covs)){
	stop(paste("Covariate ",covs[which(covs %in% supported == FALSE)], " not found in list of supported covariates"))
}
marg.ex=NULL
fixed=as.formula("~1")
if(interaction==TRUE){
        global.model= as.formula(paste("exp_count~", paste(covs, collapse="*")))
}else if(interaction==FALSE){
        global.model= as.formula(paste("exp_count~", paste(covs, collapse="+")))
}

all.terms <- getAllTerms(global.model)
n.vars=length(all.terms)
m.max <- n.vars
gterms <- tryCatch(terms(formula(global.model)),
		error=function(...) terms(global.model))

if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1]] != "~" || length(fixed) != 2)
				warning("'fixed' should be a formula of form: ",
						"~ a + b + c")
			fixed <- c(getAllTerms(fixed))
		} else if (!is.character(fixed)) {
			stop ("'fixed' should be either a character vector with"
				  + " names of variables or a one-sided formula")
		}
		if (!all(fixed %in% all.terms)) {
			warning("Not all terms in 'fixed' exist in 'global.model'")
			fixed <- fixed[fixed %in% all.terms]
		}
		m.max <- m.max - length(fixed)
	}

	if (m.max > 0) {
		num.opt.vars <- 1:n.vars
		if (!is.null(fixed))
			num.opt.vars <- num.opt.vars[!(all.terms %in% fixed)]

		all.comb <- lapply(seq(m.max), combn, x = num.opt.vars)
		all.comb <- unlist(lapply(all.comb, function(.x) split(.x, col(.x))),
						   recursive = FALSE)
		all.comb <- c(`0` = list(0), all.comb)

	}else{
		all.comb <- list(0)
	}

	if (!is.null(fixed))
		all.comb <- lapply(all.comb, append, (1:n.vars)[all.terms %in% fixed])

	int.term <- 1

	formulas <- lapply(all.comb,
		function(.x) reformulate(c(all.terms[.x], int.term), response = "exp_count" ))

	env <- attr(gterms, ".Environment")
	formulas <- lapply(formulas, `attr<-`, ".Environment", env)

	ss <- sapply(formulas, formulaAllowed, except = marg.ex)

	all.comb <- all.comb[ss]
	formulas <- formulas[ss]

	names(formulas) <- seq(formulas)
	if (any(inherits(global.model, c("mer", "lmer", "glmer")))) {
          formulas <- lapply(formulas, update, attr(all.terms, "random"))
	}

print(formulas)
if(selection=="complete"){
	index=start:(length(formulas)^3)
		write.table(t(c("i","formula", "formulaE", "formulaZ","BIC","AIC","ICL","k","ll", "fail")),loc,quote=F, append=F,  row.names=F, 	col.names=F, sep="\t")
 		registerDoMC(numProc)
 		mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
 		getDoParWorkers()
 		result <- foreach(i=index,.combine='rbind',.inorder=FALSE,.errorhandling="remove",.options.multicore = mcoptions) %dopar%
                     covanal(file=file,formula=formulas[[ceiling(i/length(formulas)^2)]],formulaE=formulas[[ceiling(i/length(formulas))-((ceiling(i/length(formulas)^2))-1)*length(formulas)]],formulaZ=formulas[[i-ceiling(i/length(formulas)-1)*length(formulas)]],i=i,loc=loc, size=length(formulas)^3)

		
}else if(selection=="dirty"){
	if(is.null(formulaZ)) formulaZ=exp_count~1 
	index=start:(length(formulas)^2)
	write.table(t(c("i","formula", "formulaE", "formulaZ","BIC","AIC","ICL","k","ll","fail")),loc,quote=F, append=F, row.names=F, 	col.names=F, sep="\t")
 	registerDoMC(numProc)
 	mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
 	getDoParWorkers()
 	result <- foreach(i=index,.combine='rbind',.inorder=FALSE,.errorhandling="remove",.options.multicore = mcoptions) %dopar%{
                     covanal(file=file,formula=formulas[[ceiling(i/length(formulas))]],
											formulaE=formulas[[i-(ceiling(i/length(formulas))-1)*length(formulas)]],
											formulaZ=formulaZ,i=i,loc=loc, size=length(formulas)^2)
			}
	
	size=file.info(loc)$size
	fields=count.fields(loc, sep="\t")
	lines=countLines(loc)

	#zero file size?
	if(size == 0) cat("Intermediate (dirty) model file", loc, "has size 0.  Please check if you have run out of disk space")

	#unequal number of columns?
	if(length(table(fields)) != 1){
		cat("Warning: intermediate (dirty) model file", loc, "has", length(skip), "lines with less than expected columns.  Removing problematic lines")
	}
	skip=which(fields != max(fields))

	#no lines written to output?
	if(lines == 0) cat("Intermediate (dirty) model file", loc, "has 0 lines written to it.  This suggests anomaly in model selection, please contact the package administrator")

	#now pick best zero-inflated formula given bg and enrichd covariates
	if(length(skip)==0){
		final=read.table(loc, header=T, sep="\t")
	}else{
		final_lines=readLines(loc)[-skip]
		final = read.table(textConntection(final_lines), header=T, sep="\t") 
	}

	#make sure that all models did not fail
	if(all(final$fail==1))  stop("covselect.R:  all models failed after model selection intermediate step")

	final=final[final$fail==0,]
	bestBIC=which.min(final$BIC)
	
	a=try(as.formula(paste("exp_count~",final$formula[bestBIC])), silent=T)
	b=try(as.formula(paste("exp_count~",final$formulaE[bestBIC])), silent=T)
	c=try(as.formula(paste("exp_count~",final$formulaZ[bestBIC])), silent=T)

	#check for proper form of each formula, errors would not result from .model file having 0 size or 0 lines
	if(inherits(a, "try-error") ==T | inherits(b, "try-error") ==T | inherits(b, "try-error") ==T) stop("Intermediate (dirty) formulas not read properly from .model file, check .model file in _files/directory for errors")

	formvector=c(
		as.formula(paste("exp_count~",final$formula[bestBIC])),
		as.formula(paste("exp_count~",final$formulaE[bestBIC])),
		as.formula(paste("exp_count~",final$formulaZ[bestBIC]))
	)
	#write.table(t(c("i","formula", "formulaE", "formulaZ","BIC","AIC","ICL","k","ll")),loc,quote=F, append=F, row.names=F, 	col.names=F, sep="\t")
	index=start:(length(formulas))	
	result <- foreach(i=index,.combine='rbind',.inorder=FALSE,.errorhandling="remove",.options.multicore = mcoptions) %dopar%
                     covanal(file=file,formula=formvector[[1]],formulaE=formvector[[2]],formulaZ=formulas[[i]],i=i,loc=loc, size=length(formulas))
	
	
}else{
	stop("selection parameter must be specified as either 'dirty' or 'complete'")
} 

#check for errors in writing process
#file size zero due to running out of disk space

size=file.info(loc)$size
fields=count.fields(loc, sep="\t")
lines=countLines(loc)

#zero file size?
if(size == 0) cat("Final model file", loc, "has size 0.  Please check if you have run out of disk space")

#unequal number of columns?

if(length(table(fields)) != 1){
	cat("Warnings: final model file", loc, "has", length(skip), "lines with less than expected columns.  Removing problematic lines")
}
skip=which(fields != max(fields))

#no lines written to output?
if(lines == 0) cat("Final model file", loc, "has 0 lines written to it.  This suggests anomaly in model selection, please contact the package administrator")


if(length(skip)==0){
	final=read.table(loc, header=T, sep="\t")
}else{
	final_lines=readLines(loc)[-skip]
	final = read.table(textConntection(final_lines), header=T, sep="\t") 
}

if(all(final$fail==1))  stop("covselect.R:  all models failed after model selection final step")
final=final[final$fail==0,]
bestBIC=which.min(final$BIC)

a=try(as.formula(paste("exp_count~",final$formula[bestBIC])), silent=T)
b=try(as.formula(paste("exp_count~",final$formulaE[bestBIC])), silent=T)
c=try(as.formula(paste("exp_count~",final$formulaZ[bestBIC])), silent=T)


#check for proper form of each formula, errors would not result from .model file having 0 size or 0 lines
if(inherits(a, "try-error") ==T | inherits(b, "try-error") ==T | inherits(b, "try-error") ==T) stop("Final formulas not read properly from .model file, check .model file in _files/directory for errors")

formvector=c(
	as.formula(paste("exp_count~",final$formula[bestBIC])),
	as.formula(paste("exp_count~",final$formulaE[bestBIC])),
	as.formula(paste("exp_count~",final$formulaZ[bestBIC]))
)

return(formvector)
}

