collapsewindows=function(winlist,thresholds,method='pscl',printFullOut=0, winGap=0, FDR=NULL, output=NULL){
	if(!file.exists(winlist)){stop("Specified winlist doesnt exist or is incorrect")}
	if(is.null(FDR)) stop("FDR usage to threshold windows needs to be specified as TRUE or FALSE")
	if(is.null(output))	print("no output path specified for unrefined peak file, output place in _files/ directory")
	cReturn <- .C("collapse_windows",as.character(winlist),as.character(method),as.integer(printFullOut),as.integer(length(thresholds)),as.double(thresholds), winGap=as.integer(winGap), as.integer(FDR^2), 
as.character(output),PACKAGE="zinba")
}
