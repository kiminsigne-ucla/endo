hgfakeagp=function(faFile,outFile){
    .C("hgFakeAgp",as.character(faFile),as.character(outFile),PACKAGE="zinba")
}
	
