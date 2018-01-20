twobitinfo=function(infile, outfile){
	.C("twoBitInfo",as.character(infile), as.character(outfile),PACKAGE="zinba")
}
	
