signalnoise=function(inputfile,twoBitFile,winSize=100000){
	if(!file.exists(inputfile)){
    		stop(paste("Input file not found,",inputfile,sep=" "))
    	}
    	if(!file.exists(twoBitFile)){
        	stop(paste("twoBit file not found,",twoBitFile,sep=" "))
    	}
	output='tempsignalnoise.out'
	c=.C("signalnoise",as.character(inputfile),as.character(output),as.character(twoBitFile),as.integer(winSize),PACKAGE="zinba")
	a=read.table(output)
	unlink(output)
	median=a[,1]
	max=a[,2]	
	print("Summary Information of ratio of Max Window Count vs Median Window Count")
	print(summary((max/median)[median>0]))
}
	
