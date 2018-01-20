rc.wins=function(alignreads,coordfile,outputfile,twobitfile,extension){
    cReturn <- .C("readcountWindows",as.character(alignreads),as.character(coordfile),as.character(outputfile),as.character(twobitfile),as.integer(extension),PACKAGE="zinba")
}
