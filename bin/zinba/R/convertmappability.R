convertmappability=function(inputfile, outputfile){
    .C("mapout2alignability",as.character(inputfile),as.character(outputfile),PACKAGE="zinba")
}
