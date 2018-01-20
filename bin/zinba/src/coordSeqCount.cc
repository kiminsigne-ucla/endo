#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "canalysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void coordSeqCount(char **Rsbpcfile,char ** Rcoordfile,char **Routputfile, char **Rtwobitfile){
	const char*  sbpcFile = Rsbpcfile[0];
	const char* coordfile = Rcoordfile[0];
	const char* outputFile = Routputfile[0];
	const char* twobitfile = Rtwobitfile[0];
	
	canalysis newAnalysis;// = new analysis;
	int ret=newAnalysis.importCoords(coordfile);
	
	if(ret == 0){
		Rprintf("Getting basecount data from %s\n",sbpcFile);
		ret = newAnalysis.processCoords(sbpcFile,outputFile,twobitfile);
	}
	
	if(ret != 0){
		Rprintf("\ncoordSeqCount EXITING with ERRORS\n");		
	}else if (ret == 0){
		Rprintf("\ncoordSeqCount SUCCESSFULLY COMPLETE\n");
	}

}
}
