#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "winCounts.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void readcountWindows(char **RseqFile,char **RcoordFile,char **Routfile,char **RtwoBitFile,int *Rextension){

	string seqFile = RseqFile[0];
	string coordFile = RcoordFile[0];
	string outfile = Routfile[0];
	const char * twoBitFile = RtwoBitFile[0];
	int extension = Rextension[0];
	
	int ret;
	winCounts newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s \n", seqFile.c_str());
	ret=newAnalysis.importRawSignalWC(seqFile.c_str(),extension,twoBitFile);

	if(ret == 0){
		Rprintf("\nImporting coordinates from file %s \n", coordFile.c_str());		
		ret = newAnalysis.importCoordsWC(coordFile.c_str());
		if(ret == 0){
			Rprintf("\nGetting reads counts in windows\n");
			ret = newAnalysis.processSignalsWC(twoBitFile,outfile.c_str());
		}
	}

	if(ret != 0){
		Rprintf("\n\n--------READ COUNT IN WINDOWS EXITING WITH ERRORS----------------\n\n");
	}else if (ret == 0){
		Rprintf("\n\n--------READ COUNT IN WINDOWS COMPLETE-------\n\n");
	}
}
}
