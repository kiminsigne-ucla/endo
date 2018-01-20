#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "calcCovs.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void getWindowCounts(char **RexpSeqFile,char **RtwoBitFile,int *RzWinSize,int *RzOffsetSize, char **Rfiletype,int *Rextension, double *RNthresh){

	string expSeqFile = RexpSeqFile[0];
	const char * twoBitFile = RtwoBitFile[0];
	int zWinSize = RzWinSize[0];
	int zOffsetSize = RzOffsetSize[0];
	const char * filetype = Rfiletype[0];
	int extension = Rextension[0];
	double Nthresh = *RNthresh;
	int ret;
	
	Rprintf("\nRunning getWindowCounts\n");
	Rprintf("\nWindow size is %i\n",zWinSize);
	Rprintf("\nOffset size is %i\n",zOffsetSize);	
	Rprintf("\nExtension is %i\n",extension);
	Rprintf("\nN threshold is %f\n",Nthresh);
	Rprintf("\nFiletype is %s\n",filetype);
	
	calcCovs newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s \n", expSeqFile.c_str());
	import b;
	ret=newAnalysis.importRawSignal(expSeqFile.c_str(),extension,filetype,0,twoBitFile);

	if(ret == 0){
		Rprintf("\nBuilding window data\n");
		size_t found = expSeqFile.find_last_of(".");
		string outfile_prefix = expSeqFile.substr(0,found);
		ret = newAnalysis.processWinSignal(zWinSize,zOffsetSize,twoBitFile,outfile_prefix,extension,filetype, Nthresh);
	}
	
	if(ret != 0){
		Rprintf("ERROR: building windows was unsuccssful\n");
	}else if (ret == 0){
		Rprintf("\n\n--------GET WINDOW COUNTS COMPLETE-------\n\n");
	}
}
}
