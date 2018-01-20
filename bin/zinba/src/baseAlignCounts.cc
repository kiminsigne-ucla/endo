#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "bcanalysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C"{
void baseAlignCounts(char **RinputFile,char **RoutputFile, char **Rtwobitfile,int *RextendLength,char **Rfiletype, int *Rbinary){
	const char * inputFile = RinputFile[0];
	const char * outputFile = RoutputFile[0];
	const char * twobitfile = Rtwobitfile[0];
	int extendLength = RextendLength[0];
	int binary = Rbinary[0];
	const char * filetype = Rfiletype[0];
	
	bcanalysis newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s ....\n",inputFile);
	Rprintf("\tReads are formatted as %s ....\n",filetype);
	Rprintf("\tExtending reads by %i bp....\n",extendLength);
	int ret=newAnalysis.importRawSignal(inputFile,extendLength,filetype,twobitfile);
	if(ret == 0){
		Rprintf("\nCalculating counts at each base\n");
		ret = newAnalysis.processSignals(outputFile,extendLength, binary);
		if(ret == 0)
			Rprintf("-------- BASE ALIGN COUNTS COMPLETE SUCCESSFULLY --------\n");
	}
		
	if(ret != 0)
		Rprintf("-------- BASE ALIGN COUNTS EXITING WITH ERRORS --------\n");
}
}
