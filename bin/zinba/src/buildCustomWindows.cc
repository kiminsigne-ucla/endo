/*#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "calcCovs.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025
#define max(x,y) ((x) > (y) ? (x) : (y))

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void buildCustomWindows(char **RexpSeqFile,char **RinSeqFile,char **RtwoBitFile,int *RzWinSize,int *RzOffsetSize,char **Rfiletype,char **Rfilelist,int *Rextension,char **Routdir){

	string expSeqFile = RexpSeqFile[0];
	const char * inSeqFile = RinSeqFile[0];
	const char * twoBitFile = RtwoBitFile[0];
	int zWinSize = RzWinSize[0];
	int zOffsetSize = RzOffsetSize[0];
	const char * filelist = Rfilelist[0];
	const char * outdir = Routdir[0];
	const char * defaultdir = "default";
	int extension = Rextension[0];
	const char * filetype = Rfiletype[0];

	int ret;
	calcCovs newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s \n", expSeqFile.c_str());
	Rprintf("\tFiletype is custom \n");
	Rprintf("\tExtension is %i \n", extension);
	ret=newAnalysis.importRawSignal(expSeqFile.c_str(),extension,filetype,0,twoBitFile);

	if(ret == 1){
		Rprintf("\nEXITING due to error\n");
	}else{
		size_t found;
		size_t found2;
		string outfile_prefix;
		Rprintf("\nBuilding custom window data\n");
		if(strcmp(defaultdir,outdir)==0){
			found = expSeqFile.find_last_of(".");
			outfile_prefix = expSeqFile.substr(0,found);
		}else{
			found = expSeqFile.find_last_of(".");
			found2 = expSeqFile.find_last_of("/");
			outfile_prefix = outdir+expSeqFile.substr(found2+1,found-found2-1);
		}		
		ret = newAnalysis.processCustomSignal(zWinSize,zOffsetSize,twoBitFile,inSeqFile,outfile_prefix,filelist,extension);
		if(ret == 1){
			Rprintf("ERROR: building windows was unsuccssful\n");
		}
	}
	Rprintf("\n\n--------BUILD WINDOWS COMPLETE-------\n\n");
}
}*/
