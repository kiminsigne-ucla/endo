#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "calcCovs.h"
#include "import.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025
#define max(x,y) ((x) > (y) ? (x) : (y))

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void buildWindowsfast(char **RexpSeqFile,char **RinSeqFile,char **RtwoBitFile,int *RzWinSize,int *RzOffsetSize, int *RcWinSize, int *RcOffsetSize, char **Rfiletype,char **Rfilelist,int *Rextension,char **Routdir){

	string expSeqFile = RexpSeqFile[0];
	const char * inSeqFile = RinSeqFile[0];
	//string alignDir = RalignDir[0];
	const char * twoBitFile = RtwoBitFile[0];
	int zWinSize = RzWinSize[0];
	int zOffsetSize = RzOffsetSize[0];
	int cWinSize = RcWinSize[0];
	int cOffsetSize = RcOffsetSize[0];
	//int binary = Rbinary[0];
	const char * filetype = Rfiletype[0];
	const char * filelist = Rfilelist[0];
	const char * outdir = Routdir[0];
	const char * defaultdir = "default";
	int extension = Rextension[0];
	
	int ret;
	calcCovs newAnalysis;// = new analysis;
	Rprintf("\nImporting reads from file %s \n", expSeqFile.c_str());
	Rprintf("\tFiletype is %s \n", filetype);
	Rprintf("\tExtension is %i \n", extension);
	ret=newAnalysis.importRawSignal(expSeqFile.c_str(),extension,filetype,0,twoBitFile);

	if(ret == 1){
		error("\nEXITING due to error in importing mapped reads\n");
	}else{ 
		size_t found;
		size_t found2;
		string outfile_prefix;
		Rprintf("\nBuilding window data\n");
		if(strcmp(defaultdir,outdir)==0){
			found = expSeqFile.find_last_of(".");
			outfile_prefix = expSeqFile.substr(0,found);
		}else{
			found = expSeqFile.find_last_of(".");
			found2 = expSeqFile.find_last_of("/");
			outfile_prefix = outdir+expSeqFile.substr(found2+1,found-found2-1);
		}		
		//ret = newAnalysis.processSignals(zWinSize,zOffsetSize,cWinSize,cOffsetSize,alignDir,twoBitFile,inSeqFile,outfile_prefix,filelist,extension,filetype, binary);
      ret = newAnalysis.processSignalsfast(zWinSize,zOffsetSize, cWinSize, cOffsetSize, twoBitFile,inSeqFile,outfile_prefix,filelist,extension,filetype);
		if(ret == 1){
			error("ERROR: building windows was unsuccssful\n");
		}
	}
	Rprintf("\n\n--------BUILD WINDOWS COMPLETE-------\n\n");
}
}
