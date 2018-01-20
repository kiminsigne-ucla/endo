#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "frpanalysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C"{
void full_refinepeaks(char **RreadFile,char **RbcountFile,int *RrextendLength,char **Rrfiletype,char **Rwinlist,double *Rthreshold,char **Rmethod,int *Rwformat,char **Rchromosome,int *RwinGap,int *RFDR,char **Rtwobitfile,char **RoutputFile){

	const char * readFile = RreadFile[0];
	const char * bcountFile = RbcountFile[0];
	const char * twobitfile = Rtwobitfile[0];
	int extendLength = RrextendLength[0];
	const char * filetype = Rrfiletype[0];

	const char* winlist = Rwinlist[0];	
	double threshold = Rthreshold[0];
	const char* method = Rmethod[0];
	int wformat = Rwformat[0];
	int winGap = RwinGap[0];
	int FDR=*RFDR;
	const char* chromosome = Rchromosome[0];	
	
	const char * outputFile = RoutputFile[0];
	frpanalysis newAnalysis;// = new analysis;

	Rprintf("\nGetting significant windows from %s\n",winlist);
	Rprintf("\tThreshold is %f\n",threshold);
	if(FDR == 0)
		Rprintf("\tUsing posterior probability as threshold\n");
	else if (FDR == 1)
		Rprintf("\tUsing qvalue as threshold\n");
	
	Rprintf("\tData generated using %s\n",method);
	
	if(wformat == 0)
		Rprintf("\tData format is compact\n");
	else if (wformat == 1)
		Rprintf("\tData format is expanded\n");
	
	Rprintf("\tDistance to collapse windows is %i\n",winGap);
	
	int ret=newAnalysis.importCoords(winlist,threshold,method,wformat,winGap,FDR,twobitfile);
	
	if(ret == 0){
		string bcflag = "none";
		if(strcmp(bcountFile,bcflag.c_str()) == 0){
			Rprintf("\nImporting basecount data from %s ....\n",bcountFile);
		}else{
			Rprintf("\nImporting reads from file %s ....\n",readFile);
			Rprintf("\tReads are formatted as %s ....\n",filetype);
			Rprintf("\tExtending reads by %i bp....\n",extendLength);
			ret=newAnalysis.importRawSignal(readFile,extendLength,filetype);
		}

		if(ret == 0){
			Rprintf("\nGetting basecount data for significant coords\n");
			ret = newAnalysis.processSigCoords(bcountFile,extendLength,twobitfile,chromosome,outputFile);
			if(ret == 0)
				Rprintf("-------- REFINE PEAKS COMPLETED SUCCESSFULLY --------\n");
		}
	}
		
	if(ret != 0)
		Rprintf("-------- REFINE PEAKS EXITING WITH ERRORS --------\n");
}
}
