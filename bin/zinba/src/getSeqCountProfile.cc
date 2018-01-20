#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include "analysis.h"
#include <ext/slist>
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

extern "C" {
void getSeqCountProfile(char **Rinputfile,char **Rwinlist,double *Rthreshold,char **Rmethod,int *Rwformat,char **Routputfile, char **Rtwobitfile,char **Rchromosome, int *RwinGap, int *RFDR, int *Rbinary){

	const char*  inputFile = Rinputfile[0];
	const char* winlist = Rwinlist[0];
	const char* outputFile = Routputfile[0];
	double threshold = Rthreshold[0];
	const char* method = Rmethod[0];
	int wformat = Rwformat[0];
	int winGap = RwinGap[0];
	int FDR=*RFDR;
	const char* twobitfile = Rtwobitfile[0];
	const char* chromosome = Rchromosome[0];
	int binary = Rbinary[0];
	
	analysis newAnalysis;// = new analysis;
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
	
	int ret=newAnalysis.importCoords(winlist,threshold,method,wformat,winGap,FDR);
	
	if(ret == 0){
		
		Rprintf("\nGetting basecount data for %s\n",chromosome);
		ret = newAnalysis.processCoords(inputFile,outputFile,twobitfile,chromosome, binary);
	}
	
	if(ret != 0)
		error("\ngetSeqCountProfile exiting with ERRORS\n");
	else if (ret == 0)
		Rprintf("\ngetSeqCountProfile COMPLETED SUCCESSFULLY\n");

}
}
