#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <fstream>
#include <vector>
#include "process_bin_scores.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#define MAX_LEN 1025

using namespace std;

extern "C" {
void binAlignAdjust(char **Rfilelist, char **Routdir,char **Rtwobitfile,int *RaThresh, int *RadjustSize){

	const char * filelist = Rfilelist[0];
	string outDir = Routdir[0];
	const char * twoBit = Rtwobitfile[0];
	int aThresh = RaThresh[0];
	int adjustSize = RadjustSize[0];

	process_bin_scores newProcess;
	
	int ret = newProcess.adjustCoords(filelist,outDir,twoBit,aThresh,adjustSize);
	if(ret != 0)
		error("\n---------------- ALIGN ADJUST EXITING WITH ERRORS ----------------\n");
	//else
		//Rprintf("\n---------------- ALIGN ADJUST COMPLETED SUCCESSFULLY ----------------\n");
}
}
