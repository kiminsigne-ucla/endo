#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "process_scores.h"
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

using namespace std;

extern "C" {
void alignAdjust(char **Rinputfile, char **Routdir,char **Rtwobitfile,int *RaThresh, int *RadjustSize){
	string inputFile = Rinputfile[0];
	string outDir = Routdir[0];
	const char * twoBit = Rtwobitfile[0];
	int aThresh = RaThresh[0];
	int adjustSize = RadjustSize[0];

	process_scores newProcess;
	int ret = newProcess.adjustCoords(inputFile,outDir,twoBit,aThresh,adjustSize);
	if(ret != 0)
		Rprintf("\n---------------- ALIGN ADJUST EXITING WITH ERRORS ----------------\n");
	else
		Rprintf("\n---------------- ALIGN ADJUST COMPLETED SUCCESSFULLY ----------------\n");
}
}
