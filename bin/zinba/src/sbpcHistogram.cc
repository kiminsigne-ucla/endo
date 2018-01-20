#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "histsbpc.h"
#include <R.h>
#include <Rmath.h>
#define MAX_LEN 1025

using namespace std;

extern "C" {
void sbpchistogram(char **Rinput, char **Routput,char **Rtwobitfile,int *Rmax, int *Rmin, int *RnBin){
	const char * input = Rinput[0];
	const char * output = Routput[0];
	const char * twoBit = Rtwobitfile[0];
	int max = Rmax[0];
	int min = Rmin[0];
	int numBin = RnBin[0];

	histsbpc newHist;
	int ret = newHist.hist_data(input,output,twoBit,max,min,numBin);
	if(ret != 0)
		Rprintf("\n---------------- ALIGN ADJUST EXITING WITH ERRORS ----------------\n");
	else
		Rprintf("\n---------------- ALIGN ADJUST COMPLETED SUCCESSFULLY ----------------\n");
}
}
