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
void signalnoise(char **Rinput, char **Routput,char **Rtwobitfile, int *RnBin){
	const char * input = Rinput[0];
	const char * output = Routput[0];
	const char * twoBit = Rtwobitfile[0];
	int numBin = RnBin[0];

	histsbpc newHist;
	int ret = newHist.signalnoise(input,output,twoBit,numBin);
	if(ret != 0)
		Rprintf("\n---------------- SIGNAL TO NOISE ANALYSIS EXITING WITH ERRORS ----------------\n");
	else
		Rprintf("\n---------------- SIGNAL TO NOISE ANALYSIS COMPLETED SUCCESSFULLY ----------------\n");
}
}
