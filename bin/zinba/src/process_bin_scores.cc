#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short	      /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int	      /* Wants to be signed 32 bits. */
#define boolean bool	      /* Wants to be signed 32 bits. */
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "process_bin_scores.h"
#include <sstream>
#include <string>
#include <string.h>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include "R.h"
extern "C"{
#include "twoBit.h"
}
using namespace std;

process_bin_scores::process_bin_scores(){}

process_bin_scores::~process_bin_scores(){}

int process_bin_scores::adjustCoords(const char * filelist,string outDir,const char* twoBitFile,int aThresh,int adjustSize){
	
	FILE * tempTB;
	struct tm *timeinfo;
	time_t rtime;
	//const char * tChrSize = "tempChromSize.txt";
	//cout << "\nGetting chrm size from " << twoBitFile << endl;
	/*tempTB = fopen(tInfo,"w");
	fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
	fclose (tempTB);
	
	//	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
	int s = system("R CMD BATCH tempInfo.txt /dev/null");
	if(s != 0){
		cout << "\nERROR: failed to generate info file using twoBitInfo, twoBit file is: " << twoBitFile << endl;
		return 1;
	}
	remove(tInfo);	
	*/
	
	char tChrSize[128];
	char tChrSize2[512];
	char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
	time(&rtime);
	timeinfo=localtime(&rtime);
	strftime(tChrSize,128,"tempChromSize_%H_%M_%S",timeinfo);
		struct timeval time;
	     	gettimeofday(&time,NULL);
	     	srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
		int random = rand();
		char rand[1024];
		sprintf(rand,"_%d.txt", random);
	//Rprintf("%s %d %s", filelist, random, rand);
	strcat(tChrSize2, tChrSize);
	//Rprintf(" %s %s ", tChrSize, tChrSize2);
	strcat(tChrSize2, rand);
	//Rprintf(" %s\n", tChrSize2);

	twoBitInfo2(twoBitFile2, tChrSize2);
	free(twoBitFile2);

	tempTB = fopen(tChrSize2,"r");
	char tbChrom[128];
	unsigned long int tbStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",tbChrom,&tbStart);
		if(ret == 2){
		//	cout << "\tFor " << tbChrom << " length is " << tbStart << endl;
			string sChr(tbChrom);
			chr_size[sChr] = tbStart;
		}
	}
	fclose(tempTB);
	remove(tChrSize2);
	
	//cout << "\nProcessing chromosome files......" << endl;
	FILE * flist;
	flist = fopen(filelist,"r");
	if(flist == NULL){
		cout << "ERROR: opening file list file" << filelist << endl;
		return 1;
	}
	
	char mapfile[128];
	char mapChar[1];
	unsigned int map_count;
	unsigned long int countBases = 0;
	unsigned short int * align_count = NULL;
	int collectdata = 0;
	unsigned long int startOffset = 0;
	string chrm;
	
	while(!feof(flist)){
		fscanf(flist,"%s",mapfile);
		//cout << "\tProcessing " << mapfile << endl;
		chrm = string(mapfile);
		size_t bout = chrm.find("b.out");
		chrm.erase(bout);
	        size_t  last = chrm.find_last_of('/');
		if(chrm.npos != last)  chrm = chrm.substr(last+1);

		//cout << "Chromosome is " << chrm << endl;
		align_count = new unsigned short int[chr_size[chrm]+1];
		for(int c = chr_size[chrm];c--;)
			align_count[c] = 0;
		
		ifstream mfile (mapfile,ios::in|ios::binary);
		countBases = 0;
		if(mfile.is_open()){
			while(!mfile.eof()){
				mfile.read(mapChar,1);
				unsigned char mchr = mapChar[0];
				countBases++;
				if((int) mchr <= aThresh && (int) mchr > 0){
					align_count[countBases] = 1;
				}
			}
		}else{
			cout << "ERROR: Unable to open chromosome file " << mapfile << endl;
			return 1;
		}
		mfile.close();
		
		string outputFile = outDir + chrm + ".bwig";
		ofstream ofile (outputFile.c_str(),ios::out|ios::binary);
		
		if(!ofile.is_open()){
			cout << "ERROR: Can't open output file: " << outputFile << endl;
			return 1;
		}
		
		int alignVal = 0;
		for(int pos = 1; pos < adjustSize;pos++)
			ofile.write((const char*)&alignVal,1);
		
		for(int i = adjustSize;i <= (countBases - adjustSize); i++){
			alignVal = 0;
			alignVal += align_count[(i-adjustSize+1)];
			alignVal += align_count[(i+adjustSize)];
			ofile.write((const char*)&alignVal,1);
		}

		alignVal = 0;
		for(int pos = (countBases - adjustSize+1); pos <= countBases;pos++)
			ofile.write((const char*)&alignVal,1);
		
		if(countBases < chr_size[chrm])
			for(int pos = countBases+1;pos <= chr_size[chrm];pos++)
				ofile.write((const char*)&alignVal,1);

		ofile.close();
		delete [] align_count;
		align_count = NULL;
	}
	fclose(flist);
	return 0;
}
