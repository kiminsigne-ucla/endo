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

#include <iostream>
#include <fstream>
#include <stdio.h>
extern "C"{
#include "twoBit.h"
}

#include "process_scores.h"
#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "string.h"
using namespace std;

process_scores::process_scores(){}

process_scores::~process_scores(){}

int process_scores::adjustCoords(string alignFile,string outDir,const char* twoBitFile,int aThresh,int adjustSize){
	
	FILE * tempTB;
	//const char * tInfo = "tempInfo.txt"; 
	char * tChrSize = "tempChromSize.txt";
	/*tempTB = fopen(tInfo,"w");
	fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
	fclose (tempTB);
	
	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
	int s = system("R CMD BATCH tempInfo.txt /dev/null");
	if(s != 0){
		cout << "\nERROR: failed to generate info file using twoBitInfo, twoBit file is: " << twoBitFile << endl;
		return 1;
	}
	remove(tInfo);	
	*/

		char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
		twoBitInfo2(twoBitFile2, tChrSize);
		free(twoBitFile2);

	tempTB = fopen(tChrSize,"r");
	char tbChrom[128];
	unsigned long int tbStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",tbChrom,&tbStart);
		if(ret == 2){
			string sChr(tbChrom);
			chr_size[sChr] = tbStart;
		}
	}
	fclose(tempTB);
	remove(tChrSize);
	
	unsigned int map_count;
	unsigned long int countBases = 0;
	unsigned short int * align_count = NULL;
	int collectdata = 0;
	unsigned long int startOffset = 0;
	string line;
	string field;
	string chrom;
	FILE * fh;

	cout << "Loading align count data from " << alignFile.c_str() << endl;
	ifstream seqfile(alignFile.c_str());
	while (getline(seqfile, line)){
		if (line[0] != 't' && line[0] != 'f'){
			countBases++;
			if(atoi(line.c_str()) <= aThresh && atoi(line.c_str()) > 0){
				align_count[countBases] = 1;
			}
			collectdata = 1;
		}else if (line[0] == 'f'){
			if(collectdata == 1){
				string outputFile = outDir + chrom + ".wig";
				cout << "\t\tPrinting output to: " << outputFile.c_str() << endl;
				fh = fopen(outputFile.c_str(),"w");
				if(fh == NULL){
					cout << "\nERROR: Can't open output file: " << outputFile.c_str() << endl;
					return 1;
				}
				int alignVal = 0;
				for(int pos = 1; pos < adjustSize;pos++)
					fprintf(fh,"%i\n",alignVal);
				for(int i = adjustSize;i <= (countBases - adjustSize); i++){
					alignVal = 0;
					alignVal += align_count[(i-adjustSize+1)];
					alignVal += align_count[(i+adjustSize)];
					fprintf(fh,"%i\n",alignVal);
				}
				alignVal = 0;
				for(int pos = (countBases - adjustSize+1); pos <= countBases;pos++)
					fprintf(fh,"%i\n",alignVal);
				if(countBases < chr_size[chrom])
					for(int pos = countBases+1;pos <= chr_size[chrom];pos++)
						fprintf(fh,"%i\n",alignVal);
				fclose (fh);
				delete [] align_count;
				align_count = NULL;
			}
			istringstream iss(line);
			while(iss >> field){
				if (field[0] == 'c'){
					chrom = "";
					for (int it= 6; it < field.length(); it++ )
						chrom += field.at(it);
					align_count = new unsigned short int[chr_size[chrom]+1];
					for(int c = chr_size[chrom]; c--;)
						align_count[c] = 0;
				}else if (field[0] == 's' && field[2] == 'a'){
					string startVal;
					for (int it= 6; it < field.length(); it++)
						startVal += field.at(it);
					startOffset = atoi(startVal.c_str());
					countBases = startOffset - 1;
				}
			}
			cout << "\tProcessing " << chrom.c_str() << endl;
		}
	}
	if(collectdata == 1){
		string outputFile = outDir + chrom + ".wig";
		cout << "\t\tPrinting output to: " << outputFile.c_str() << endl;
		fh = fopen(outputFile.c_str(),"w");
		if(fh == NULL){
			cout << "\nERROR: Can't open output file: " << outputFile.c_str() << endl;
			return 1;
		}
		int alignVal = 0;
		for(int pos = 1; pos < adjustSize;pos++)
			fprintf(fh,"%i\n",alignVal);
		for(int i = adjustSize;i <= (countBases - adjustSize); i++){
			alignVal = 0;
			alignVal += align_count[(i-adjustSize+1)];
			alignVal += align_count[(i+adjustSize)];
			fprintf(fh,"%i\n",alignVal);
		}
		alignVal = 0;
		for(int pos = (countBases - adjustSize+1); pos <= countBases;pos++)
			fprintf(fh,"%i\n",alignVal);
		if(countBases < chr_size[chrom])
			for(int pos = countBases+1;pos <= chr_size[chrom];pos++)
				fprintf(fh,"%i\n",alignVal);
		fclose (fh);
		delete [] align_count;
		align_count = NULL;
	}	
	seqfile.close();
	return 0;
}
