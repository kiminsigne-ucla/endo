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
#include "analysis.h"
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include <algorithm>
#include <R.h>
#include "string.h"

using namespace std;

analysis::analysis(){
	chromCounter = 0;
}

analysis::~analysis(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int analysis::processCoords(const char* inputFile,const char* outputFile,const char* twoBitFile,const char* chrmSearch, int binary){

	FILE * tempTB;
	time_t rtime;
	struct tm *timeinfo;
	
	char tInfo[128];// = "tempInfo.txt";
	char tChrSize[128];// = "tempChromSize.txt";
	char sysCall[256];
	time(&rtime);
	timeinfo=localtime(&rtime);
//	strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
	strftime(tChrSize,128,"tempChromSize_%H_%M_%S.txt",timeinfo);
	
	char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
	twoBitInfo2(twoBitFile2, tChrSize);
	free(twoBitFile2);

	/*tempTB = fopen(tInfo,"w");
	fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
	fclose (tempTB);
	
	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
	int s = 1;
	sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
	while(s != 0){
		s = system(sysCall);
		if(s != 0)
			cout << "Trying twoBitInfo again, s is" << s << endl;
	}
	remove(tInfo);
	*/

	tempTB = fopen(tChrSize,"r");
	char tbChrom[128];
	unsigned long int tbStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",tbChrom,&tbStart);
		if(ret == 2){
			cout << "\tFor " << tbChrom << " length is " << tbStart << endl;
			string sChr(tbChrom);
			unsigned short int chromIntTB = getHashValue(sChr.c_str());
			chr_size[chromIntTB] = tbStart;
		}
	}
	fclose(tempTB);
	remove(tChrSize);
	
	unsigned short int chromInt;
	string line;
	string field;
//////////////////////////////
	int profile_extend = 500;
//	int profile_extend = 0;
//////////////////////////////
	int printFlag = 0;
	int getChrmData = 0;
	unsigned short int collectData = 0;
	list<coord2>::iterator i = coordOUT_slist.begin();
	unsigned int * basepair = NULL;
	unsigned int * profile = NULL;
	unsigned long int countBases = 250000000;
	unsigned long int startOffset = 0;
	cout << "\nGetting basecount data from " << inputFile << endl;
	ifstream seqfile(inputFile);
	if(!seqfile.is_open()){
		cout << "ERROR opening input file" << inputFile << ", exiting" << endl;
		return 1;
	}

	while (getline(seqfile, line)){
		if (line[0] == 'f'){
			if(collectData == 1 && getChrmData == 1){
				i = coordOUT_slist.begin();
				while(i!=coordOUT_slist.end()){
					if(i->chrom == chromInt){
						int pIndex = 0;
						unsigned long int startPos = 1;
						if(i->start > profile_extend)
							startPos = i->start-profile_extend;
						unsigned long int stopPos = chr_size[chromInt];
						if((i->end+profile_extend) < chr_size[chromInt])
							stopPos = i->end+profile_extend;
						profile = new unsigned int[(stopPos-startPos)+1];
						for(int s = startPos; s <= stopPos; s++){
							profile[pIndex] = basepair[s];
							pIndex++;
						}
						if(outputData( outputFile,printFlag,i->chrom,startPos,stopPos,i->sigVal,i->qVal,pIndex,profile) != 0){
							cout << "Error printing output to file, exiting" << endl;
							return 1;
						}
						delete [] profile;
						profile = NULL;
						printFlag = 1;
						coordOUT_slist.erase(i++);
					}else{
						i++;
					}
				}
				delete [] basepair;
				basepair = NULL;

				if(coordOUT_slist.empty()){
					seqfile.close();
					return 0;
				}
			}
			
			cout << "\tProcessing: " << line << endl;
			istringstream iss(line);
			while(iss >> field){
				if (field[0] == 'c'){
					string chrom;
					for ( int it= 6; it < field.length(); it++ ){
						chrom += field.at(it);
					}
					chromInt = getHashValue(chrom.c_str());					
					string allChrm = "all";
					if(strcmp(chrmSearch,allChrm.c_str()) == 0){
						getChrmData = 1;
					}else if (strcmp(chrmSearch,chrom.c_str()) == 0){
						getChrmData = 1;
					}else{
						getChrmData = 0;
						cout << "\t\tSkipping " << chrom.c_str() << endl;
					}
					
					if(getChrmData == 1){
						cout << "\t\tLoading data for " << chrom.c_str() << endl;
						cout << "\t\tInitializing length of " << chr_size[chromInt] << endl;
						basepair = new unsigned int[chr_size[chromInt]+1];
						for(int c = chr_size[chromInt];c--;)
							basepair[c] = 0;
						countBases = 0;
					}
				}else if (field[0] == 's' && field[2] == 'a' && getChrmData == 1){
					string startVal;
					for ( int it= 6; it < field.length(); it++ ){
						startVal += field.at(it);
					}
					startOffset = atoi(startVal.c_str());
					cout << "\t\tStarting position in chromosome is " << startOffset << endl;
					countBases = startOffset - 1;
				}
			}
		}else if(getChrmData == 1 && line[0] != 't'){
				collectData = 1;
				countBases++;
				if(countBases > chr_size[chromInt]){
					cout << "\nERROR: Current position is " << countBases << " and chromosome length is " << chr_size[chromInt] << endl;
					return 1;
				}
				basepair[countBases] = atoi(line.c_str());
		}
	}

	if(getChrmData == 1){
		i = coordOUT_slist.begin();
		while(i!=coordOUT_slist.end()){
			if(i->chrom == chromInt){
				int pIndex = 0;
				unsigned long int startPos = 1;
				if(i->start > profile_extend)
					startPos = i->start-profile_extend;
				unsigned long int stopPos = i->end+profile_extend;
				profile = new unsigned int[(stopPos-startPos)+1];
				for(int s = startPos; s <= stopPos; s++){
					profile[pIndex] = basepair[s];
					pIndex++;
				}
				if(outputData( outputFile,printFlag,i->chrom,startPos,stopPos,i->sigVal,i->qVal,pIndex,profile) != 0){
					cout << "Error printing output to file, exiting" << endl;
					return 1;
				}
				delete [] profile;
				profile = NULL;
				printFlag = 1;
				coordOUT_slist.erase(i++);
			}else{
				i++;
			}
		}
	}
	seqfile.close();
	delete [] basepair;
	basepair = NULL;
	return 0;
}

int analysis::outputData(const char * outputFile,int pFlag,unsigned short int pChrom,unsigned long int pStart,unsigned long int pStop,double pSigVal, double pqVal,int printStop,unsigned int pProfile[]){
	FILE * fh;
	if(pFlag == 0){
		fh = fopen(outputFile,"w");
		if(fh==NULL){
			cout << "Unable to open output file: " << outputFile << endl;
			return 1;
		}
		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND\tSIGVAL\tQVAL");
		for(int p = 1;p <= printStop;p++){
			fprintf(fh,"\tPosition%i",p);
		}
		fprintf(fh,"\n");
	}else if (pFlag == 1){
		fh = fopen(outputFile,"a");
		if(fh==NULL){
			cout << "Unable to open output file: " << outputFile << endl;
			return 1;
		}
	}
	const char * chromName = getKey(pChrom);
	char strand[] = "+";
	char winID[255];
	char start[10];
	char stop[10];
	sprintf( start,"%lu", pStart);
	sprintf( stop,"%lu", pStop);
	strcpy(winID,chromName);strcat(winID,":");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
	fprintf(fh,"%s\t%s\t%lu\t%lu\t%s\t%.14f\t%.14f",winID,chromName,pStart,pStop,strand,pSigVal,pqVal);
	for(int posP = 0; posP < printStop;posP++){
		fprintf(fh,"\t%i",pProfile[posP]);
	}
	fprintf(fh,"\n");
	fclose (fh);
	return 0;
}

unsigned short int analysis::getHashValue(const char *currChrom){
	map<const char*, int>::iterator i;
	i = chroms.find(currChrom);
	if(i == chroms.end()){
		char * chromosome = new char[128];
		strcpy(chromosome,currChrom);
		chroms[chromosome] = chromCounter;
		intsToChrom[chromCounter] = chromosome;
		return(chromCounter++);
	}else{
		return i->second;	
	}
}

const char * analysis::getKey(unsigned short int chrom){
	map<int, const char*>::iterator i;
	i = intsToChrom.find(chrom);
	if(i == intsToChrom.end()){
		cout << chrom << endl;
		cout << "REALLY REALLY BAD ERROR!" << endl;
		exit(1);
	}else{
		return i->second;	
	}
}

int analysis::importCoords(const char *winlist,double threshold,const char *method,int wformat, int winGap, int FDR){
	
	const char *pscl = "pscl";
	const char *mixture = "mixture";
	int rval;
	if(strcmp(method,pscl) == 0){
		cout << "Getting significant regions from pscl run" << endl;
		rval = importPscl(winlist,threshold,wformat);
	}else if(strcmp(method,mixture) == 0){
		cout << "Getting significant regions from mixture run" << endl;
		rval = importMixture(winlist,threshold,wformat,FDR);
	}else{
		cout << "ERROR: method must be either pscl OR mixture" << endl;
		return 1;
	}

	if(rval == 0){
		cout << "\nImported " << coordIN_slist.size() << " coordinates" << endl;
		coordIN_slist.sort();
		list<coord2>::iterator back =  coordIN_slist.begin();
		///////////////////////////
		unsigned long int winSizeThresh = (back->end - back->start) * 10;
		//////////////////////////
		coord2 tempCoord = *back;
		coordIN_slist.erase(back++);
		int flagFinish = 0;
		
		while(flagFinish == 0){
			if(coordIN_slist.empty())
				flagFinish = 1;
			if(flagFinish == 0 && back->chrom == tempCoord.chrom && back->start <= tempCoord.end+1+winGap){
				tempCoord.end = back->end;
				if(strcmp(method,pscl) == 0 && tempCoord.sigVal > back->sigVal){
					tempCoord.sigVal = back->sigVal;
				}else if(strcmp(method,mixture) == 0 && tempCoord.sigVal < back->sigVal & FDR==0){
					tempCoord.sigVal = back->sigVal;
					tempCoord.qVal = back->qVal;
				}else if(strcmp(method,mixture) == 0 && tempCoord.qVal > back->qVal & FDR==1){
					tempCoord.sigVal = back->sigVal;
					tempCoord.qVal = back->qVal;
				}

				coordIN_slist.erase(back++);
			}else{
//				if((tempCoord.end-tempCoord.start) <= winSizeThresh){
					coordOUT_slist.push_back(tempCoord);
//				}else{
					//REMOVE ONCE C PEAKBOUND IS RUNNING
//					const char * cName = getKey(tempCoord.chrom);
//					cout << "Excluding " << cName << ":" << tempCoord.start << "-" << tempCoord.end << " SIZE=" << (tempCoord.end-tempCoord.start) << endl;
//				}

				if(flagFinish == 0){
					tempCoord = *back;
					coordIN_slist.erase(back++);
				}
			}
		}
		coordOUT_slist.sort();
		cout << "\nCollapsed to " << coordOUT_slist.size() << " non-overlapping regions" << endl;
		return 0;
	}else{
		cout << "ERROR: problem importing windows" << endl;
		return 1;
	}
}

int analysis::importPscl(const char *winlist,double threshold,int wformat){

	FILE * wlist;
	wlist = fopen(winlist,"r");
	if(wlist == NULL){
		cout << "ERROR: opening window list file" << winlist << endl;
		return 1;
	}
	char cChrom[128];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	double sigVal;
	char sigFile [256];
	int readResult = 0;
	char firstline [256];
	int rwline;
	double sres;int ec;double ic;double gc;double ap;double cnv;

	while(!feof(wlist)){
		int rline = fscanf(wlist,"%s",sigFile);
		if(rline == 1){
			string signalFile(sigFile);
			FILE * fh;
			fh = fopen(signalFile.c_str(),"r");
			if(fh != NULL){
				cout << "\tImporting windows from " << signalFile.c_str() << "..." << endl;
				fgets(firstline,256,fh);
				while(!feof(fh)){
					readResult = 0;
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%lf%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal,&sres);
						if(sigVal <= threshold && rwline == 6){
							readResult = 1;
						}
					}else if(wformat == 1){
						rwline = fscanf(fh,"%s%lu%lu%d%lf%lf%lf%lf%hu%lf%lf",cChrom,&iStart,&iEnd,&ec,&ic,&gc,&ap,&cnv,&qFlag,&sigVal,&sres);
						if(sigVal <= threshold && rwline == 11){
							readResult = 1;
						}
					}

					if(readResult == 1){
						if(qFlag == 1){
							string chromIn(cChrom);
							unsigned short int chromInt = getHashValue(chromIn.c_str());
							coord2 c(chromInt,iStart,iEnd,qFlag,sigVal,sigVal);
							coordIN_slist.push_back(c);
						}
					}
				}
				fclose(fh);
			}else{
				cout << "ERROR: unable to open .wins file " << signalFile.c_str() << endl;
				return 1;
			}
		}else{
			cout << "\tSkipping line, " << sigFile << endl;
		}
	}
	fclose(wlist);
	return 0;
}

int analysis::importMixture(const char *winlist,double threshold,int wformat, int FDR){

	FILE * wlist;
	wlist = fopen(winlist,"r");
	if(wlist == NULL){
		cout << "ERROR: opening window list file" << winlist << endl;
		return 1;
	}
	char cChrom[128];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	double sigVal;	
	double qVal;
	char sigFile [256];
	int readResult = 0;
	char firstline [256];
	int rwline;
	int ec;double ic;double gc;double ap;double cnv;
	
	while(!feof(wlist)){
		int rline = fscanf(wlist,"%s",sigFile);
		if(rline == 1){
			string signalFile(sigFile);
			FILE * fh;
			fh = fopen(signalFile.c_str(),"r");
			if(fh != NULL){
				cout << "\tImporting windows from " << signalFile.c_str() << "..." << endl;
				fgets(firstline,256,fh);
				while(!feof(fh)){
					readResult = 0;
					if(wformat == 0 & FDR==0){
						rwline = fscanf(fh,"%s%lu%lu%hu%lf%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal,&qVal);
						if(sigVal >= threshold && rwline == 6){
							readResult = 1;
						}
					}else if(wformat == 1 & FDR==0){
						rwline = fscanf(fh,"%s%lu%lu%d%lf%lf%lf%lf%hu%lf%lf",cChrom,&iStart,&iEnd,&ec,&ic,&gc,&ap,&cnv,&qFlag,&sigVal,&qVal);
						if(sigVal >= threshold && rwline == 11){
							readResult = 1;
						}
					}else if(wformat == 0 & FDR==1){
						rwline = fscanf(fh,"%s%lu%lu%hu%lf%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal,&qVal);
						if(qVal <= threshold && rwline == 6){
							readResult = 1;
						}
					}else if(wformat == 1 & FDR==1){
						rwline = fscanf(fh,"%s%lu%lu%d%lf%lf%lf%lf%hu%lf%lf",cChrom,&iStart,&iEnd,&ec,&ic,&gc,&ap,&cnv,&qFlag,&sigVal,&qVal);
						if(qVal <= threshold && rwline == 11){
							readResult = 1;
						}
					}
					
					//Rprintf("chrom %s start %lu stop %lu qflag %d pp %lf qval %lf thresh %f readres %d FDR %d\n", cChrom,iStart,iEnd,qFlag,sigVal,qVal, threshold, readResult, FDR);
					if(readResult == 1){
						if(qFlag == 1){
							string chromIn(cChrom);
							unsigned short int chromInt = getHashValue(chromIn.c_str());
							coord2 c(chromInt,iStart,iEnd,qFlag,sigVal,qVal);
							coordIN_slist.push_back(c);
						}
					}
				}
				fclose(fh);
			}else{
				cout << "ERROR: unable to open .wins file " << signalFile.c_str() << endl;
				return 1;
			}
		}else{
			cout << "\tSkipping line, " << sigFile << endl;
		}
	}
	fclose(wlist);
	return 0;
}

