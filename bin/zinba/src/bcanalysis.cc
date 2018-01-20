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
#include "bcanalysis.h"
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <R.h>
#include "import.h"
#include "import.h"
using namespace std;

/*
bcanalysis::bcanalysis(){
	chromCounter = 0;
	tbSizeFlag = 0;
}

bcanalysis::~bcanalysis(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}
*/

int bcanalysis::processSignals(const char* outputFile,int extend, int binary){
	//import b;
	unsigned short int currchr;
	int i;
	const char * chromReport;
//	unsigned short int * basepair = NULL;
	unsigned int * basepair = NULL;
	unsigned short int printFLAG = 0;
	unsigned long int aStart,aStop;
	cout << "Processing " << b.signal_slist.size() << " reads" << endl;
	
	while(!b.signal_slist.empty()){
		i = 0;
		currchr = b.signal_slist[0].chrom;
		chromReport = b.getKey(currchr);
		cout << "\tProcessing " << chromReport << ".........." << endl;
		cout << "\t\tInitializing length to " << b.chr_size[currchr] << endl;
//		basepair = new unsigned short int[b.chr_size[currchr]+1];
		basepair = new unsigned int[b.chr_size[currchr]+1];
		for(int in = b.chr_size[currchr]; in--;)
			basepair[in] = 0;

		while(b.signal_slist[i].chrom==currchr && i < (int) b.signal_slist.size()){
			aStart = b.signal_slist[i].pos;
			aStop = b.chr_size[currchr];
			if((aStart+extend-1) <= b.chr_size[currchr])
				aStop = aStart + extend - 1;
			for(unsigned int rpos = aStart; rpos <= aStop;rpos++)
				basepair[rpos]++;
			i++;
		}
		cout << "\t\t" << i << " reads mapped to " << chromReport << endl;

		if(outputData(outputFile,currchr,printFLAG,basepair,binary) != 0){
			return 1;
		}else{
			cout << "\t\tPrinted results to " << outputFile << endl;
		}
		printFLAG = 1;
		delete [] basepair;
		basepair = NULL;
		b.signal_slist.erase(b.signal_slist.begin(),b.signal_slist.begin()+i);
	}
	return 0;
}

//int bcanalysis::outputData(const char * outputFile, unsigned short int currChr,unsigned short int pFLAG,unsigned short int basepair[]){
int bcanalysis::outputData(const char * outputFile, unsigned short int currChr,unsigned short int pFLAG,unsigned int basepair[], int binary){
	//import b;
	if(binary ==0){
		FILE * fh;
		const char * chrom = b.getKey(currChr);
		if(pFLAG == 0){
			fh = fopen(outputFile,"w");
			fprintf(fh,"track type=wiggle_0 name=\"%s\" desc=\"%s\" visibility=full\n", 	
							outputFile,outputFile);
		}else{
			fh = fopen(outputFile,"a");
		}
		if(fh == NULL){
			cout << "\nERROR: Can't open output file: " << outputFile << endl;
			return 1;
		}
		fprintf(fh,"fixedStep chrom=%s start=1 step=1\n",chrom);
	
		for(int posP = 1; posP <= b.chr_size[currChr];posP++)
			fprintf(fh,"%u\n",basepair[posP]);
			//fprintf(fh,"%hu\n",basepair[posP]);
		fclose(fh);

	}else{
		const char * chrom = b.getKey(currChr);
		ofstream ofile;
		string outputFilestring = string(outputFile)+"_"+ string(chrom);		
		ofile.open (outputFile,ios::out|ios::binary);
		if(!ofile.is_open()){
			cout << "ERROR: Can't open output file: " << outputFile << endl;
			return 1;
		}
		for(int posP = 1; posP <= b.chr_size[currChr];posP++)
			ofile.write((const char*) &basepair[posP], sizeof(unsigned int));
		ofile.close();
	}
		return 0;
}

/*
unsigned short int bcanalysis::getHashValue(char *currChrom){
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

const char * bcanalysis::b.getKeyunsigned short int chrom){
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
*/
int bcanalysis::importRawSignal(const char * signalFile,int extension,const char * filetype,const char * twoBitFile){
	//import b;
	if(b.tbSizeFlag == 0){
		FILE * tempTB;
		time_t rtime;
		struct tm *timeinfo;
		
		char tInfo[128];// = "tempInfo.txt";
		char tChrSize[128];// = "tempChromSize.txt";
		char sysCall[256];
		time(&rtime);
		timeinfo=localtime(&rtime);
		//strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tChrSize,128,"tempChromSize_%H_%M_%S.txt",timeinfo);

		char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
		twoBitInfo2(twoBitFile2, tChrSize);
		free(twoBitFile2);		
		/*tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
		fclose (tempTB);
		
		cout << "\tGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
		int s = 1;
		int twobitCount = 0;
		sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitInfo again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "TwoBitInfo failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);
		*/

		tempTB = fopen(tChrSize,"r");
		char cChrom[128];
		unsigned long int cStart;
		while(!feof(tempTB)){
			int ret = fscanf(tempTB,"%s%lu",cChrom,&cStart);
			cout << "\t\tFor " << cChrom << " length is " << cStart << endl;
			unsigned short int chromInt = b.getHashValue(cChrom);
			b.chr_size[chromInt] = cStart;
		}
		fclose(tempTB);
		remove(tChrSize);
		b.tbSizeFlag = 1;
	}
	
	const char * bowtie = "bowtie";
	const char * bed = "bed";
	const char * tagAlign = "tagAlign";
	int rval = 0;
	
	if(strcmp(filetype,bed) == 0){
		rval = b.importBed(signalFile,extension, 0);
	}else if (strcmp(filetype,bowtie) == 0){
		rval = b.importBowtie(signalFile,extension, 0);
	}else if(strcmp(filetype,tagAlign) == 0){
		rval = b.importTagAlign(signalFile,extension, 0);
	}else{
		cout << "Unrecognized type of file " << filetype << ", must be either bowtie, bed, or tagAlign" << endl;
		return 1;
	}
	
	if(rval == 0){
		cout << "\tLoaded " << b.signal_slist.size() << " reads" << endl;
		cout << "\tSorting reads ...";
		sort (b.signal_slist.begin(), b.signal_slist.end());
		cout << "COMPLETE" << endl;
	}else{
		cout << "Stopping basealigncounts" << endl;
		return 1;
	}
	return 0;
}

/*
int bcanalysis::importBowtie(const char * signalFile,int extension){
	
	cout << "\nLoading bowtie formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "\nERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	char line[512];char seq[128];
	char name[128];char sscore[128];int ival;
	int rval;
	int num_skip = -1;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%s%s%lu%s%s%i",name,strand,cChrom,&pos,seq,sscore,&ival);
		fgets(line,512,fh);
		if(rval == 7){
			if(strcmp(strand,minus) == 0){
				if((pos + strlen(seq)) >= extension)
					pos = (pos + strlen(seq)) - extension + 1;
				else
					pos = 1;
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			b.signal_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}
	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}

int bcanalysis::importTagAlign(const char * signalFile,int extension){
	
	cout << "\nLoading tagAlign formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char seq[128];int score;
	int rval;
	int num_skip = -1;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,seq,&score,strand);
		if(rval == 6){
			if(strcmp(strand,minus) == 0){
				if(stop >= extension)
					pos = stop - extension + 1;
				else
					pos = 1;
			}else{
				pos = start;
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			b.signal_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}
	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}

int bcanalysis::importBed(const char * signalFile,int extension){
	
	cout << "\nLoading bed formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	unsigned short int sval = 1;
	char strand[1];
	char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char name[128];int bscore;
	int rval;
	int num_skip = -1;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,name,&bscore,strand);
		if(rval == 6){
			if(strcmp(strand,minus) == 0){
				if(stop >= extension)
					pos = stop - extension + 1;
				else
					pos = 1;
			}else{
				pos = start;
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			b.signal_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}
	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}
*/
