#include <iostream>
#include <fstream>
#include <stdio.h>
#include "winCounts.h"
//#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <Rmath.h>

#include <algorithm>
#include <time.h>

using namespace std;

winCounts::winCounts(){
	chromCounter = 0;
	tbSizeFlag = 0;
}

winCounts::~winCounts(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int winCounts::processSignalsWC(const char * twoBitFile,const char * outfile){

	FILE * tempTB;
	time_t rtime;
	struct tm *timeinfo;
	char tInfo[128];// = "tempInfo.txt";
	char sysCall[256];
	
	unsigned short int currchr = 999;
	unsigned int * basepair = NULL;
	int i;
	int wcount = 0;
	int printflag = 0;
	
	while(!signal_slist.empty()){
		i = 0;
		currchr = signal_slist[0].chrom;
		const char * chromReport = getKey(currchr);
		cout << "\nProcessing " << chromReport << endl;
		cout << "\tLength is " << chr_size[currchr] << endl;
		
		basepair = new unsigned int[chr_size[currchr]+1];
		for(int ch = chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(signal_slist[i].chrom==currchr && i < (int) signal_slist.size()){
			basepair[signal_slist[i].pos]++;
			i++;
		}
		signal_slist.erase(signal_slist.begin(),signal_slist.begin()+i);

		cout << "\tGetting counts for windows.........." << endl;
		list<coord>::iterator c = coord_slist.begin();
		while(c!=coord_slist.end()){
			wcount = 0;
			if(c->chrom == currchr){
				for(int s = c->start; s <= c->end; s++){
					wcount += basepair[s];
				}
				c->sigVal = (double) wcount;
			}
			c++;
		}
		if(outputDataWC(outfile,currchr,printflag) != 0){
			cout << "Error printing output to file, exiting" << endl;
			exit(1);
		}
		printflag = 1;
		delete [] basepair;
		basepair = NULL;
	}
	return 0;
}

int winCounts::outputDataWC(const char * outputFile, unsigned short int currChr,int pflag){
	FILE * fh;
	if(pflag == 0){
		fh = fopen(outputFile,"w");
		fprintf(fh,"REGID\tCHROM\tSTART\tSTOP\tSTRAND\tSCORE\tWIN_COUNT\n");
	}else{
		fh = fopen(outputFile,"a");
	}
	if(fh == NULL){
		cout << "\nERROR: Can't open output file: " << outputFile << endl;
		return 1;
	}
	const char * chrom = getKey(currChr);
	char winID[255];
	char start[10];
	char stop[10];
	list<coord>::iterator c = coord_slist.begin();
	while(c != coord_slist.end()){
		if(c->chrom == currChr){
			sprintf( start,"%lu", c->start);
			sprintf( stop,"%lu", c->end);
			strcpy(winID,chrom);strcat(winID,"-");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
			
			fprintf(fh,"%s\t%s\t%lu\t%lu\t+\t0.05\t%i\n",winID,chrom,c->start,c->end,(int) c->sigVal);
			coord_slist.erase(c++);
		}else{
			c++;
		}
	}
	fclose (fh);
	return 0;
}

unsigned short int winCounts::getHashValue(char *currChrom){
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

const char * winCounts::getKey(unsigned short int chrom){
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

int winCounts::importRawSignalWC(const char * signalFile,int extension,const char * twoBitFile){

	if(tbSizeFlag == 0){
		FILE * tempTB;
		time_t rtime;
		struct tm *timeinfo;

		char tInfo[128];// = "tempInfo.txt";
		char tChrSize[128];// = "tempChromSize.txt";
		char sysCall[256];
		time(&rtime);
		timeinfo=localtime(&rtime);
		strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tChrSize,128,"tempChromSize_%H_%M_%S.txt",timeinfo);

		tempTB = fopen(tInfo,"w");
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

		tempTB = fopen(tChrSize,"r");
		char cChrom[128];
		unsigned long int cStart;
		while(!feof(tempTB)){
			int ret = fscanf(tempTB,"%s%lu",cChrom,&cStart);
			cout << "\tFor " << cChrom << " length is " << cStart << endl;
			unsigned short int chromInt = getHashValue(cChrom);
			chr_size[chromInt] = cStart;
		}
		fclose(tempTB);
		remove(tChrSize);
		tbSizeFlag = 1;
	}

	cout << "\nImporting reads from " << signalFile << endl;
	int rval = importBowtieWC(signalFile,extension);
	if(rval == 0){
		cout << "\tImported " << signal_slist.size() << " reads" << endl;
		cout << "\tSorting reads ...";
		sort (signal_slist.begin(), signal_slist.end());
		cout << "COMPLETE" << endl;
	}else{
		return 1;
	}
	return 0;
}

int winCounts::importBowtieWC(const char * signalFile,int extension){
	
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
	char line[512];char seq[128];
	char name[128];char sscore[128];int ival;
	int extend = (int)(extension/2);
	int rval;
	int num_skip = -1;

	while(!feof(fh)){
		rval = fscanf(fh,"%s%s%s%lu%s%s%i",name,strand,cChrom,&pos,seq,sscore,&ival);
		fgets(line,512,fh);
		if(rval == 7){
			if(strcmp(strand,minus) == 0){
				if((pos + strlen(seq)) >= extend)
					pos = (pos + strlen(seq)) - extend + 1;
				else
					pos = 1;
			}else{
				if((pos + extend-1) <= chr_size[getHashValue(cChrom)])
					pos += extend - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			signal_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}
	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}
	
int winCounts::importCoordsWC(const char *coordfile){
	char cChrom[128];
	char cid[128];
	char strand[1];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	double score;
	char firstline [256];
	int rwline;
	
	FILE * cfile;
	cout << "\tOpening coord file " << coordfile << endl;
	cfile = fopen(coordfile,"r");
	if(cfile == NULL){
		cout << "ERROR: opening window list file" << endl;
		return 1;
	}
	
	fgets(firstline,256,cfile);
	while(!feof(cfile)){
		rwline = fscanf(cfile,"%s%s%lu%lu%s%lf",cid,cChrom,&iStart,&iEnd,strand,&score);
		if(rwline == 6){
			unsigned short int chromInt = getHashValue(cChrom);
//			iStart += 500;iEnd -= 500;
			coord c(chromInt,iStart,iEnd,1,score);
			coord_slist.push_back(c);
		}
	}
	fclose(cfile);
	
	cout << "\nImported " << coord_slist.size() << " coordinates.....\nSorting coords......." << endl;
	coord_slist.sort();
	cout << "COMPLETE" << endl;
	return 0;
}

