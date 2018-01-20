#include <iostream>
#include <fstream>
#include <stdio.h>
#include "histsbpc.h"
#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>

using namespace std;

histsbpc::histsbpc(){
	chromCounter = 0;
}

histsbpc::~histsbpc(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int histsbpc::hist_data(const char * sbpcFile,const char * outfile,const char* twoBitFile,int max,int min,int numBin){

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
	
	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
	int s = 1;
	sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
	while(s != 0){
		s = system(sysCall);
		if(s != 0)
			cout << "Trying twoBitInfo again, s is" << s << endl;
	}
	remove(tInfo);
	
	tempTB = fopen(tChrSize,"r");
	char tbChrom[128];
	unsigned long int tbStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",tbChrom,&tbStart);
		if(ret == 2){
			string sChr(tbChrom);
			unsigned short int chromIntTB = getHashValue(sChr.c_str());
			chr_size[chromIntTB] = tbStart;
		}
	}
	fclose(tempTB);
	remove(tChrSize);
	
	int dataRange = max - min;
	int binSize;
	if((dataRange % numBin) != 0)
		numBin++;

	if(dataRange > numBin)
		binSize = (int)dataRange/numBin;
	else
		binSize = 1;

	unsigned long int * sbpc_hist = NULL;
	sbpc_hist = new unsigned long int [numBin];
	for(int b = numBin;b--;)
		sbpc_hist[b] = 0;
	
	int ltMin = 0;
	int gtMax = 0;

	unsigned short int * sbpc_count = NULL;
	string line;string field;
	string chrom;int chromInt;
	int countBases = 0;
	cout << "Getting data for histogram from " << sbpcFile << endl;
	cout << "\n\nCHRM\tMEDIAN\t75thPERCENTILE\t90thPERCENTILE\t95thPERCENTILE\tMAX" << endl;
	ifstream seqfile(sbpcFile);
	while (getline(seqfile, line)){
		if (line[0] != 't' && line[0] != 'f'){
			countBases++;
			sbpc_count[countBases] = atoi(line.c_str());
		}else if (line[0] == 'f'){
			if(countBases > 0){
				sort(sbpc_count, sbpc_count+(chr_size[chromInt]+1));
				
				int nonzeroInd = 0;
				while(sbpc_count[nonzeroInd] == 0)
					nonzeroInd++;
				
				int sInd = 0;
				int bMin = min;
				int bMax = bMin + binSize;
				for(int b = 0; b < numBin; b++){
					if(b > 0){
						bMin = bMax + 1;
						bMax = bMin + binSize - 1;
					}
					while (sbpc_count[sInd] < min && sInd <= chr_size[chromInt]){
						ltMin++;
						sInd++;
					}
					while(sbpc_count[sInd] >= bMin && sbpc_count[sInd] <= bMax && sInd <= chr_size[chromInt]){
						sbpc_hist[b]++;
						sInd++;
					}
					while(sbpc_count[sInd] > max && sInd <= chr_size[chromInt]){
						gtMax++;
						sInd++;
					}
				}
				
				int nzlength = chr_size[chromInt] - nonzeroInd;
				int medInd = (int) nzlength/2;
				int sfInd = (int) nzlength * 0.75;
				int nInd = (int) nzlength * 0.9;
				int nfInd = (int) nzlength * 0.95;
				cout << chrom << "\t" << sbpc_count[nonzeroInd+medInd] << "\t" << sbpc_count[nonzeroInd+sfInd] << "\t" << sbpc_count[nonzeroInd+nInd] << "\t" << sbpc_count[nonzeroInd+nfInd] << "\t" << sbpc_count[chr_size[chromInt]] << endl;
				delete [] sbpc_count;
				sbpc_count = NULL;
			}
			
			istringstream iss(line);
			while(iss >> field){
				if (field[0] == 'c'){
					chrom = "";
					for (int it= 6; it < field.length(); it++ )
						chrom += field.at(it);
					chromInt = getHashValue(chrom.c_str());
					sbpc_count = new unsigned short int[chr_size[chromInt]+1];
					for(int c = chr_size[chromInt]; c--;)
						sbpc_count[c] = 0;
				}else if (field[0] == 's' && field[2] == 'a'){
					string startVal;
					for (int it= 6; it < field.length(); it++)
						startVal += field.at(it);
					int startOffset = atoi(startVal.c_str());
					countBases = startOffset - 1;
				}
			}
		}
	}
	seqfile.close();

	if(countBases > 0){
		sort(sbpc_count, sbpc_count+(chr_size[chromInt]+1));
		
		int nonzeroInd = 0;
		while(sbpc_count[nonzeroInd] == 0)
			nonzeroInd++;
		
		int sInd = 0;
		int bMin = min;
		int bMax = bMin + binSize;
		for(int b = 0; b < numBin; b++){
			if(b > 0){
				bMin = bMax + 1;
				bMax = bMin + binSize - 1;
			}
			while (sbpc_count[sInd] < min && sInd <= chr_size[chromInt]){
				ltMin++;
				sInd++;
			}
			while(sbpc_count[sInd] >= bMin && sbpc_count[sInd] <= bMax && sInd <= chr_size[chromInt]){
				sbpc_hist[b]++;
				sInd++;
			}
			while(sbpc_count[sInd] > max && sInd <= chr_size[chromInt]){
				gtMax++;
				sInd++;
			}
		}
		
		int nzlength = chr_size[chromInt] - nonzeroInd;
		int medInd = (int) nzlength/2;
		int sfInd = (int) nzlength * 0.75;
		int nInd = (int) nzlength * 0.9;
		int nfInd = (int) nzlength * 0.95;
		cout << chrom << "\t" << sbpc_count[nonzeroInd+medInd] << "\t" << sbpc_count[nonzeroInd+sfInd] << "\t" << sbpc_count[nonzeroInd+nInd] << "\t" << sbpc_count[nonzeroInd+nfInd] << "\t" << sbpc_count[chr_size[chromInt]] << endl;
		delete [] sbpc_count;
		sbpc_count = NULL;
	}
	
	
	cout << "\n\n" << ltMin << " bp less than " << min << endl;
	cout << gtMax << " bp greater than " << max << endl;
	
	FILE * fh;
	fh = fopen(outfile,"w");
	if(fh==NULL){
		cout << "Unable to open output file: " << outfile << endl;
		return 1;
	}
	fprintf(fh,"BIN_START\tBIN_STOP\tCOUNT\n");

	int bMin = min;
	int bMax = bMin + binSize;
	for(int b = 0; b < numBin; b++){
		if(b > 0){
			bMin = bMax + 1;
			bMax = bMin + binSize - 1;
		}
		fprintf(fh,"%i\t%i\t%lu\n",bMin,bMax,sbpc_hist[b]);
	}
	fclose (fh);
	return 0;
}


unsigned short int histsbpc::getHashValue(const char *currChrom){
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

const char * histsbpc::getKey(unsigned short int chrom){
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


int histsbpc::signalnoise(const char * sbpcFile,const char * outfile,const char* twoBitFile,int binSize){

	FILE * tempTB;
	time_t rtime;
	struct tm *timeinfo;
	int max=500;
	int min=0;
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
	
	cout << "\nGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
	int s = 1;
	sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
	while(s != 0){
		s = system(sysCall);
		if(s != 0)
			cout << "Trying twoBitInfo again, s is" << s << endl;
	}
	remove(tInfo);
	tempTB = fopen(tChrSize,"r");
	char tbChrom[128];
	unsigned long int tbStart;
	while(!feof(tempTB)){
		int ret = fscanf(tempTB,"%s%lu",tbChrom,&tbStart);
		if(ret == 2){
			string sChr(tbChrom);
			unsigned short int chromIntTB = getHashValue(sChr.c_str());
			chr_size[chromIntTB] = tbStart;
		}
	}
	fclose(tempTB);
	remove(tChrSize);
	int dataRange = max - min;
	//int binSize;
	//if((dataRange % numBin) != 0)
	int numBin=1;


	//if(dataRange > numBin)
	//	binSize = (int)dataRange/numBin;
	//else
	//	binSize = 1;
	
	int ltMin = 0;
	int gtMax = 0;
	FILE * fh;
	fh = fopen(outfile,"w");
	if(fh==NULL){
		cout << "Unable to open output file: " << outfile << endl;
		return 1;
	}
	fclose (fh);
	unsigned short int * sbpc_count = NULL;
	string line;string field;
	string chrom;int chromInt;
	int countBases = 0;
	cout << "Getting data for signal to noise analysis  from " << sbpcFile << endl;
	ifstream seqfile(sbpcFile);
	while (getline(seqfile, line)){
		if (line[0] != 't' && line[0] != 'f'){
			countBases++;
			sbpc_count[countBases] = atoi(line.c_str());
		}else if (line[0] == 'f'){
			if(countBases > 0){
				//sort(sbpc_count, sbpc_count+(chr_size[chromInt]+1));
				
				//int nonzeroInd = 0;
				//while(sbpc_count[nonzeroInd] == 0)
				//	nonzeroInd++;
				
				int sInd = 0;
				int bMin = min;
				int bMax = bMin + binSize;
				int window[binSize];
				int binNum=(int) chr_size[chromInt]/binSize;
				unsigned long int * sbpc_hist = NULL;
				unsigned long int * sbpc_hist2 = NULL;
				sbpc_hist = new unsigned long int [binNum];
				for(int b = binNum;b--;)
					sbpc_hist[b] = 0;
				sbpc_hist2 = new unsigned long int [binNum];
				for(int b = binNum;b--;)
					sbpc_hist2[b] = 0;
				//for(int b = 0; b < numBin; b++){
				while(bMax<countBases && bMax <= chr_size[chromInt]){
					for(int i=bMin;i<bMax;i++){
						window[i-bMin]=sbpc_count[i];
					}	
					sort(window, window+binSize);				
					int nonzeroInd = 0;
					while(window[nonzeroInd] == 0) nonzeroInd++;
					
					int nzlength = binSize - nonzeroInd;
					int medInd = (int) nzlength/2;
					sbpc_hist[sInd]=window[nonzeroInd+medInd-1];					
					sbpc_hist2[sInd]=window[binSize-1];
					sInd++;
					bMin = bMax ;
					bMax = bMin + binSize;
					
				}
				
				fh = fopen(outfile,"a");
				if(fh==NULL){
					cout << "Unable to open output file: " << outfile << endl;
					return 1;
				}
				for(int b = 0; b < binNum; b++){
					fprintf(fh,"%lu\t%lu\n",sbpc_hist[b],sbpc_hist2[b]);
				}
				fclose (fh);
				countBases = 0;
				delete [] sbpc_count;
				sbpc_count = NULL;
			}
			
			istringstream iss(line);
			while(iss >> field){
				if (field[0] == 'c'){
					chrom = "";
					for (int it= 6; it < field.length(); it++ )
						chrom += field.at(it);
					chromInt = getHashValue(chrom.c_str());
					sbpc_count = new unsigned short int[chr_size[chromInt]+1];
					for(int c = chr_size[chromInt]; c--;)
						sbpc_count[c] = 0;
				}else if (field[0] == 's' && field[2] == 'a'){
					string startVal;
					for (int it= 6; it < field.length(); it++)
						startVal += field.at(it);
					int startOffset = atoi(startVal.c_str());
					countBases = startOffset - 1;
				}
			}
		}
	}
	seqfile.close();

	if(countBases > 0){
				//sort(sbpc_count, sbpc_count+(chr_size[chromInt]+1));
				
				//int nonzeroInd = 0;
				//while(sbpc_count[nonzeroInd] == 0)
				//	nonzeroInd++;
				
				int sInd = 0;
				int bMin = min;
				int bMax = bMin + binSize;
				int window[binSize];
				int binNum=(int) chr_size[chromInt]/binSize;
				unsigned long int * sbpc_hist = NULL;
				unsigned long int * sbpc_hist2 = NULL;
				sbpc_hist = new unsigned long int [binNum];
				for(int b = binNum;b--;)
					sbpc_hist[b] = 0;
				sbpc_hist2 = new unsigned long int [binNum];
				for(int b = binNum;b--;)
					sbpc_hist2[b] = 0;
				//for(int b = 0; b < numBin; b++){
				while(bMax<countBases && bMax <= chr_size[chromInt]){
					for(int i=bMin;i<bMax;i++){
						window[i-bMin]=sbpc_count[i];
					}	
					sort(window, window+binSize);				
					int nonzeroInd = 0;
					while(window[nonzeroInd] == 0) nonzeroInd++;
					
					int nzlength = binSize - nonzeroInd;
					int medInd = (int) nzlength/2;
					sbpc_hist[sInd]=window[nonzeroInd+medInd-1];					
					sbpc_hist2[sInd]=window[binSize-1];
					sInd++;
					bMin = bMax ;
					bMax = bMin + binSize;
					
				}
				
				fh = fopen(outfile,"a");
				if(fh==NULL){
					cout << "Unable to open output file: " << outfile << endl;
					return 1;
				}
				for(int b = 0; b < binNum; b++){
					fprintf(fh,"%lu\t%lu\n",sbpc_hist[b],sbpc_hist2[b]);
				}
				fclose (fh);
				countBases = 0;
				delete [] sbpc_count;
				sbpc_count = NULL;
			
			}
	
	

	
	
	return 0;
}

