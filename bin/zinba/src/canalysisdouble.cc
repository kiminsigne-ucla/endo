#include <iostream>
#include <fstream>
#include <stdio.h>
#include "canalysisdouble.h"
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

using namespace std;

canalysisdouble::canalysisdouble(){
	chromCounter = 0;
}

canalysisdouble::~canalysisdouble(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int canalysisdouble::processCoords(const char* inputFile,const char* outputFile,const char* twoBitFile){

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
	
	unsigned short int chromInt;
	string line;
	string field;

	int printFlag = 0;
	unsigned short int collectData = 0;
	list<coord>::iterator i = coord_slist.begin();
	double * basepair = NULL;
	double * profile = NULL;
	unsigned long int countBases = 250000000;
	unsigned long int startOffset = 0;
	cout << "Getting basecount data from " << inputFile << endl;
	ifstream seqfile(inputFile);
	if(!seqfile.is_open()){
		cout << "Error opening input file" << inputFile << ", exiting" << endl;
		return 1;
	}

	while (getline(seqfile, line)){
		if (line[0] == 'f'){
			if(collectData == 1){
				i = coord_slist.begin();
				while(i!=coord_slist.end()){
					if(i->chrom == chromInt){
						int pIndex = 0;
						profile = new double[(i->end-i->start)+1];
						for(int s = i->start; s <= i->end; s++){
							profile[pIndex] = basepair[s];
							pIndex++;
						}
						if(outputData( outputFile,printFlag,i->chrom,i->start,i->end,i->qFlag,i->sigVal,pIndex,profile) != 0){
							cout << "Error printing output to file, exiting" << endl;
							return 1;
						}
						delete [] profile;
						profile = NULL;
						printFlag = 1;
						coord_slist.erase(i++);
					}else{
						i++;
					}
				}
				delete [] basepair;
				basepair = NULL;

				if(coord_slist.empty()){
					seqfile.close();
					return 0;
				}
			}
			
			istringstream iss(line);
			while(iss >> field){
				if (field[0] == 'c'){
					string chrom;
					for ( int it= 6; it < field.length(); it++ ){
						chrom += field.at(it);
					}
					chromInt = getHashValue(chrom.c_str());
					cout << "Loading data for " << chrom.c_str() << endl;
					basepair = new double[chr_size[chromInt]+1];
					for(int c = chr_size[chromInt];c--;)
						basepair[c] = 0.0;
					countBases = 0;

				}else if (field[0] == 's' && field[2] == 'a'){
					string startVal;
					for ( int it= 6; it < field.length(); it++ ){
						startVal += field.at(it);
					}
					startOffset = atoi(startVal.c_str());
					countBases = startOffset - 1;
				}
			}
		}else if(line[0] != 't'){
				collectData = 1;
				countBases++;
				basepair[countBases] = atof(line.c_str());
				if(countBases > chr_size[chromInt]){
					cout << "Print some error, adding more data than basepairs in chrom\n";
				}
		}
	}

	i = coord_slist.begin();
	while(i!=coord_slist.end()){
		if(i->chrom == chromInt){
			int pIndex = 0;
			profile = new double[(i->end-i->start)+1];
			for(int s = i->start; s <= i->end; s++){
				profile[pIndex] = basepair[s];
				pIndex++;
			}
			if(outputData( outputFile,printFlag,i->chrom,i->start,i->end,i->qFlag,i->sigVal,pIndex,profile) != 0){
				cout << "Error printing output to file, exiting" << endl;
				return 1;
			}
			delete [] profile;
			profile = NULL;
			printFlag = 1;
			coord_slist.erase(i++);
		}else{
			i++;
		}
	}
	seqfile.close();
	delete [] basepair;
	basepair = NULL;
	return 0;
}

int canalysisdouble::outputData(const char * outputFile,int pFlag,unsigned short int pChrom,unsigned long int pStart,unsigned long int pStop,unsigned short int sFlag,double pSigVal,int printStop,double pProfile[]){
	FILE * fh;
	if(pFlag == 0){
		fh = fopen(outputFile,"w");
		if(fh==NULL){
			cout << "Unable to open output file: " << outputFile << endl;
			return 1;
		}
		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND\tSIGVAL");
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
	string strand;
	if(sFlag == 1)
		strand = "+";
	else
		strand = "-";
	char winID[255];
	char start[10];
	char stop[10];
	sprintf( start,"%lu", pStart);
	sprintf( stop,"%lu", pStop);
	strcpy(winID,chromName);strcat(winID,":");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
	fprintf(fh,"%s\t%s\t%lu\t%lu\t%s\t%.5f",winID,chromName,pStart,pStop,strand.c_str(),pSigVal);
	if(sFlag == 1){
		for(int posP = 0; posP < printStop;posP++){
			fprintf(fh,"\t%f",pProfile[posP]);
		}
	}else{
		for(int posP = printStop-1; posP >= 0;posP--){
			fprintf(fh,"\t%f",pProfile[posP]);
		}
	}
	fprintf(fh,"\n");
	fclose (fh);
	return 0;
}

unsigned short int canalysisdouble::getHashValue(const char *currChrom){
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

const char * canalysisdouble::getKey(unsigned short int chrom){
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

int canalysisdouble::importCoords(const char *coordfile){
	char cChrom[128];
	char cid[128];
	char strand[1];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	double score;
	char firstline [256];
	int rwline;
	char pstrand[] = "+";

	FILE * cfile;
	cfile = fopen(coordfile,"r");
	if(cfile == NULL){
		cout << "ERROR: opening window list file" << endl;
		return 1;
	}
	
//	fgets(firstline,256,cfile);
	while(!feof(cfile)){
		rwline = fscanf(cfile,"%s%s%lu%lu%s%lf",cid,cChrom,&iStart,&iEnd,strand,&score);
		string chromIn(cChrom);
		unsigned short int chromInt = getHashValue(chromIn.c_str());
		if(strcmp(strand,pstrand) == 0)
			qFlag = 1;
		else
			qFlag = 0;
		coord c(chromInt,iStart,iEnd,qFlag,score);
		coord_slist.push_back(c);
	}
	fclose(cfile);
	
	coord_slist.sort();
	cout << "\nImported " << coord_slist.size() << " coordinates....." << endl;
	return 0;
}
