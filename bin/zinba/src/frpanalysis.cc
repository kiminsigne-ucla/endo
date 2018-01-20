#include <iostream>
#include <fstream>
#include <stdio.h>
#include "frpanalysis.h"
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
using namespace std;

frpanalysis::frpanalysis(){
	chromCounter = 0;
}

frpanalysis::~frpanalysis(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int frpanalysis::processSigCoords(const char* bcountFile,int extend,const char* twoBitFile,const char* chrmSearch,const char* outputFile){

	string nflag = "none";
	int profile_extend = 500;
	const char * chromReport;
	unsigned int * basepair = NULL;
	unsigned int * profile = NULL;
	unsigned short int printFLAG = 0;
	unsigned short int currchr;
	unsigned long int countBases = 250000000;
	unsigned long int startOffset = 0;
	int getChrmData = 0;
	int collectData = 0;
	string allChrm = "all";
	string line;
	string field;
	list<coord2>::iterator i = coordOUT_slist.begin();
	int rnum;
	ifstream seqfile;
	
	if(strcmp(bcountFile,nflag.c_str()) == 0){
		cout << "\nGetting basecount data from reads ....." << endl;
		cout << "\tProcessing " << signal_slist.size() << " reads" << endl;
	}else{
		cout << "\nGetting basecount data from " << bcountFile << endl;
		ifstream seqfile(bcountFile);
		if(!seqfile.is_open()){
			cout << "ERROR opening input file" << bcountFile << ", exiting" << endl;
			return 1;
		}
	}
	
	while(!coordOUT_slist.empty()){
	
		if(strcmp(bcountFile,nflag.c_str()) == 0){
			unsigned long int aStart,aStop;
			rnum = 0;
			
			while(rnum < (int) signal_slist.size()){
				chromReport = getKey(signal_slist[rnum].chrom);
				if(strcmp(chrmSearch,allChrm.c_str()) == 0){
					getChrmData = 1;
					currchr = signal_slist[rnum].chrom;
				}else if (strcmp(chrmSearch,chromReport) == 0){
					getChrmData = 1;
					currchr = signal_slist[rnum].chrom;
				}else{
					getChrmData = 0;
				}

				if(getChrmData == 1){
					cout << "\nProcessing " << chromReport << ".........." << endl;
					cout << "\tInitializing length to " << chr_size[currchr] << endl;
					basepair = new unsigned int[chr_size[currchr]+1];
					for(int in = chr_size[currchr]; in--;)
						basepair[in] = 0;

					while(signal_slist[rnum].chrom==currchr){
						aStart = signal_slist[rnum].pos;
						aStop = chr_size[currchr];
						if((aStart+extend-1) <= chr_size[currchr])
							aStop = aStart + extend - 1;
						for(unsigned int rpos = aStart; rpos <= aStop;rpos++)
							basepair[rpos]++;
						rnum++;
					}
					cout << "\t" << rnum << " reads mapped to " << chromReport << endl;
					signal_slist.erase(signal_slist.begin(),signal_slist.begin()+rnum);
					rnum = (int) signal_slist.size();
				}else{
					rnum++;
				}
			}
			
		}else{

			getChrmData = 0;
			if(collectData == 1){
				istringstream iss(line);
				while(iss >> field){
					if (field[0] == 'c'){
						string chrom;
						for ( int it= 6; it < field.length(); it++ ){
							chrom += field.at(it);
						}
						if(strcmp(chrmSearch,allChrm.c_str()) == 0){
							getChrmData = 1;
							currchr = getHashValue(chrom.c_str());
						}else if (strcmp(chrmSearch,chrom.c_str()) == 0){
							getChrmData = 1;
							currchr = getHashValue(chrom.c_str());
						}else{
							getChrmData = 0;
						}
						
						if(getChrmData == 1){
							cout << "\nProcessing for " << chrom.c_str() << endl;
							cout << "\tInitializing length of " << chr_size[currchr] << endl;
							basepair = new unsigned int[chr_size[currchr]+1];
							for(int c = chr_size[currchr];c--;)
								basepair[c] = 0;
							countBases = 0;
						}
					}else if (field[0] == 's' && field[2] == 'a' && getChrmData == 1){
						string startVal;
						for ( int it= 6; it < field.length(); it++ ){
							startVal += field.at(it);
						}
						startOffset = atoi(startVal.c_str());
						cout << "\tStarting position is " << startOffset << endl;
						countBases = startOffset - 1;
					}
				}				
			}
			collectData = 0;
			
			while (getline(seqfile,line) && collectData == 0){
				if (line[0] == 'f' && getChrmData == 0){
					istringstream iss(line);
					while(iss >> field){
						if (field[0] == 'c'){
							string chrom;
							for ( int it= 6; it < field.length(); it++ ){
								chrom += field.at(it);
							}
							if(strcmp(chrmSearch,allChrm.c_str()) == 0){
								getChrmData = 1;
								currchr = getHashValue(chrom.c_str());
							}else if (strcmp(chrmSearch,chrom.c_str()) == 0){
								getChrmData = 1;
								currchr = getHashValue(chrom.c_str());
							}else{
								getChrmData = 0;
							}
							
							if(getChrmData == 1){
								cout << "\nProcessing for " << chrom.c_str() << endl;
								cout << "\tInitializing length of " << chr_size[currchr] << endl;
								basepair = new unsigned int[chr_size[currchr]+1];
								for(int c = chr_size[currchr];c--;)
									basepair[c] = 0;
								countBases = 0;
							}
						}else if (field[0] == 's' && field[2] == 'a' && getChrmData == 1){
							string startVal;
							for ( int it= 6; it < field.length(); it++ ){
								startVal += field.at(it);
							}
							startOffset = atoi(startVal.c_str());
							cout << "\tStarting position is " << startOffset << endl;
							countBases = startOffset - 1;
						}
					}
				}else if(line[0] == 'f' && getChrmData == 1){
					collectData = 1;
				}else if(getChrmData == 1 && collectData == 0){
					countBases++;
					if(countBases > chr_size[currchr]){
						cout << "\nERROR: Current position is " << countBases << " and chromosome length is " << chr_size[currchr] << endl;
						return 1;
					}
					basepair[countBases] = atoi(line.c_str());
				}
			}
		}

		i = coordOUT_slist.begin();
		while(i!=coordOUT_slist.end()){
			if(i->chrom == currchr){
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// THE PROFILE ARRAY CONTAINS THE BASECOUNT DATA FOR SIGREGION //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
//				if(outputData( outputFile,printFlag,i->chrom,startPos,stopPos,i->sigVal,i->qVal,pIndex,profile) != 0){
//					cout << "Error printing output to file, exiting" << endl;
//					return 1;
//				}
				delete [] profile;
				profile = NULL;
//				printFlag = 1;
				coordOUT_slist.erase(i++);
			}
		}

		delete [] basepair;
		basepair = NULL;
	}

	if(strcmp(bcountFile,nflag.c_str()) != 0){
		seqfile.close();
	}
	return 0;
}

int frpanalysis::outputData(const char * outputFile, unsigned short int currChr,unsigned short int pFLAG,unsigned int basepair[]){
	FILE * fh;
	const char * chrom = getKey(currChr);
	if(pFLAG == 0){
		fh = fopen(outputFile,"w");
		fprintf(fh,"track type=wiggle_0 name=\"%s\" desc=\"%s\" visibility=full\n",outputFile,outputFile);
	}else{
		fh = fopen(outputFile,"a");
	}
	if(fh == NULL){
		cout << "\nERROR: Can't open output file: " << outputFile << endl;
		return 1;
	}
	fprintf(fh,"fixedStep chrom=%s start=1 step=1\n",chrom);

	for(int posP = 1; posP <= chr_size[currChr];posP++)
		fprintf(fh,"%u\n",basepair[posP]);
	fclose (fh);
	return 0;
}

unsigned short int frpanalysis::getHashValue(const char *currChrom){
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

const char * frpanalysis::getKey(unsigned short int chrom){
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

int frpanalysis::importRawSignal(const char * readFile,int extension,const char * filetype){
	
	const char * bowtie = "bowtie";
	const char * bed = "bed";
	const char * tagAlign = "tagAlign";
	int rval = 0;
	
	if(strcmp(filetype,bed) == 0){
		rval = importBed(readFile,extension);
	}else if (strcmp(filetype,bowtie) == 0){
		rval = importBowtie(readFile,extension);
	}else if(strcmp(filetype,tagAlign) == 0){
		rval = importTagAlign(readFile,extension);
	}else{
		cout << "Unrecognized type of file " << filetype << ", must be either bowtie, bed, or tagAlign" << endl;
		return 1;
	}

	if(rval == 0){
		cout << "\tLoaded " << signal_slist.size() << " reads" << endl;
		cout << "\tSorting reads ...";
		sort (signal_slist.begin(), signal_slist.end());
		cout << "COMPLETE" << endl;
	}else{
		return 1;
	}

	return 0;
}

int frpanalysis::importBowtie(const char * signalFile,int extension){
	
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

int frpanalysis::importTagAlign(const char * signalFile,int extension){
	
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

int frpanalysis::importBed(const char * signalFile,int extension){
	
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

int frpanalysis::importCoords(const char *winlist,double threshold,const char *method,int wformat, int winGap, int FDR, const char * twoBitFile){
	
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
		cout << "\t\tFor " << cChrom << " length is " << cStart << endl;
		unsigned short int chromInt = getHashValue(cChrom);
		chr_size[chromInt] = cStart;
	}
	fclose(tempTB);
	remove(tChrSize);
	
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
				coordOUT_slist.push_back(tempCoord);
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

int frpanalysis::importPscl(const char *winlist,double threshold,int wformat){
	
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

int frpanalysis::importMixture(const char *winlist,double threshold,int wformat, int FDR){
	
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
