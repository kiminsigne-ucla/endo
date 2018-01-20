#include <iostream>
#include <fstream>
#include <stdio.h>
#include "coord2.h"
#include "pc.h"
#include <sstream>
#include <string>
#include <string.h>
#include <cstring>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <list>
#include <ext/hash_map>

extern "C" {
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 

extern "C" {
void collapse_windows(const char ** Rwinlist, const char ** Rmethod, int * Rwformat, int * 
RlengthThresholds,double thresholds[], int *RwinGap, int *RFDR, const char ** Routput){
	hash_map<string, int, hash<string>,equal_to<string> > cind_map;
	vector<string> cnames;

	const char * winlist = Rwinlist[0];
	const char * method = Rmethod[0];
	const char * output = Routput[0];
	int wformat = Rwformat[0];
	int winGap = RwinGap[0];
	int FDR=RFDR[0];
	int lengthThresholds = RlengthThresholds[0];
	
	string winstring = string(winlist);
	//size_t found = winstring.find_last_of(".");
	//string outfile = winstring.substr(0,found);
	string output2= string(output);

	FILE * wlist;
	wlist = fopen(winlist,"r");
	if(wlist == NULL){error("Unable to open list file: %s\n", winlist);}
	char cChrom[128];
	unsigned long int iStart;
	unsigned long int iEnd;
	unsigned short int qFlag;
	double sigVal,qVal;
	const char *pscl = "pscl";
	const char *mixture = "mixture";
	char sigFile [256];
	int readResult = 0;
	char firstline [256];
	int rwline;
	double sres;int ec;double ic;double gc;double ap;double cnv;
	string chr;
	///////////////////////////
	unsigned long int winSizeThresh;
	int winSize;
	//////////////////////////
	list<coord2> coordIN_slist;
	list<coord2> coordOUT_slist;
	
	double highThresh;
	if(strcmp(method,pscl) == 0){
		highThresh = 0;
		for(int t = 0; t < lengthThresholds; t++){
			if(thresholds[t] > highThresh)
				highThresh = thresholds[t];
		}
	}else if(strcmp(method,mixture) == 0 & FDR==0){
		highThresh = 1;
		for(int t = 0; t < lengthThresholds; t++){
			if(thresholds[t] < highThresh)
				highThresh = thresholds[t];
		}
	}else if(strcmp(method,mixture) == 0 & FDR==1){
                highThresh = 0;
                for(int t = 0; t < lengthThresholds; t++){
                        if(thresholds[t] > highThresh)
                                highThresh = thresholds[t];
                }
        }

	cout << "There are " << lengthThresholds << " thresholds: " << endl;
	for(int t = 0; t < lengthThresholds; t++){
		cout << thresholds[t] << " ";
	}
	cout << "\nHigh threshold= " << highThresh << endl;
	
	while(!feof(wlist)){
		int rline = fscanf(wlist,"%s",sigFile);
		if(rline == 1){
			string signalFile(sigFile);
			FILE * fh;
			fh = fopen(signalFile.c_str(),"r");
			if(fh == NULL){error("Unable to open input file: %s\n", signalFile.c_str());}
			cout << "\tImporting windows from " << signalFile.c_str() << "..." << endl;
			fgets(firstline,256,fh);
			while(!feof(fh)){
				readResult = 0;
				if(strcmp(method,pscl) == 0){
					if(wformat == 0){
						rwline =fscanf(fh,"%s%lu%lu%hu%lf%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal,&sres);
						if(sigVal <= highThresh && rwline == 6)
							readResult = 1;
					}else if(wformat == 1){
						rwline = fscanf(fh,"%s%lu%lu%d%lf%lf%lf%lf%hu%lf%lf",cChrom,&iStart,&iEnd,&ec,&ic,&gc,&ap,&cnv,&qFlag,&sigVal,&sres);
						if(sigVal <= highThresh && rwline == 11)
							readResult = 1;
					}
					qVal=sigVal;
				}else if(strcmp(method,mixture) == 0 & FDR==0){
					if(wformat == 0){
						rwline = fscanf(fh,"%s%lu%lu%hu%lf%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal,&qVal);
						if(sigVal >= highThresh && rwline == 6)
							readResult = 1;
					}else if(wformat == 1){
						int exp;double inp,gc,ap,ecl;
						rwline = fscanf(fh,"%s%lu%lu%d%lf%lf%lf%lf%hu%lf%lf",cChrom,&iStart,&iEnd,&ec,&ic,&gc,&ap,&cnv,&qFlag,&sigVal,&qVal);
						if(sigVal >= highThresh && rwline== 11)
							readResult = 1;
					}
				}else if(strcmp(method,mixture) == 0 & FDR==1){
                                        if(wformat == 0){
                                                rwline = fscanf(fh,"%s%lu%lu%hu%lf%lf",cChrom,&iStart,&iEnd,&qFlag,&sigVal,&qVal);
                                                if(qVal <= highThresh && rwline == 6)
                                                        readResult = 1;
                                        }else if(wformat == 1){
                                                int exp;double inp,gc,ap,ecl;
                                                rwline = fscanf(fh,"%s%lu%lu%d%lf%lf%lf%lf%hu%lf%lf",cChrom,&iStart,&iEnd,&ec,&ic,&gc,&ap,&cnv,&qFlag,&sigVal,&qVal);
                                                if(qVal <= highThresh && rwline== 11)
                                                        readResult = 1;
                                        }
                                }

				if(readResult == 1){
					if(qFlag == 1){
						// determine the chromosome index
						chr = string(cChrom);
						hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
						int chromInt=-1;
						if(li==cind_map.end()) {
							// register new chromosome
							chromInt=cnames.size();
							cnames.push_back(chr);
							cind_map[chr]=chromInt;
						}else{
							chromInt=li->second;
						}
						coord2 c(chromInt,iStart,iEnd,qFlag,sigVal,qVal);
						coordIN_slist.push_back(c);
					}
				}
			}
			fclose(fh);
		}
	}
	fclose(wlist);
	winSize = iEnd - iStart;
	winSizeThresh = winSize * 10;
	cout << "\nImported " << coordIN_slist.size() << " coordinates" << endl;
	coordIN_slist.sort();
	list<coord2> tempCoords;	
	list<coord2>::iterator back;
	
	for(int t = 0; t < lengthThresholds; t++){
		cout << "Processing threshold " << thresholds[t] << endl;
		back =  coordIN_slist.begin();
		while(back != coordIN_slist.end()){
			if(strcmp(method,pscl) == 0){
				if(back->sigVal <= thresholds[t])
					tempCoords.push_back(*back);
			}else if(strcmp(method,mixture) == 0){
				if(back->sigVal >= thresholds[t])
					tempCoords.push_back(*back);
			}
			back++;
		}
		tempCoords.sort();
		cout << "\t" << tempCoords.size() << " coords at this threshold" << endl;
		back = tempCoords.begin();
		coord2 tempCoord = *back;
		int flagEnd = 0;
		back++;
		while(flagEnd == 0){
			if(back == tempCoords.end())
				flagEnd = 1;
			if(back->chrom == tempCoord.chrom && back->start <= tempCoord.end+1+winGap && flagEnd == 0){
				tempCoord.end = back->end;
				if(strcmp(method,pscl) == 0 && tempCoord.sigVal > back->sigVal){
					tempCoord.sigVal = back->sigVal;
				}else if(strcmp(method,mixture) == 0 && tempCoord.sigVal < back->sigVal && FDR==0){
					tempCoord.sigVal = back->sigVal;
					tempCoord.qVal = back->qVal;
				}else if(strcmp(method,mixture) == 0 && tempCoord.qVal < back->sigVal & FDR==1){
                                        tempCoord.sigVal = back->sigVal;
                                        tempCoord.qVal = back->qVal; 
                                }
			}else{
				//if((tempCoord.end-tempCoord.start) <= winSizeThresh){
					coordOUT_slist.push_back(tempCoord);
				//}else{
					//REMOVE ONCE C PEAKBOUND IS RUNNING
				//	string cName = cnames[tempCoord.chrom];
				//	cout << "Excluding " << cName.c_str() << ":" << tempCoord.start << "-" << tempCoord.end << " SIZE=" << (tempCoord.end-tempCoord.start) << endl;
				//}
				tempCoord = *back;
			}
			if(flagEnd == 0)
				back++;
		}
		coordOUT_slist.sort();
		cout << "\tCollapsed to " << coordOUT_slist.size() << " regions" << endl;
	
		///output data
		FILE * fh;
		char tval[128];
		sprintf(tval,"%.14lf",thresholds[t]);
		string tv = string(tval);
//		string outputFile = outfile + "_" + tv + ".coords";
//  		string outputFile = outfile + "_" + tv + ".coords.bed";
		fh = fopen(output2.c_str(),"w");
		if(fh==NULL){error("Unable to open output file: %s\n", output2.c_str());}
//		fprintf(fh,"COORDID\tCHROM\tSTART\tSTOP\tSTRAND\tSIGVAL\n");
//		fprintf(fh,"CHROM\tSTART\tSTOP\tSTRAND\tSIGVAL\n");
		string chromName;
		char strand[] = "+";
//		char winID[255];
		char start[10];
		char stop[10];
		back = coordOUT_slist.begin();
		while(back != coordOUT_slist.end()){
			chromName = cnames[back->chrom];
			sprintf( start,"%lu", back->start);
			sprintf( stop,"%lu", back->end);
//			strcpy(winID,chromName.c_str());strcat(winID,":");strcat(winID,start);strcat(winID,"-");strcat(winID,stop);
//			fprintf(fh,"%s\t%s\t%lu\t%lu\t%s\t%.14f\n",winID,chromName.c_str(),back->start,back->end,strand,back->sigVal);
			fprintf(fh,"%s\t%lu\t%lu\t%.14f\t%.14f\%s\n",chromName.c_str(),back->start,back->end,back->sigVal,back->qVal,strand);
			back++;
		}
		fclose (fh);
		tempCoords.clear();
		coordOUT_slist.clear();
	}
}
}
