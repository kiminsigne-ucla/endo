#ifndef CALCCOVS_H_
#define CALCCOVS_H_

#include <string>
//#include <vector>
//#include "bwRead.h"
//#include "bwRead2.h"
#include "dataWins.h"
#include "import.h"
#include "cnvWins.h"
#include <cstring>
//#include "base.h"
//#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class calcCovs{
	
	public:
	
		//calcCovs(); //Implemented
		//~calcCovs(); //Implemented
		import b;
		int importRawSignal(const char *,int,const char *,int,const char *);//Implemented
		int processSignals(int zWinSize, int zOffsetSize, int cWinSize, int cOffsetSize, string alignDir,const char * twoBitFile,const char * inputFile,string outfile,const char * flist,int extension,const char * filetype, int binary);//Implemented
    int processSignalsfast(int zWinSize, int zOffsetSize,int cWinSize, int cOffsetSize, const char * twoBitFile,const char * inputFile,string outfile,const char * flist,int extension,const char * filetype);
		int outputData(const char *, unsigned short int);//Implemented
		int outputDataWinCount(const char *, unsigned short int);
/*		int importBowtie(const char *,int,int);
		int importTagAlign(const char *,int,int);
		int importBed(const char *,int,int);*/
		int processWinSignal(int , int ,const char * ,string outfile,int ,const char *, double);
		int	processCustomSignal(int zWinSize, int zOffsetSize,const char * twoBitFile,const char * inputFile,string outfile, const char* flist, int extension);
		//int importCustomBed(const char * signalFile,int extension);
		int outputCustomData(const char * outputFile, unsigned short int currChr);
		/*struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};
	
		struct ltint{
			bool operator()(const int i1, const int i2) const
			{
				return i1<i2;
		 	}
		};*/
	private:
		//commented out items have moved to import.h
/*	vector<bwRead> signal_slist;//Implemented
		vector<bwRead> input_slist;
		vector<bwRead2> custom_slist;
//	list<bwRead> signal_slist;//Implemented
//	slist<bwRead> input_slist;*/
		list<cnvWins> cnv_wins;
		slist<dataWins> peak_wins;
		slist<dataWinsCount> peak_wins2;
		slist<dataWinsCustom> peak_wins3;
		//unsigned short int tbSizeFlag;
	
		//unsigned short int chromCounter;//Implemented
		//unsigned short int getHashValue(char *);//Implemented
		//const char * getKey(unsigned short int);//Implemented
		//map<const char*, int, ltstr> chroms;//Implemented
		//map<int, const char*> intsToChrom;//Implemented
	
//		map<unsigned short int,unsigned long int> chr_size;
		 //map<unsigned short int,unsigned long int> chr_size;
};

#endif /*CALCCOVS_H_*/
