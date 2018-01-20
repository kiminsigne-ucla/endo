#ifndef FRPANALYSIS_H_
#define FRPANALYSIS_H_

#include <string>
#include <cstring>
#include <vector>
#include "bwRead.h"
#include "coord2.h"
//#include "base.h"
#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class frpanalysis{
	
	public:
	
		frpanalysis(); //Implemented
		~frpanalysis(); //Implemented
		
		int importCoords(const char *,double,const char *,int, int, int,const char *);//Implemented
		int importPscl(const char *,double,int);
		int importMixture(const char *,double,int,int);

		int importRawSignal(const char *,int,const char *);//Implemented
		int importTagAlign(const char *,int);
		int importBowtie(const char *,int);
		int importBed(const char *,int);

		int processSigCoords(const char *,int,const char *,const char *,const char *);//Implemented


		int outputData(const char *, unsigned short int,unsigned short int,unsigned int[]);//Implemented

	
	
		struct ltstr{
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
		};

	private:

		list<coord2> coordIN_slist;//Implemented
		list<coord2> coordOUT_slist;//Implemented
	
//		slist<bwRead> signal_slist;//Implemented
		vector<bwRead> signal_slist;//Implemented
		unsigned short int tbSizeFlag;
		
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(const char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
};

#endif /*FRPANALYSIS_H_*/
