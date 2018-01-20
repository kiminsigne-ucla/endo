#ifndef WINCOUNTS_H_
#define WINCOUNTS_H_

#include <string>
#include <vector>
#include "bwRead.h"
#include "coord.h"
#include <cstring>
//#include "base.h"
#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class winCounts{
	
	public:
	
		winCounts(); //Implemented
		~winCounts(); //Implemented
		
		int importRawSignalWC(const char *,int,const char *);//Implemented
		int importCoordsWC(const char *);//Implemented
		int processSignalsWC(const char *,const char *);//Implemented
		int outputDataWC(const char *, unsigned short int,int);//Implemented
		int importBowtieWC(const char *,int);

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

		vector<bwRead> signal_slist;//Implemented
		list<coord> coord_slist;//Implemented
		unsigned short int tbSizeFlag;
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
};

#endif /*WINCOUNTS_H_*/
