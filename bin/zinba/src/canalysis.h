#ifndef CANALYSIS_H_
#define CANALYSIS_H_

#include <string>
#include <cstring>
#include <vector>
#include "coord.h"
#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class canalysis{
	
	public:
	
		canalysis(); //Implemented
		~canalysis(); //Implemented
		int importCoords(const char *);//Implemented
		int processCoords(const char *,const char *,const char *);//Implemented
		int outputData(const char *,int,unsigned short int,unsigned long int,unsigned long int,unsigned short int,double,int,unsigned int[]);//Implemented
	
		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};

	private:
		
		list<coord> coord_slist;//Implemented

		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(const char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
		map<unsigned short int,unsigned long int> chr_med;
};

#endif /*CANALYSIS_H_*/
