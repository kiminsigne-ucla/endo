#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <string>
#include <cstring>
#include <vector>
#include "coord2.h"
#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class analysis{
	
	public:
	
		analysis(); //Implemented
		~analysis(); //Implemented
		int importCoords(const char *,double,const char *,int, int, int);//Implemented
		int processCoords(const char *,const char *,const char *,const char *, int);//Implemented
		int outputData(const char *,int,unsigned short int,unsigned long int,unsigned long int,double,double, int,unsigned int[]);//Implemented
		int importPscl(const char *,double,int);
		int importMixture(const char *,double,int,int);
	
		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};

	private:
		
		list<coord2> coordIN_slist;//Implemented
		list<coord2> coordOUT_slist;//Implemented
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(const char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
		map<unsigned short int,unsigned long int> chr_med;
};

#endif /*ANALYSIS_H_*/
