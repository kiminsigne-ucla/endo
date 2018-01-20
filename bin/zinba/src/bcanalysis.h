#ifndef BCANALYSIS_H_
#define BCANALYSIS_H_

#include <string>
#include <cstring>
//#include <vector>
//#include "bwRead.h"
//#include "base.h"
//#include <map>
#include <ext/slist>
#include "import.h"
namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class bcanalysis{
	
	public:
	
		//bcanalysis(); //Implemented
		//~bcanalysis(); //Implemented
		import b;
		int importRawSignal(const char *,int,const char *,const char *);//Implemented
		int processSignals(const char *,int, int);//Implemented
//		int outputData(const char *, unsigned short int,unsigned short int,unsigned short int[]);//Implemented
		int outputData(const char *, unsigned short int,unsigned short int,unsigned int[], int binary);//
		/*int importTagAlign(const char *,int);
		//int importBowtie(const char *,int);
		//int importBed(const char *,int);
	
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
		};*/

	private:
		
//		slist<bwRead> signal_slist;//Implemented
		 //vector<bwRead> signal_slist;//Implemented
		//unsigned short int tbSizeFlag;
		
		//unsigned short int chromCounter;//Implemented
		//unsigned short int getHashValue(char *);//Implemented
		//const char * getKey(unsigned short int);//Implemented
		//map<const char*, int, ltstr> chroms;//Implemented
		//map<int, const char*> intsToChrom;//Implemented
	
//		map<unsigned short int,unsigned long int> chr_size;
		//map<unsigned short int,unsigned int> chr_size;
};

#endif /*BCANALYSIS_H_*/
