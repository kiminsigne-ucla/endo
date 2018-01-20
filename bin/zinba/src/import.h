#ifndef IMPORT_H_
#define IMPORT_H_

#include <string>
#include <vector>
#include "bwRead.h"
#include "bwRead2.h"
#include "dataWins.h"
#include "cnvWins.h"
#include <cstring>
//#include "base.h"
#include <map>
#include <list>
#include <ext/slist>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;
class import{
	
	public:

		import(); //Implemented
		~import(); //Implemented
		
	
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

		int importBowtie(const char *,int,int);
		int importTagAlign(const char *,int,int);
		int importBed(const char *,int,int);
		int importCustomBed(const char * signalFile,int extension);
		vector<bwRead> signal_slist;//Implemented
		vector<bwRead> input_slist;
		vector<bwRead2> custom_slist;
	
	/*	list<cnvWins> cnv_wins;
		slist<dataWins> peak_wins;
		slist<dataWinsCount> peak_wins2;
		slist<dataWinsCustom> peak_wins3;*/
		unsigned short int tbSizeFlag;
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;
};

#endif /*IMPORT_H_*/
