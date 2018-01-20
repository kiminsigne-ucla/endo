#ifndef HISTSBPC_H_
#define HISTSBPC_H_

#include <string>
#include <cstring>
#include <vector>
#include <map>

namespace sgi = ::__gnu_cxx;
using namespace sgi;
using namespace std;

class histsbpc{
	
	public:
	
		histsbpc(); //Implemented
		~histsbpc(); //Implemented
		int hist_data(const char *,const char *,const char *,int,int,int);//Implemented
		int signalnoise(const char *,const char *,const char *,int);//Implemented	
		struct ltstr{
			bool operator()(const char* s1, const char* s2) const
			{
				return strcmp(s1, s2) < 0;
		 	}
		};

	private:

//		vector<unsigned long int> sbpc_hist;
	
		unsigned short int chromCounter;//Implemented
		unsigned short int getHashValue(const char *);//Implemented
		const char * getKey(unsigned short int);//Implemented
		map<const char*, int, ltstr> chroms;//Implemented
		map<int, const char*> intsToChrom;//Implemented
	
		map<unsigned short int,unsigned long int> chr_size;

};

#endif /*HISTSBPC_H_*/
