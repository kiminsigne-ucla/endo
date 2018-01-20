#ifndef PROCESS_SCORES_H_
#define PROCESS_SCORES_H_

#include <string>
#include <map>
using namespace std;

class process_scores{
	
	public:
	
		process_scores();
		~process_scores();
		int adjustCoords(string,string,const char *,int,int);

	private:
	
		map<string,unsigned long int> chr_size;
		
};

#endif
