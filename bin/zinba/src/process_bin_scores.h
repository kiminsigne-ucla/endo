#ifndef PROCESS_BIN_SCORES_H_
#define PROCESS_BIN_SCORES_H_

#include <string>
#include <map>
using namespace std;

class process_bin_scores{
	
	public:
	
		process_bin_scores();
		~process_bin_scores();
		int adjustCoords(const char *,string,const char *,int,int);

	private:
	
		map<string,unsigned long int> chr_size;
		
};

#endif
