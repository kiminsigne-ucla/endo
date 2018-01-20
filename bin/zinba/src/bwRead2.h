#ifndef bwRead2_H_
#define bwRead2_H_

#include <iostream>

class bwRead2 {
	
	public:
			bwRead2();
			~bwRead2();
//			bwRead2(unsigned short int, unsigned long int,unsigned short int);
			bwRead2(unsigned short int, unsigned long, int);
			bwRead2(const bwRead2&);
			bool operator <(bwRead2) const;
			
			unsigned short int chrom;//2
			unsigned long int pos;//8
			int score;
//			unsigned short int strand;//Plus = 1;Minus = 0
	
	private:
			
};

#endif /*bwRead2_H_*/
