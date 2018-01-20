#ifndef BWREAD_H_
#define BWREAD_H_

#include <iostream>

class bwRead {
	
	public:
			bwRead();
			~bwRead();
//			bwRead(unsigned short int, unsigned long int,unsigned short int);
			bwRead(unsigned short int, unsigned long int);
			bwRead(const bwRead&);
			bool operator <(bwRead) const;
			
			unsigned short int chrom;//2
			unsigned long int pos;//8
//			unsigned short int strand;//Plus = 1;Minus = 0
	
	private:
			
};

#endif /*BWREAD_H_*/
