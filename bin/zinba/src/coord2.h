#ifndef COORD_H_
#define COORD_H_

#include <iostream>

class coord2 {
	
	public:
			coord2();
			~coord2();
			coord2(unsigned short int, unsigned long int, unsigned long int,unsigned short int,double,double);
			coord2(const coord2&);
			bool operator <(coord2) const;

			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int end;//8
			unsigned short int qFlag;
			double sigVal;
			double qVal;
	private:
			
};

#endif /*COORD_H_*/
