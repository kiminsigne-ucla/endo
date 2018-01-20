#ifndef COORD_H_
#define COORD_H_

#include <iostream>

class coord {
	
	public:
			coord();
			~coord();
			coord(unsigned short int, unsigned long int, unsigned long int,unsigned short int,double);
			coord(const coord&);
			bool operator <(coord) const;

			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int end;//8
			unsigned short int qFlag;
			double sigVal;
	
	private:
			
};

#endif /*COORD_H_*/
