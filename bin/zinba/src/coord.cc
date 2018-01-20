#include "coord.h"

coord::coord(){}

coord::~coord(){}

coord::coord(const coord& c){
	chrom = c.chrom;
	start = c.start;
	end = c.end;
	qFlag = c.qFlag;
	sigVal = c.sigVal;
}

coord::coord(unsigned short int _chrom, unsigned long int _start, unsigned long int _end, unsigned short int _qFlag,double _sigVal){
	chrom = _chrom;
	start = _start;
	end = _end;
	qFlag = _qFlag;
	sigVal = _sigVal;
}

bool coord::operator <(coord other) const{
	if(chrom==other.chrom){
		if( start < other.start ){
			return true;
		}else{
			return false;	
		}
	}else{
		return (chrom<other.chrom);	
	}
}
