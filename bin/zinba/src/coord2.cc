#include "coord2.h"

coord2::coord2(){}

coord2::~coord2(){}

coord2::coord2(const coord2& c){
	chrom = c.chrom;
	start = c.start;
	end = c.end;
	qFlag = c.qFlag;
	sigVal = c.sigVal;
	qVal=c.qVal;
}

coord2::coord2(unsigned short int _chrom, unsigned long int _start, unsigned long int _end, unsigned short int _qFlag,double _sigVal, double _qVal){
	chrom = _chrom;
	start = _start;
	end = _end;
	qFlag = _qFlag;
	sigVal = _sigVal;
	qVal=_qVal;
}

bool coord2::operator <(coord2 other) const{
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
