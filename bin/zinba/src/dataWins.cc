#include "dataWins.h"

dataWins::~dataWins(){}

dataWins::dataWins(){}

dataWins::dataWins(const dataWins& s){
	chrom = s.chrom;
	start = s.start;
	stop = s.stop;
	eCount = s.eCount;
	iCount = s.iCount;
	gcPerc = s.gcPerc;
	alignPerc = s.alignPerc;
	cnvScore = s.cnvScore;
}

dataWins::dataWins(unsigned short int _chrom, unsigned long int _start, unsigned long int _stop,int _eCount, double _iCount, double _gcPerc, double _alignPerc, double _cnvScore){
	chrom = _chrom;
	start = _start;
	stop = _stop;
	eCount = _eCount;
	iCount = _iCount;
	gcPerc = _gcPerc;
	alignPerc = _alignPerc;
	cnvScore = _cnvScore;
}
bool dataWins::operator <(dataWins other) const{
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

dataWinsCount::~dataWinsCount(){}

dataWinsCount::dataWinsCount(){}

dataWinsCount::dataWinsCount(const dataWinsCount& s){
	chrom = s.chrom;
	start = s.start;
	stop = s.stop;
	eCount = s.eCount;

}

dataWinsCount::dataWinsCount(unsigned short int _chrom, unsigned long int _start, unsigned long int _stop,int _eCount){
	chrom = _chrom;
	start = _start;
	stop = _stop;
	eCount = _eCount;

}
bool dataWinsCount::operator <(dataWinsCount other) const{
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

dataWinsCustom::~dataWinsCustom(){}

dataWinsCustom::dataWinsCustom(){}

dataWinsCustom::dataWinsCustom(const dataWinsCustom& s){
	chrom = s.chrom;
	start = s.start;
	stop = s.stop;
	eCount = s.eCount;
	iCount = s.iCount;
}

dataWinsCustom::dataWinsCustom(unsigned short int _chrom, unsigned long int _start, unsigned long int _stop,int _eCount, double _iCount){
	chrom = _chrom;
	start = _start;
	stop = _stop;
	eCount = _eCount;
	iCount = _iCount;
}
bool dataWinsCustom::operator <(dataWinsCustom other) const{
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
