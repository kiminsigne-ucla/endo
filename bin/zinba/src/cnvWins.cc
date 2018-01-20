#include "cnvWins.h"

cnvWins::~cnvWins(){}

cnvWins::cnvWins(){}

cnvWins::cnvWins(const cnvWins& s){
	start = s.start;
	stop = s.stop;
	cnvScore = s.cnvScore;
	varWin = s.varWin;
	chiSq = s.chiSq;
	pGap = s.pGap;
}

cnvWins::cnvWins(unsigned long int _start,unsigned long int _stop,double _cnvScore,double _varWin,double _chiSq,double _pGap){
	start = _start;
	stop = _stop;
	cnvScore = _cnvScore;
	varWin = _varWin;
	chiSq = _chiSq;
	pGap = _pGap;
}

bool cnvWins::operator <(cnvWins other) const{
	if( start < other.start ){
		return true;
	}else{
		return false;	
	}
}
