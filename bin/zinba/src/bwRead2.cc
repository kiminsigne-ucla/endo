#include "bwRead2.h"

bwRead2::bwRead2(){}

bwRead2::~bwRead2(){}

bwRead2::bwRead2(const bwRead2& s){
	chrom = s.chrom;
	pos = s.pos;
	score =s.score;
//	strand = s.strand;
}

bwRead2::bwRead2(unsigned short int _chrom, unsigned long int _pos, int _score){
	chrom = _chrom;
	pos = _pos;
	score = _score;
//	strand = _strand;
}


bool bwRead2::operator <(bwRead2 other) const{
		return (chrom<other.chrom);
}

