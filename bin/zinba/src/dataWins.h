#ifndef DATAWINS_H_
#define DATAWINS_H_

#include <iostream>

class dataWins {
	
	public:
			dataWins();
			~dataWins();
			dataWins(unsigned short int, unsigned long int, unsigned long int,int,double,double,double,double);
			dataWins(const dataWins&);
			bool operator <(dataWins) const;
			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int stop;
			int eCount;
			double iCount;
			double gcPerc;
			double alignPerc;
			double cnvScore;
	
	private:
			
};

#ifndef DATAWINSC_H_
#define DATAWINSC_H_

#include <iostream>
class dataWinsCount {
	
	public:
			dataWinsCount();
			~dataWinsCount();
			dataWinsCount(unsigned short int, unsigned long int, unsigned long int,int);
			dataWinsCount(const dataWinsCount&);
			bool operator <(dataWinsCount) const;
			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int stop;
			int eCount;
	private:
			
};

#ifndef DATAWINSCu_H_
#define DATAWINSCu_H_

#include <iostream>
class dataWinsCustom {
	
	public:
			dataWinsCustom();
			~dataWinsCustom();
			dataWinsCustom(unsigned short int, unsigned long int, unsigned long int,int,double);
			dataWinsCustom(const dataWinsCustom&);
			bool operator <(dataWinsCustom) const;
			unsigned short int chrom;//2
			unsigned long int start;//8
			unsigned long int stop;
			int eCount;
			double iCount;
	private:
			
};

#endif /*DATAWINSCu_H_*/

#endif /*DATAWINSC_H_*/

#endif /*DATAWINS_H_*/

