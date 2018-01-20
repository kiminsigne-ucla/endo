#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short	      /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int	      /* Wants to be signed 32 bits. */
#define boolean bool	      /* Wants to be signed 32 bits. */

#include <iostream>
#include <fstream>
#include <stdio.h>
extern "C"{
#include "twoBit.h"
}
#include "import.h"
#include "calcCovs.h"
//#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <Rmath.h>
#include <R.h>
#include <algorithm>
#include <time.h>
#include "string.h"

using namespace std;

/*
calcCovs::calcCovs(){
	chromCounter = 0;
	tbSizeFlag = 0;
}

calcCovs::~calcCovs(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}
*/

int calcCovs::processSignals(int zWinSize, int zOffsetSize, int cWinSize, int cOffsetSize, string alignDir,const char * twoBitFile,const char * inputFile,string outfile,const char * flist,int extension,const char * filetype, int binary){
	FILE * tempTB;
	time_t rtime;
	struct tm *timeinfo;
	char tInfo[128];// = "tempInfo.txt";
	char sysCall[256];

	unsigned short int currchr = 999;
	unsigned short int * basepair = NULL;
	unsigned short int * ibasepair = NULL;
	unsigned char * gcContent = NULL;	
	unsigned char * alignability = NULL;
	int i;
	int printflag = 0;

	int readInput = 0;
	const char* noneVal = "none";
	char mapChar[1];

	while(!b.signal_slist.empty()){
		i = 0;
		currchr = b.signal_slist[0].chrom;
		const char * chromReport = b.getKey(currchr);
		char * chromReport2=(char*) malloc(strlen(chromReport)+1) ; strcpy(chromReport2,chromReport);
		char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
		cout << "\nProcessing " << chromReport << endl;
		cout << "\tInitializing to length " << b.chr_size[currchr] << endl;
		basepair = new unsigned short int[b.chr_size[currchr]+1];
		for(int ch = b.chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(b.signal_slist[i].chrom==currchr && i < (int) b.signal_slist.size()){		
			basepair[b.signal_slist[i].pos]++;
			i++;
		}
		b.signal_slist.erase(b.signal_slist.begin(),b.signal_slist.begin()+i);

		alignability = new unsigned char[b.chr_size[currchr] + 1];
		for(int ch = b.chr_size[currchr]; ch--;)
			alignability[ch] = (unsigned char) 0;
			unsigned long int pos = 1;
		if(binary == 0){

			string alignFileS = alignDir + chromReport + ".wig";
			char * alignFile = new char[alignFileS.size() + 1];
			strcpy(alignFile, alignFileS.c_str());
			cout << "\tGetting alignability info from:\n\t\t" << alignFile << endl;

			tempTB = fopen(alignFile,"r");
			if(tempTB == NULL){
				cout << "Unable to open alignability file " << alignFile << endl;
				return 1;
			}
			unsigned short int aScore;
			while(!feof(tempTB)){
				int ret = fscanf(tempTB,"%hu",&aScore);
				alignability[pos] = (unsigned char) aScore;
				pos++;
			}
			fclose(tempTB);

		}else if(binary == 1){
			string alignFileS = alignDir + chromReport + ".bwig";
			char * alignFile = new char[alignFileS.size() + 1];
			strcpy(alignFile, alignFileS.c_str());
			cout << "\tGetting alignability info from:\n\t\t" << alignFile << endl;
			int intmapchr;
			ifstream mfile (alignFile,ios::in|ios::binary);
			if(mfile.is_open()){
				unsigned long int pos = 1;							
				
				while(!mfile.eof()){
					mfile.read((char *)&intmapchr,1);
					unsigned char mchr = mapChar[0];
					alignability[pos] = mchr;
					pos++;
				}
			}else{
				cout << "ERROR: Unable to open chromosome file " << alignFile << endl;
				return 1;
			}
			mfile.close();
			
		}else{
			cout << "ERROR: BINARY FLAG MUST BE PASS AS EITHER 0 OR 1 IN BUILDWINDOWS.R" << endl;
			exit(1);
		}

		
		
		cout << "\tGetting sequence from .2bit file:\n\t\t" << twoBitFile << endl;
		gcContent = new unsigned char[b.chr_size[currchr] + 1];
		for(int ch = b.chr_size[currchr]; ch--;)
			gcContent[ch] = (unsigned char) 0;

		char tSeq[128];
		time(&rtime);
		timeinfo=localtime(&rtime);
		//strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo); 
		strftime(tSeq,128,"tempSeq_%H_%M_%S.txt",timeinfo);  //need to add random number

    string tSeq2= outfile+tSeq;
		char *tSeq3 = new char[tSeq2.size() + 1];
		strcpy(tSeq3, tSeq2.c_str());   
		twoBitToFa2(chromReport2,1,b.chr_size[currchr],twoBitFile2,tSeq3);

//		twoBitToFa2(chromReport2,1,b.chr_size[currchr],twoBitFile2,tSeq);
		free(chromReport2);
		free(twoBitFile2);
		// opens R file tempTB = fopen(tInfo,"w");
		// prints command to R file fprintf(tempTB,"library(zinba);\ntwobittofa(chrm=\"%s\",start=1,end=%lu,twoBitFile=\"%s\",gcSeq=\"%s\");\n",chromReport,b.chr_size[currchr],twoBitFile,tSeq);
		// fclose (tempTB);
		// command to run r file in batch sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		/* this part runs the system command and retries if failure
		int s = 1;
		int twobitCount = 0;
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitToFa again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "twoBitToFa failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);
		*/

		ifstream seqfile(tSeq3);
		string line;
		pos = 1;
		unsigned long int nStart = 0;
		unsigned long int nStop = 0;
		unsigned short int prevEnt = 0;
		while(getline(seqfile,line)){
			if(line[0] != '>'){
				for(int s = 0; s < line.length();s++){
					if(line[s] == 'G' || line[s] == 'C'){
						gcContent[pos] = (unsigned char) 1;
					}else if (line[s] == 'N'){
						gcContent[pos] = (unsigned char) 2;
					}else{
						gcContent[pos] = (unsigned char) 0;
					}
					pos++;
				}
			}
		}
		seqfile.close();
		remove(tSeq3);

		cout << "\tGetting counts for " << cWinSize << "bp windows.........." << endl;
		int numOffsets = 1;
		if(cOffsetSize > 0){
			numOffsets = (int) cWinSize/cOffsetSize;
		}
		for(int o = 0; o < numOffsets; o++){
			unsigned long int cWinStart = (cOffsetSize * o) + 1;
			unsigned long int cWinStop = cWinStart + cWinSize - 1;
			if(cWinStop > b.chr_size[currchr])
				cWinStop = b.chr_size[currchr];
			while(cWinStop <= b.chr_size[currchr]){
				double cnvCount = 0.0;
				int alignCount = 0;
				double nCount = 0.0;
				for(int b = cWinStart; b <= cWinStop; b++){
					if((int) gcContent[b] == 2){
						nCount++;						
					}else{
						cnvCount += basepair[b];
						alignCount += (int) alignability[b];
					}
				}
				double cScore = 0.0;
				if(alignCount > 0)
					cScore = (double) cnvCount/alignCount;
				double percGap = (double) nCount/cWinSize;
				cnvWins cnv(cWinStart,cWinStop,cScore,0,0,percGap);
				cnv_wins.push_back(cnv);
				cWinStart += cWinSize;
				cWinStop += cWinSize;
			}
		}
			
		cout << "\t\tRefining boundaries...." << endl;
		cnv_wins.sort();
		
//const char * cWinFile = "cnv_wins.txt";
//tempTB = fopen(cWinFile,"w");
//list<cnvWins>::iterator cf = cnv_wins.begin();
//while(cf != cnv_wins.end()){
//	if(cf->pGap < 0.01){
//		long unsigned int cwpos = (long unsigned int) (cf->stop+cf->start)/2;
//		fprintf(tempTB,"%lu\t%f\n",cwpos,cf->cnvScore);
//	}
//	cf++;
//}
//fclose (tempTB);
				
		int numCnvWins = 0;
		double slideWinSize = numOffsets * 2.0;
		double globalSum = 0.0;
		double globalSumX2 = 0.0;
		double gapThresh = 0.0;
		int firstVar = cnv_wins.size();
		int countWin = 0;
		list<double> localSum;
		list<double> localSumX2;

//const char * varFile = "var_wins.txt";
//tempTB = fopen(varFile,"w");
		
		list<cnvWins>::iterator c = cnv_wins.begin();
		while(c != cnv_wins.end()){
			if(c->pGap <= gapThresh){
				globalSum += c->cnvScore;
				globalSumX2 += pow((c->cnvScore),2);
				numCnvWins++;
				if(localSum.size() < (slideWinSize - 1.0)){
					localSum.push_back(c->cnvScore);
					localSumX2.push_back(pow((c->cnvScore),2));
				}else{
					localSum.push_back(c->cnvScore);
					localSumX2.push_back(pow((c->cnvScore),2));
					if(localSum.size() == (slideWinSize + 1)){
						localSum.pop_front();
						localSumX2.pop_front();
					}

					list<double>::iterator x = localSum.begin();
					list<double>::iterator x2 = localSumX2.begin();
					double sum = 0.0;
					double sumX2 = 0.0;
					while(x != localSum.end()){
						sum += *x;
						sumX2 += *x2;
						x++;x2++;
					}
					c->varWin = (sumX2 -((sum * sum)/slideWinSize))/(slideWinSize-1);
//long unsigned int vpos = (unsigned long int) ((c->start+c->stop)/2)-cWinSize;
//fprintf(tempTB,"%lu\t%e\n",vpos,c->varWin);
					if(countWin < firstVar)
						firstVar = countWin;
				}
				c++;
			}else{
				if(localSum.size() > 0){
					localSum.clear();
					localSumX2.clear();
				}
				while(c->pGap != 0 && c!= cnv_wins.end())
					c++;
			}
			countWin++;
		}
		
//fclose (tempTB);

//const char * pvalFile = "pval_wins.txt";
//tempTB = fopen(pvalFile,"w");
		
		localSum.clear();
		localSumX2.clear();
		double globalVar = (globalSumX2 -((globalSum * globalSum)/numCnvWins))/(numCnvWins-1);
		double degfree = (slideWinSize - 1);
		cout << "\t\t\tGlobal variance is " << globalVar << endl;
		list<cnvWins> sigBoundary;
		c = cnv_wins.begin();
		int cnvpos = 0;
		while(cnvpos < firstVar){
			c++;
			cnvpos++;
		}
		while(c != cnv_wins.end()){
			if(c->pGap <= gapThresh){
				c->chiSq = (degfree*pow(c->varWin,2))/pow(globalVar,2);
				double cPval = dchisq(c->chiSq, degfree, 0);
//				c->pval = chisquarecdistribution(degfree,c->chiSq);
				
//long unsigned int vpos = (unsigned long int) ((c->start+c->stop)/2)-cWinSize;
//fprintf(tempTB,"%lu\t%e\n",vpos,cPval);
				
				if(cPval >= 0.000001){
					sigBoundary.push_back(*c);
				}
				c++;
			}else{
				c++;
			}
		}
		
//fclose (tempTB);
		
		sigBoundary.sort();
		list<cnvWins>::iterator sb = sigBoundary.begin();
		cnvWins lastCNV = *sb;
		sb++;
		while(sb != sigBoundary.end()){
			if(sb->start >= lastCNV.start && sb->start <= lastCNV.stop){
				long unsigned int lcStop = sb->stop;
				sigBoundary.erase(sb--);
				sb->stop = lcStop;
			}
			lastCNV = *sb;
			sb++;
		}
		cout << "\t\t\tRefining " << sigBoundary.size() << " boundaries" << endl;		
		sb = sigBoundary.begin();
		
//while(sb != sigBoundary.end()){
//	cout << "\t\t\tSIG BOUND is " << sb->start << " " << sb->stop << endl;
//	sb++;
//}
//sb = sigBoundary.begin();
		
		list<long unsigned int> transPts;
		long unsigned int refOffset = 2 * cWinSize;
		long unsigned int leftWinStart = sb->start;
		if(sb->start > refOffset)
			leftWinStart = sb->start - refOffset;

		while(!sigBoundary.empty()){
			
			int searchLen = cWinSize;
			leftWinStart = sb->start - refOffset;

			long unsigned int leftWinStop;
			if((b.chr_size[currchr] - leftWinStart) > refOffset){
				leftWinStop = leftWinStart + cWinSize;
			}else{
				searchLen = (int) (b.chr_size[currchr] - leftWinStart)/4;
				leftWinStop = leftWinStart + searchLen;
			}
			
			double lowPval = 1;
			double maxRatio = 0;
			long unsigned int transPoint = 0;
			
			while(leftWinStop <= (sb->stop-searchLen) && (leftWinStop+searchLen+1) <= b.chr_size[currchr]){
				int leftWinCount = 0;
				int rightWinCount = 0;

				for(int s = leftWinStart; s <= leftWinStop; s++){
					leftWinCount += basepair[s];
					rightWinCount += basepair[(s+searchLen+1)];
				}

				double maxCount = (double) leftWinCount;
				if(rightWinCount > leftWinCount)
					maxCount = (double) rightWinCount;
				double sumCount = (double) leftWinCount + rightWinCount;
				double ratioDiff = (double) maxCount/sumCount;
				double bPval = dbinom(maxCount,sumCount,0.5,0);
//				double bPval = binomialcdistribution(maxCount,sumCount,0.5);
				if(bPval <= lowPval){
					if(ratioDiff > maxRatio){
						lowPval = bPval;
						maxRatio = ratioDiff;
						transPoint = leftWinStop + 1;
					}
				}
				leftWinStart += zWinSize;
				leftWinStop += zWinSize;
//				leftWinStart += zOffsetSize;
//				leftWinStop += zOffsetSize;
			}
			transPts.push_back(transPoint);
			sigBoundary.erase(sb++);
		}
		transPts.sort();
		list<long unsigned int>::iterator tp = transPts.begin();
		c = cnv_wins.begin();
		while(c != cnv_wins.end() && tp != transPts.end()){
			if(c->start >= *tp)
				tp++;
			if(c->start < *tp && c->stop > *tp && tp != transPts.end())
				cnv_wins.erase(c++);
			else
				c++;
		}
		tp = transPts.begin();
		
////////////////////////////////////////////
//
// need to make sure that new cnv windows
// don't overlap other trans points
//
////////////////////////////////////////////
		
		while(!transPts.empty()){
//cout << "\t\t\t\tTransition point at " << *tp << endl;
			double leftCnvCount = 0;
			int leftAlignCount = 0;
			double rightCnvCount = 0;
			int rightAlignCount = 0;
			
			unsigned long int tpStart = 1;
			if(*tp > (cWinSize+1))
				tpStart = (*tp-cWinSize);
			unsigned long int tpStop = (unsigned long int) (b.chr_size[currchr]-tpStart)/2;
			if(b.chr_size[currchr] > (*tp+cWinSize)){
				tpStop = *tp;
			}
			for(int b = tpStart; b <= tpStop; b++){
				leftCnvCount += basepair[b];
				leftAlignCount += (int) alignability[b];
				rightCnvCount += basepair[b+cWinSize];
				rightAlignCount += (int) alignability[b+cWinSize];
			}
			double cScore = 0;
			if(leftAlignCount > 0)
				cScore = leftCnvCount/leftAlignCount;
			cnvWins cnv(tpStart,tpStop,cScore,0,0,0);
			cnv_wins.push_back(cnv);
			cScore = 0;
			if(rightAlignCount > 0)
				cScore = rightCnvCount/rightAlignCount;
			cnvWins cnvR((tpStart+cWinSize),(tpStop+cWinSize),cScore,0,0,0);
			cnv_wins.push_back(cnvR);
			transPts.erase(tp++);
		}
		cnv_wins.sort();

		if(strcmp(inputFile,noneVal)!=0){
			if(readInput == 0){
				cout << "\tLoading reads from input file " << inputFile << "........." << endl;
				int rVal = importRawSignal(inputFile,extension,filetype,1,twoBitFile);
				if(rVal == 1){
					cout << "Unable to open file with input reads" << endl;
					return 1;
				}
				readInput = 1;
			}
			cout << "\tMapping input tags to the genome........." << endl;
			i = 0;
			ibasepair = new unsigned short int[b.chr_size[currchr]+1];
			for(int ch = b.chr_size[currchr]; ch--;)
				ibasepair[ch] = 0;

			while(b.input_slist[i].chrom==currchr && i < (int) b.input_slist.size()){
				ibasepair[b.input_slist[i].pos]++;
				i++;
			}
			b.input_slist.erase(b.input_slist.begin(),b.input_slist.begin()+i);
		}
		
		cout << "\tGetting counts for zinba windows.........." << endl;
		if(printflag==0){
			tempTB = fopen(flist,"w");
			printflag = 1;
		}else{
			tempTB = fopen(flist,"a");
		}
		numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = (int) zWinSize/zOffsetSize;
		}
		string outfileDATA;
		slist<dataWins>::iterator z;
		for(int o = 0; o < numOffsets; o++){		
			z = peak_wins.previous(peak_wins.end());	
			list<cnvWins>::iterator cnvBegin = cnv_wins.begin();
			list<cnvWins>::iterator cnvEnd = cnv_wins.begin();
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			char offset[128];
			sprintf(offset,"%d",(zOffsetSize * o));
			char winsize[128];
			sprintf(winsize,"%d",zWinSize);			
			outfileDATA = outfile + "_" + chromReport + "_win" + string(winsize) + "bp_offset" + string(offset) + "bp.txt";
			
			if(o == (numOffsets-1))
				fprintf(tempTB,"%s\n",outfileDATA.c_str());
			else
				fprintf(tempTB,"%s;",outfileDATA.c_str());
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;
			while(zWinStop <= b.chr_size[currchr]){
				int peakCount = 0;
				double alignCount = 0;
				double gcCount = 0;
				double nCount = 0;
				double inCount = 0;
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];
					alignCount += (int) alignability[b];
					
					if((int) gcContent[b] == 1)
						gcCount++;
					else if((int) gcContent[b] == 2)
						nCount++;
					
					if(readInput == 1)
						inCount += ibasepair[b];
				}
				if((nCount/zWinSize) < 0.1){
					double cnvSum = 0;
					int cnvCount = 0;
					while(cnvBegin->stop < zWinStart && cnvBegin != cnv_wins.end())
						cnvBegin++;
					cnvEnd = cnvBegin;
					while( ((cnvEnd->start <= zWinStart && cnvEnd->stop >= zWinStart) || (cnvEnd->start <= zWinStop && cnvEnd->stop >= zWinStop)) && cnvEnd != cnv_wins.end() ){
						cnvSum += cnvEnd->cnvScore;
						cnvCount++;
						cnvEnd++;
					}
					
					double cnvLogScore = 0.0;
					if(cnvCount > 0)
						cnvLogScore = log(((cnvSum/cnvCount)*(zWinSize*2.0))+1.0);
					inCount = log((inCount+1));
					double gcPerc = gcCount/zWinSize;
					double aPerc = alignCount/(zWinSize*2.0);
					dataWins zwin(currchr,zWinStart,zWinStop,peakCount,inCount,gcPerc,aPerc,cnvLogScore);
					z = peak_wins.insert_after(z,zwin);
				}
				zWinStart += zWinSize;
				zWinStop += zWinSize;
			}
			
			if(outputData(outfileDATA.c_str(),currchr) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}else{
				cout << "\tPrinted out data to " << outfileDATA << endl;
			}
			
		}
		fclose(tempTB);
		cnv_wins.clear();
		delete [] basepair;
		basepair = NULL;
		delete [] gcContent;
		gcContent = NULL;
		delete [] alignability;
		alignability = NULL;
		if(readInput == 1){
			delete [] ibasepair;
			ibasepair = NULL;
		}
	}
	return 0;
}

int calcCovs::outputData(const char * outputFile, unsigned short int currChr){
	FILE * fh;
	fh = fopen(outputFile,"w");
	fprintf(fh,"chromosome\tstart\tstop\texp_count\tinput_count\tgcPerc\talign_perc\texp_cnvwin_log\n");
	const char * chrom = b.getKey(currChr);
	slist<dataWins>::iterator c = peak_wins.begin();
	while(c != peak_wins.end()){
		fprintf(fh,"%s\t%lu\t%lu\t%i\t%f\t%f\t%f\t%f\n",chrom,c->start,c->stop,c->eCount,c->iCount,c->gcPerc,c->alignPerc,c->cnvScore);
		peak_wins.erase(c++);
	}
	fclose (fh);
	return 0;
}






int calcCovs::importRawSignal(const char * signalFile,int extension,const char * filetype,int dataType,const char * twoBitFile){
	if(b.tbSizeFlag == 0){
		FILE * tempTB;
		time_t rtime;
		struct tm *timeinfo;

		char tInfo[128];// = "tempInfo.txt";
		char tChrSize[128];// = "tempChromSize.txt";
		char sysCall[256];
		char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
		time(&rtime);
		timeinfo=localtime(&rtime);
		//strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tChrSize,128,"tempChromSize_%H_%M_%S.txt",timeinfo);

		twoBitInfo2(twoBitFile2, tChrSize);
		free(twoBitFile2);
		/*tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobitinfo(infile=\"%s\",outfile=\"%s\");\n",twoBitFile,tChrSize);
		fclose (tempTB);

		cout << "\tGetting chromosome lengths from .2bit file: " << twoBitFile << endl;
		int s = 1;
		int twobitCount = 0;
		sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitInfo again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "TwoBitInfo failed, exiting" << endl;
				return 1;
			}
		}import
		remove(tInfo);*/

		tempTB = fopen(tChrSize,"r");
		char cChrom[128];
		unsigned long int cStart;
		while(!feof(tempTB)){
			int ret = fscanf(tempTB,"%s%lu",cChrom,&cStart);
			//cout << "\t\tFor " << cChrom << " length is " << cStart << endl;
			unsigned short int chromInt = b.getHashValue(cChrom);
			b.chr_size[chromInt] = cStart;
		}
		fclose(tempTB);
		remove(tChrSize);
		b.tbSizeFlag = 1;
	}

	const char * bowtie = "bowtie";
	const char * bed = "bed";
	const char * tagAlign = "tagAlign";
	int rval = 0;

	if(strcmp(filetype,bed) == 0){
		rval = b.importBed(signalFile, (int) (extension/2) ,dataType);
	}else if (strcmp(filetype,bowtie) == 0){
		rval = b.importBowtie(signalFile,(int) (extension/2),dataType);
	}else if(strcmp(filetype,tagAlign) == 0){
		rval = b.importTagAlign(signalFile,(int) (extension/2),dataType);
	}else{
		cout << "Unrecognized type of file " << filetype << ", must be either bowtie, bed, or tagAlign" << endl;
		return 1;
	}

	if(rval == 0){
		if(dataType == 0){
			cout << "\tImported " << b.signal_slist.size() << " reads" << endl;
			cout << "\tSorting reads ...";
			sort (b.signal_slist.begin(), b.signal_slist.end());
		}else if(dataType == 1){
			cout << "\t\tImported " << b.input_slist.size() << " inpput reads" << endl;
			cout << "\t\tSorting reads ...";
			sort (b.input_slist.begin(), b.input_slist.end());
		}
		cout << "COMPLETE" << endl;
	}else{
		cout << "Stopping buildwindows" << endl;
		return 1;
	}
	return 0;
}
/*
int calcCovs::b.importBowtie(const char * signalFile,int extension,int dataType){
	
	cout << "\tImporting bowtie formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
		
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	char line[512];char seq[128];
	char name[128];char sscore[128];int ival;
	int extend = (int)(extension/2);
	int rval;
	int num_skip = -1;
  int formatflag = 0 ;
	while(!feof(fh)){
		rval = fscanf(fh,"%s%s%s%lu%s%s%i",name,strand,cChrom,&pos,seq,sscore,&ival);
		if(rval == 8){
			formatflag = 1;
		}
		fgets(line,512,fh);
		if(rval == 7 | 8){
			if(strcmp(strand,minus) == 0){
				if((pos + strlen(seq)) >= extend)
					pos = (pos + strlen(seq)) - extend + 1;
				else
					pos = 1;
			}else{
				if((pos + extend-1) <= b.chr_size[getHashValue(cChrom)])
					pos += extend - 1;
				else
					pos = b.chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				b.signal_slist.push_back(sig);
			else if(dataType == 1)
				b.input_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}

	if(formatflag ==1 ){
		cout << "WARNING:  8 columns detected in bowtie file, ignoring last column pertaining to mismatch descriptiors (see Input Files in ZINBA tutorial on the ZINBA website)" << endl;
	}
	if((b.signal_slist.size()==0 & dataType==0) | (b.input_slist.size()==0 & dataType==1) ){
			cout << "error:  0 reads imported, check formatting of reads on zinba website" << endl;
			cout << "error: number of columns found is"<< rval << endl;
			return(1);
	}

	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}

int calcCovs::b.importTagAlign(const char * signalFile,int extension,int dataType){
	
	cout << "\tImporting tagAlign formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char seq[128];int score;
	int extend = (int)(extension/2);
	int rval;
	int num_skip = -1;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,seq,&score,strand);
		if(rval == 6){
			if(strcmp(strand,minus) == 0){
				if(stop >= extend)
					pos = stop - extend + 1;
				else
					pos = 1;
			}else{
				if((start + extend-1) <= b.chr_size[getHashValue(cChrom)])
					pos = start + extend - 1;
				else
					pos = b.chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				b.signal_slist.push_back(sig);
			else if(dataType == 1)
				b.input_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}

	if((b.signal_slist.size()==0 & dataType==0) | (b.input_slist.size()==0 & dataType==1) ){
			cout << "error:  0 reads imported, check formatting of reads on zinba website" << endl;
			cout << "error: number of columns found is"<< rval << endl;
			return(1);
	}

	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}

int calcCovs::importBed(const char * signalFile,int extension,int dataType){
	
	cout << "\tImporting bed formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[1];
	char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char name[128];int bscore;
	int extend = (int)(extension/2);
	int rval;
	int num_skip = -1;
	
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,name,&bscore,strand);
		if(rval == 6){
			if(strcmp(strand,minus) == 0){
				if(stop >= extend)
					pos = stop - extend + 1;
				else
					pos = 1;
			}else{
				if((start + extend-1) <= b.chr_size[getHashValue(cChrom)])
					pos = start + extend - 1;
				else
					pos = b.chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				b.signal_slist.push_back(sig);
			else if(dataType == 1)
				b.input_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}

	if((b.signal_slist.size()==0 & dataType==0) | (b.input_slist.size()==0 & dataType==1) ){
			cout << "error:  0 reads imported, check formatting of reads on zinba website" << endl;
			cout << "error: number of columns found is"<< rval << endl;
			return(1);
	}

	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}
*/

int calcCovs::outputDataWinCount(const char * outputFile, unsigned short int currChr){
	FILE * fh;
	fh = fopen(outputFile,"w");
	fprintf(fh,"chromosome\tstart\tstop\texp_count\n");
	const char * chrom = b.getKey(currChr);
	slist<dataWinsCount>::iterator c = peak_wins2.begin();
	while(c != peak_wins2.end()){
		fprintf(fh,"%s\t%lu\t%lu\t%i\n",chrom,c->start,c->stop,c->eCount);
		peak_wins2.erase(c++);
	}
	fclose (fh);
	return 0;
}

int calcCovs::processWinSignal(int zWinSize, int zOffsetSize,const char * twoBitFile,string outfile,int extension,const char * filetype, double Nthresh){
	time_t rtime;
	struct tm *timeinfo;
	FILE * tempTB;
	char tInfo[128];// = "tempInfo.txt";
	char sysCall[256];
	
	unsigned short int currchr = 999;
	unsigned int * basepair = NULL;
	unsigned int * ibasepair = NULL;
	unsigned char * gcContent = NULL;	
	//unsigned char * alignability = NULL;
	int i;
	int printflag = 0;

	int readInput = 0;
	const char* noneVal = "none";
	
	while(!b.signal_slist.empty()){
		i = 0;
		currchr = b.signal_slist[0].chrom;
		const char * chromReport = b.getKey(currchr);
		char * chromReport2=(char*) malloc(strlen(chromReport)+1) ; strcpy(chromReport2,chromReport);
		char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
		cout << "\nProcessing " << chromReport << endl;
		basepair = new unsigned int[b.chr_size[currchr]+1];
		for(int ch = b.chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(b.signal_slist[i].chrom==currchr && i < (int) b.signal_slist.size()){
			basepair[b.signal_slist[i].pos]++;
			i++;
		}	
		b.signal_slist.erase(b.signal_slist.begin(),b.signal_slist.begin()+i);

		cout << "\tGetting sequence from .2bit file:\n\t\t" << twoBitFile << endl;
		gcContent = new unsigned char[b.chr_size[currchr] + 1];
		for(int ch = b.chr_size[currchr]; ch--;)
			gcContent[ch] = (unsigned char) 0;

		char tSeq[128];
		time(&rtime);
		timeinfo=localtime(&rtime);
		//strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo);
		strftime(tSeq,128,"tempSeq_%H_%M_%S.txt",timeinfo);

		twoBitToFa2(chromReport2,1,b.chr_size[currchr],twoBitFile2,tSeq);
		free(chromReport2);
		free(twoBitFile2);
		/*tempTB = fopen(tInfo,"w");
		fprintf(tempTB,"library(zinba);\ntwobittofa(chrm=\"%s\",start=1,end=%lu,twoBitFile=\"%s\",gcSeq=\"%s\");\n",chromReport,b.chr_size[currchr],twoBitFile,tSeq);
		fclose (tempTB);
		sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		int s = 1;
		int twobitCount = 0;
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitToFa again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "twoBitToFa failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);*/
		
		ifstream seqfile(tSeq);
		string line;
		int pos = 1;
		unsigned long int nStart = 0;
		unsigned long int nStop = 0;
		unsigned short int prevEnt = 0;
		while(getline(seqfile,line)){
			if(line[0] != '>'){
				for(int s = 0; s < line.length();s++){
					if(line[s] == 'G' || line[s] == 'C'){
						gcContent[pos] = (unsigned char) 1;
					}else if (line[s] == 'N'){
						gcContent[pos] = (unsigned char) 2;
					}else{
						gcContent[pos] = (unsigned char) 0;
					}
					pos++;
				}
			}
		}
		seqfile.close();
		remove(tSeq);


		cout << "\tGetting counts for zinba windows.........." << endl;
		int numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = (int) zWinSize/zOffsetSize;
		}
		string outfileDATA;
		slist<dataWinsCount>::iterator z;
		for(int o = 0; o < numOffsets; o++){
			z = peak_wins2.previous(peak_wins2.end());
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			char offset[128];
			sprintf(offset,"%d",(zOffsetSize * o));
			char winsize[128];
			sprintf(winsize,"%d",zWinSize);
			outfileDATA = outfile + "_" + chromReport + "_win" + string(winsize) + "bp_offset" + string(offset) + "bp.txt";
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;
			while(zWinStop <= b.chr_size[currchr]){
				int peakCount = 0;
				double gcCount = 0;
				double nCount = 0;
				
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];
					if((int) gcContent[b] == 1)
						gcCount++;
					else if((int) gcContent[b] == 2)
						nCount++;
				}
				if((nCount/zWinSize) < Nthresh | Nthresh==0){
					dataWinsCount zwin(currchr,zWinStart,zWinStop,peakCount);
					z = peak_wins2.insert_after(z,zwin);
				}
				zWinStart += zWinSize;
				zWinStop += zWinSize;
			}
			
			if(outputDataWinCount(outfileDATA.c_str(),currchr) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}
			
		}
		delete [] basepair;
		basepair = NULL;
		}
	return 0;
}






int calcCovs::processCustomSignal(int zWinSize, int zOffsetSize,const char * twoBitFile,const char * inputFile,string outfile, const char* flist, int extension){
	FILE * tempTB;
	unsigned short int currchr = 999;
	unsigned short int * basepair = NULL;
	unsigned short int * ibasepair = NULL;
	int i;
	int printflag = 0;

	int readInput = 0;
	const char* noneVal = "none";
	
	while(!b.signal_slist.empty()){
		//import b;
		i = 0;
		currchr = b.signal_slist[0].chrom;
		const char * chromReport = b.getKey(currchr);
		cout << "\nProcessing " << chromReport << endl;
		cout << "\tInitializing to length " << b.chr_size[currchr] << endl;
		basepair = new unsigned short int[b.chr_size[currchr]+1];
		for(int ch = b.chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(b.signal_slist[i].chrom==currchr && i < (int) b.signal_slist.size()){		
			basepair[b.signal_slist[i].pos]++;
			i++;
		}
		b.signal_slist.erase(b.signal_slist.begin(),b.signal_slist.begin()+i);

		if(strcmp(inputFile,noneVal)!=0){
			if(readInput == 0){
				cout << "\tLoading reads from input file " << inputFile << "........." << endl;
				int rVal = b.importCustomBed(inputFile,extension);
				if(rVal == 1){
					cout << "Unable to open file with input reads" << endl;
					return 1;
				}
				readInput = 1;
			}
			cout << "\tMapping tags and scores from custom file to the genome........." << endl;
			i = 0;
			ibasepair = new unsigned short int[b.chr_size[currchr]+1];
			for(int ch = b.chr_size[currchr]; ch--;)
				ibasepair[ch] = 0;

			while(b.custom_slist[i].chrom==currchr && i < (int) b.custom_slist.size()){
				ibasepair[b.custom_slist[i].pos]++;
				i++;
			}
			b.custom_slist.erase(b.custom_slist.begin(),b.custom_slist.begin()+i);
		}
		
		cout << "\tGetting counts for zinba windows.........." << endl;
		if(printflag==0){
			tempTB = fopen(flist,"w");
			printflag = 1;
		}else{
			tempTB = fopen(flist,"a");
		}
		int numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = (int) zWinSize/zOffsetSize;
		}
		string outfileDATA;
		slist<dataWinsCustom>::iterator z;
		for(int o = 0; o < numOffsets; o++){		
			z = peak_wins3.previous(peak_wins3.end());	
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			char offset[128];
			sprintf(offset,"%d",(zOffsetSize * o));
			char winsize[128];
			sprintf(winsize,"%d",zWinSize);			
			outfileDATA = outfile + "_" + chromReport + "_win" + string(winsize) + "bp_offset" + string(offset) + "bp.txt";
			
			if(o == (numOffsets-1))
				fprintf(tempTB,"%s\n",outfileDATA.c_str());
			else
				fprintf(tempTB,"%s;",outfileDATA.c_str());
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;
			while(zWinStop <= b.chr_size[currchr]){
				int peakCount = 0;
				double inCount = 0;
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];					
					if(readInput == 1)
						inCount += ibasepair[b];
				}					
					inCount = log((inCount+1)); //need to check this
					dataWinsCustom zwin(currchr,zWinStart,zWinStop,peakCount,inCount);
					z = peak_wins3.insert_after(z,zwin);
					zWinStart += zWinSize;
					zWinStop += zWinSize;
			}

			}
			
			if(outputCustomData(outfileDATA.c_str(),currchr) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}else{
				cout << "\tPrinted out data to " << outfileDATA << endl;
			}
			
		}
		fclose(tempTB);
		cnv_wins.clear();
		delete [] basepair;
		basepair = NULL;
		if(readInput == 1){
			delete [] ibasepair;
			ibasepair = NULL;
		}
	return 0;	
}
	


int calcCovs::outputCustomData(const char * outputFile, unsigned short int currChr){
	FILE * fh;
	fh = fopen(outputFile,"w");
	fprintf(fh,"chromosome\tstart\tstop\texp_count\tinput_count\n");
	const char * chrom = b.getKey(currChr);
	slist<dataWinsCustom>::iterator c = peak_wins3.begin();
	while(c != peak_wins.end()){
		fprintf(fh,"%s\t%lu\t%lu\t%i\t%f\n",chrom,c->start,c->stop,c->eCount,c->iCount);
		peak_wins3.erase(c++);
	}
	fclose (fh);
	return 0;
}



// fast functions
int calcCovs::processSignalsfast(int zWinSize, int zOffsetSize, int cWinSize, int cOffsetSize, const char * twoBitFile, const char * inputFile,string outfile, const char * flist,int extension, const char * filetype){

  FILE * tempTB;
	time_t rtime;
	struct tm *timeinfo;
	char tInfo[128];// = "tempInfo.txt";
	char sysCall[256];

	unsigned short int currchr = 999;
	unsigned short int * basepair = NULL;
	unsigned short int * ibasepair = NULL;
	unsigned char * gcContent = NULL;	
	//unsigned char * alignability = NULL;
	int i;
	int printflag = 0;

	int readInput = 0;
	const char* noneVal = "none";
	char mapChar[1];

	while(!b.signal_slist.empty()){
		i = 0;
		currchr = b.signal_slist[0].chrom;
		const char * chromReport = b.getKey(currchr);
		char * chromReport2=(char*) malloc(strlen(chromReport)+1) ; strcpy(chromReport2,chromReport);
		char * twoBitFile2=(char*) malloc(strlen(twoBitFile)+1); strcpy(twoBitFile2,twoBitFile);
		cout << "\nProcessing " << chromReport << endl;
		cout << "\tInitializing to length " << b.chr_size[currchr] << endl;
		basepair = new unsigned short int[b.chr_size[currchr]+1];
		for(int ch = b.chr_size[currchr]; ch--;)
			basepair[ch] = 0;

		cout << "\tMapping reads to chromosome......" << endl;
		while(b.signal_slist[i].chrom==currchr && i < (int) b.signal_slist.size()){		
			basepair[b.signal_slist[i].pos]++;
			i++;
		}
		b.signal_slist.erase(b.signal_slist.begin(),b.signal_slist.begin()+i);

    /*
		alignability = new unsigned char[b.chr_size[currchr] + 1];
		for(int ch = b.chr_size[currchr]; ch--;)
			alignability[ch] = (unsigned char) 0;
			unsigned long int pos = 1;
			string alignFileS = alignDir + chromReport + ".wig";
			char * alignFile = new char[alignFileS.size() + 1];
			strcpy(alignFile, alignFileS.c_str());
			cout << "\tGetting alignability info from:\n\t\t" << alignFile << endl;

			tempTB = fopen(alignFile,"r");
			if(tempTB == NULL){
				cout << "Unable to open alignability file " << alignFile << endl;
				return 1;
			}
			unsigned short int aScore;
			while(!feof(tempTB)){
				int ret = fscanf(tempTB,"%hu",&aScore);
				alignability[pos] = (unsigned char) aScore;
				pos++;
			}
			fclose(tempTB);
    */
		
		
		cout << "\tGetting sequence from .2bit file:\n\t\t" << twoBitFile << endl;
		gcContent = new unsigned char[b.chr_size[currchr] + 1];
		for(int ch = b.chr_size[currchr]; ch--;)
			gcContent[ch] = (unsigned char) 0;

		char tSeq[128];
		time(&rtime);
		timeinfo=localtime(&rtime);
		//strftime(tInfo,128,"tempInfo_%H_%M_%S.txt",timeinfo); 
		strftime(tSeq,128,"tempSeq_%H_%M_%S.txt",timeinfo);  //need to add random number

    string tSeq2= outfile+tSeq;
		char *tSeq3 = new char[tSeq2.size() + 1];
		strcpy(tSeq3, tSeq2.c_str());   
		twoBitToFa2(chromReport2,1,b.chr_size[currchr],twoBitFile2,tSeq3);

//		twoBitToFa2(chromReport2,1,b.chr_size[currchr],twoBitFile2,tSeq);
		free(chromReport2);
		free(twoBitFile2);
		// opens R file tempTB = fopen(tInfo,"w");
		// prints command to R file fprintf(tempTB,"library(zinba);\ntwobittofa(chrm=\"%s\",start=1,end=%lu,twoBitFile=\"%s\",gcSeq=\"%s\");\n",chromReport,b.chr_size[currchr],twoBitFile,tSeq);
		// fclose (tempTB);
		// command to run r file in batch sprintf(sysCall,"R CMD BATCH %s /dev/null",tInfo);
		/* this part runs the system command and retries if failure
		int s = 1;
		int twobitCount = 0;
		while(s != 0){
			s = system(sysCall);
			twobitCount++;
			if(twobitCount < 5 && s != 0){
				cout << "Trying twoBitToFa again, s is" << s << endl;
			}else if(twobitCount >= 5 && s != 0){
				cout << "twoBitToFa failed, exiting" << endl;
				return 1;
			}
		}
		remove(tInfo);
		*/

		ifstream seqfile(tSeq3);
		string line;
		int pos = 1;
		unsigned long int nStart = 0;
		unsigned long int nStop = 0;
		unsigned short int prevEnt = 0;
		while(getline(seqfile,line)){
			if(line[0] != '>'){
				for(int s = 0; s < line.length();s++){
					if(line[s] == 'G' || line[s] == 'C'){
						gcContent[pos] = (unsigned char) 1;
					}else if (line[s] == 'N'){
						gcContent[pos] = (unsigned char) 2;
					}else{
						gcContent[pos] = (unsigned char) 0;
					}
					pos++;
				}
			}
		}
		seqfile.close();
		remove(tSeq3);

		cout << "\tGetting counts for " << cWinSize << "bp windows.........." << endl;
		int numOffsets = 1;
		if(cOffsetSize > 0){
			numOffsets = (int) cWinSize/cOffsetSize;
		}
		for(int o = 0; o < numOffsets; o++){
			unsigned long int cWinStart = (cOffsetSize * o) + 1;
			unsigned long int cWinStop = cWinStart + cWinSize - 1;
			if(cWinStop > b.chr_size[currchr])
				cWinStop = b.chr_size[currchr];
			while(cWinStop <= b.chr_size[currchr]){
				double cnvCount = 0.0;
				int alignCount = 1; //0;
				double nCount = 0.0;
				for(int b = cWinStart; b <= cWinStop; b++){
					if((int) gcContent[b] == 2){
						nCount++;						
					}else{
						cnvCount += basepair[b];
						//alignCount += (int) alignability[b];
					}
				}
				double cScore = 0.0;
				if(alignCount > 0)
					cScore = (double) cnvCount/alignCount;
				double percGap = (double) nCount/cWinSize;
				cnvWins cnv(cWinStart,cWinStop,cScore,0,0,percGap);
				cnv_wins.push_back(cnv);
				cWinStart += cWinSize;
				cWinStop += cWinSize; 
			}
		}
			
		//cout << "\t\tRefining boundaries...." << endl;
		//cnv_wins.sort();
		
//const char * cWinFile = "cnv_wins.txt";
//tempTB = fopen(cWinFile,"w");
//list<cnvWins>::iterator cf = cnv_wins.begin();
//while(cf != cnv_wins.end()){
//	if(cf->pGap < 0.01){
//		long unsigned int cwpos = (long unsigned int) (cf->stop+cf->start)/2;
//		fprintf(tempTB,"%lu\t%f\n",cwpos,cf->cnvScore);
//	}
//	cf++;
//}
//fclose (tempTB);
				
		int numCnvWins = 0;
		double slideWinSize = numOffsets * 2.0;
		double globalSum = 0.0;
		double globalSumX2 = 0.0;
		double gapThresh = 0.0;
		int firstVar = cnv_wins.size();
		int countWin = 0;
		list<double> localSum;
		list<double> localSumX2;

//const char * varFile = "var_wins.txt";
//tempTB = fopen(varFile,"w");
		
		list<cnvWins>::iterator c = cnv_wins.begin();
		while(c != cnv_wins.end()){
			if(c->pGap <= gapThresh){
				globalSum += c->cnvScore;
				globalSumX2 += pow((c->cnvScore),2);
				numCnvWins++;
				if(localSum.size() < (slideWinSize - 1.0)){
					localSum.push_back(c->cnvScore);
					localSumX2.push_back(pow((c->cnvScore),2));
				}else{
					localSum.push_back(c->cnvScore);
					localSumX2.push_back(pow((c->cnvScore),2));
					if(localSum.size() == (slideWinSize + 1)){
						localSum.pop_front();
						localSumX2.pop_front();
					}

					list<double>::iterator x = localSum.begin();
					list<double>::iterator x2 = localSumX2.begin();
					double sum = 0.0;
					double sumX2 = 0.0;
					while(x != localSum.end()){
						sum += *x;
						sumX2 += *x2;
						x++;x2++;
					}
					c->varWin = (sumX2 -((sum * sum)/slideWinSize))/(slideWinSize-1);
//long unsigned int vpos = (unsigned long int) ((c->start+c->stop)/2)-cWinSize;
//fprintf(tempTB,"%lu\t%e\n",vpos,c->varWin);
					if(countWin < firstVar)
						firstVar = countWin;
				}
				c++;
			}else{
				if(localSum.size() > 0){
					localSum.clear();
					localSumX2.clear();
				}
				while(c->pGap != 0 && c!= cnv_wins.end())
					c++;
			}
			countWin++;
		}
		
//fclose (tempTB);

//const char * pvalFile = "pval_wins.txt";
//tempTB = fopen(pvalFile,"w");
		
		localSum.clear();
		localSumX2.clear();
		double globalVar = (globalSumX2 -((globalSum * globalSum)/numCnvWins))/(numCnvWins-1);
		double degfree = (slideWinSize - 1);
		cout << "\t\t\tGlobal variance is " << globalVar << endl;
		list<cnvWins> sigBoundary;
		c = cnv_wins.begin();
		int cnvpos = 0;
		while(cnvpos < firstVar){
			c++;
			cnvpos++;
		}
		while(c != cnv_wins.end()){
			if(c->pGap <= gapThresh){
				c->chiSq = (degfree*pow(c->varWin,2))/pow(globalVar,2);
				double cPval = dchisq(c->chiSq, degfree, 0);
//				c->pval = chisquarecdistribution(degfree,c->chiSq);
				
//long unsigned int vpos = (unsigned long int) ((c->start+c->stop)/2)-cWinSize;
//fprintf(tempTB,"%lu\t%e\n",vpos,cPval);
				
				if(cPval >= 0.000001){
					sigBoundary.push_back(*c);
				}
				c++;
			}else{
				c++;
			}
		}
		
//fclose (tempTB);
		
		sigBoundary.sort();
		list<cnvWins>::iterator sb = sigBoundary.begin();
		cnvWins lastCNV = *sb;
		sb++;
		while(sb != sigBoundary.end()){
			if(sb->start >= lastCNV.start && sb->start <= lastCNV.stop){
				long unsigned int lcStop = sb->stop;
				sigBoundary.erase(sb--);
				sb->stop = lcStop;
			}
			lastCNV = *sb;
			sb++;
		}
		cout << "\t\t\tRefining " << sigBoundary.size() << " boundaries" << endl;		
		sb = sigBoundary.begin();
		
//while(sb != sigBoundary.end()){
//	cout << "\t\t\tSIG BOUND is " << sb->start << " " << sb->stop << endl;
//	sb++;
//}
//sb = sigBoundary.begin();
		
		list<long unsigned int> transPts;
		long unsigned int refOffset = 2 * cWinSize;
		long unsigned int leftWinStart = sb->start;
		if(sb->start > refOffset)
			leftWinStart = sb->start - refOffset;

		while(!sigBoundary.empty()){
			
			int searchLen = cWinSize;
			leftWinStart = sb->start - refOffset;

			long unsigned int leftWinStop;
			if((b.chr_size[currchr] - leftWinStart) > refOffset){
				leftWinStop = leftWinStart + cWinSize;
			}else{
				searchLen = (int) (b.chr_size[currchr] - leftWinStart)/4;
				leftWinStop = leftWinStart + searchLen;
			}
			
			double lowPval = 1;
			double maxRatio = 0;
			long unsigned int transPoint = 0;
			
			while(leftWinStop <= (sb->stop-searchLen) && (leftWinStop+searchLen+1) <= b.chr_size[currchr]){
				int leftWinCount = 0;
				int rightWinCount = 0;

				for(int s = leftWinStart; s <= leftWinStop; s++){
					leftWinCount += basepair[s];
					rightWinCount += basepair[(s+searchLen+1)];
				}

				double maxCount = (double) leftWinCount;
				if(rightWinCount > leftWinCount)
					maxCount = (double) rightWinCount;
				double sumCount = (double) leftWinCount + rightWinCount;
				double ratioDiff = (double) maxCount/sumCount;
				double bPval = dbinom(maxCount,sumCount,0.5,0);
//				double bPval = binomialcdistribution(maxCount,sumCount,0.5);
				if(bPval <= lowPval){
					if(ratioDiff > maxRatio){
						lowPval = bPval;
						maxRatio = ratioDiff;
						transPoint = leftWinStop + 1;
					}
				}
				leftWinStart += zWinSize;
				leftWinStop += zWinSize;
//				leftWinStart += zOffsetSize;
//				leftWinStop += zOffsetSize;
			}
			transPts.push_back(transPoint);
			sigBoundary.erase(sb++);
		}
		transPts.sort();
		list<long unsigned int>::iterator tp = transPts.begin();
		c = cnv_wins.begin();
		while(c != cnv_wins.end() && tp != transPts.end()){
			if(c->start >= *tp)
				tp++;
			if(c->start < *tp && c->stop > *tp && tp != transPts.end())
				cnv_wins.erase(c++);
			else
				c++;
		}
		tp = transPts.begin();
		
////////////////////////////////////////////
//
// need to make sure that new cnv windows
// don't overlap other trans points
//
////////////////////////////////////////////
		
		while(!transPts.empty()){
//cout << "\t\t\t\tTransition point at " << *tp << endl;
			double leftCnvCount = 0;
			int leftAlignCount = 1;//0;
			double rightCnvCount = 0;
			int rightAlignCount = 1;//0;
			
			unsigned long int tpStart = 1;
			if(*tp > (cWinSize+1))
				tpStart = (*tp-cWinSize);
			unsigned long int tpStop = (unsigned long int) (b.chr_size[currchr]-tpStart)/2;
			if(b.chr_size[currchr] > (*tp+cWinSize)){
				tpStop = *tp;
			}
			for(int b = tpStart; b <= tpStop; b++){
				leftCnvCount += basepair[b];
				//leftAlignCount += (int) alignability[b];
				rightCnvCount += basepair[b+cWinSize];
				//rightAlignCount += (int) alignability[b+cWinSize];
			}
			double cScore = 0;
			if(leftAlignCount > 0)
				cScore = leftCnvCount/leftAlignCount;
			cnvWins cnv(tpStart,tpStop,cScore,0,0,0);
			cnv_wins.push_back(cnv);
			cScore = 0;
			if(rightAlignCount > 0)
				cScore = rightCnvCount/rightAlignCount;
			cnvWins cnvR((tpStart+cWinSize),(tpStop+cWinSize),cScore,0,0,0);
			cnv_wins.push_back(cnvR);
			transPts.erase(tp++);
		}
		cnv_wins.sort();

		if(strcmp(inputFile,noneVal)!=0){
			if(readInput == 0){
				cout << "\tLoading reads from input file " << inputFile << "........." << endl;
				int rVal = importRawSignal(inputFile,extension,filetype,1,twoBitFile);
				if(rVal == 1){
					cout << "Unable to open file with input reads" << endl;
					return 1;
				}
				readInput = 1;
			}
			cout << "\tMapping input tags to the genome........." << endl;
			i = 0;
			ibasepair = new unsigned short int[b.chr_size[currchr]+1];
			for(int ch = b.chr_size[currchr]; ch--;)
				ibasepair[ch] = 0;

			while(b.input_slist[i].chrom==currchr && i < (int) b.input_slist.size()){
				ibasepair[b.input_slist[i].pos]++;
				i++;
			}
			b.input_slist.erase(b.input_slist.begin(),b.input_slist.begin()+i);
		}
		
		cout << "\tGetting counts for zinba windows.........." << endl;
		if(printflag==0){
			tempTB = fopen(flist,"w");
			printflag = 1;
		}else{
			tempTB = fopen(flist,"a");
		}
		numOffsets = 1;
		if(zOffsetSize > 0){
			numOffsets = (int) zWinSize/zOffsetSize;
		}
		string outfileDATA;
		slist<dataWins>::iterator z;
		for(int o = 0; o < numOffsets; o++){		
			z = peak_wins.previous(peak_wins.end());	
			list<cnvWins>::iterator cnvBegin = cnv_wins.begin();
			list<cnvWins>::iterator cnvEnd = cnv_wins.begin();
			cout << "\t\tOffset " << (zOffsetSize * o) << "bp......" << endl;
			char offset[128];
			sprintf(offset,"%d",(zOffsetSize * o));
			char winsize[128];
			sprintf(winsize,"%d",zWinSize);			
			outfileDATA = outfile + "_" + chromReport + "_win" + string(winsize) + "bp_offset" + string(offset) + "bp.txt";
			
			if(o == (numOffsets-1))
				fprintf(tempTB,"%s\n",outfileDATA.c_str());
			else
				fprintf(tempTB,"%s;",outfileDATA.c_str());
			unsigned long int zWinStart = (zOffsetSize * o) + 1;
			unsigned long int zWinStop = zWinStart + zWinSize - 1;
			while(zWinStop <= b.chr_size[currchr]){
				int peakCount = 0;
				double alignCount = 0;
				double gcCount = 0;
				double nCount = 0;
				double inCount = 0;
				for(int b = zWinStart; b <= zWinStop; b++){
					peakCount += basepair[b];
					//alignCount += (int) alignability[b];
					
					if((int) gcContent[b] == 1)
						gcCount++;
					else if((int) gcContent[b] == 2)
						nCount++;
					
					if(readInput == 1)
						inCount += ibasepair[b];
				}
				if((nCount/zWinSize) < 0.1){
					double cnvSum = 0;
					int cnvCount = 0;
					while(cnvBegin->stop < zWinStart && cnvBegin != cnv_wins.end())
						cnvBegin++;
					cnvEnd = cnvBegin;
					while( ((cnvEnd->start <= zWinStart && cnvEnd->stop >= zWinStart) || (cnvEnd->start <= zWinStop && cnvEnd->stop >= zWinStop)) && cnvEnd != cnv_wins.end() ){
						cnvSum += cnvEnd->cnvScore;
						cnvCount++;
						cnvEnd++;
					}
					
					double cnvLogScore = 0.0;
					if(cnvCount > 0)
						cnvLogScore = log(((cnvSum/cnvCount)*(zWinSize*2.0))+1.0);
					inCount = log((inCount+1));
					double gcPerc = gcCount/zWinSize;
					double aPerc = alignCount/(zWinSize*2.0);
					dataWins zwin(currchr,zWinStart,zWinStop,peakCount,inCount,gcPerc,aPerc,cnvLogScore);
					z = peak_wins.insert_after(z,zwin);
				}
				zWinStart += zWinSize;
				zWinStop += zWinSize;
			}
			
			if(outputData(outfileDATA.c_str(),currchr) != 0){
				cout << "Error printing output to file, exiting" << endl;
				exit(1);
			}else{
				cout << "\tPrinted out data to " << outfileDATA << endl;
			}
			
		}
		fclose(tempTB);
		cnv_wins.clear();
		delete [] basepair;
		basepair = NULL;
		delete [] gcContent;
		gcContent = NULL;
		//delete [] alignability;
		//alignability = NULL;
		if(readInput == 1){
			delete [] ibasepair;
			ibasepair = NULL;
		}
	}
	return 0;
}






