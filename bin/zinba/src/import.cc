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
//#include "calcCovs.h"
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
#include "import.h"


#define min(a,b)  ((a) < (b) ? (a) : (b))
#define maxx(x,y) ((x) > (y) ? (x) : (y))
#define MAX_LEN 1000000

import::import(){
	chromCounter = 0;
	tbSizeFlag = 0;
}

import::~import(){
	map<const char*, int>::iterator i;
	for(i=chroms.begin();i!=chroms.end();++i){
		delete [] i->first;	
	}
}

int import::importBowtie(const char * signalFile,int extension,int dataType){
	cout << "\tImporting bowtie formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
		
	char cChrom[128];
	unsigned long int pos;
	char strand[8];
	//char minus[] = "-";
	char line[512];char seq[128];
	char name[128];char sscore[128];int ival;	
	//int extend = (int)(extension/2); implement this above
	char *delim = "\t";
	char *ptr= NULL;  
  char str_buf[MAX_LEN];
	int num_skip = 0;
  int formatFlag = 0 ;
	int m=0; 
	int i;

	//rval = fscanf(fh,"%s%s%s%lu%s%s%i",name,strand,cChrom,&pos,seq,sscore,&ival);
 	while(fgets(str_buf, MAX_LEN, fh) != NULL){
   i = 0;
   if ((ptr = strtok(str_buf, delim)) != NULL) {
   do {
     i++;
     if(i==1){
       //strcpy(name, ptr);
     }else if(i==2){
       strcpy(strand, ptr);
     }else if(i==3){
       strcpy(cChrom, ptr);
     }else if(i==4){
       pos = atol(ptr);
     }else if(i==5){
       strcpy(seq, ptr);
     }else if(i==6){
			//strcpy(sscore, ptr);
     }else if(i==7){
			//ival = atoi(ptr);
     }
		//anything greater is ignored 
   } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", signalFile);
  }
	if(i==7 | i==8){
			if(strcmp(strand,"-") == 0 | strcmp(strand,"-\n") == 0){
				if((pos + strlen(seq)) >= extension)
					pos = (pos + strlen(seq)) - extension + 1;
				else
					pos = 1;
			}else{
				if((pos + extension-1) <= chr_size[getHashValue(cChrom)])
					pos += extension - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				signal_slist.push_back(sig);
			else if(dataType == 1)
				input_slist.push_back(sig);
		}else{
			formatFlag=1;
			num_skip++;
		}
	}
	
	if(i != 7){
		if(i<7){
			cout << "WARNING:  Less than 7 bowtie columns detected (see ZINBA wiki for import file format)" << endl;
		}else if(i==8){ 	
			cout << "Note:  8 bowtie columns detected, ignoring 8th column" << endl;
		}else{
			cout << "WARNING:  More than 8 bowtie columns detected (see ZINBA wiki for import file format)" << endl;
		} 
	}else{
		if(signal_slist.size() <=1 & dataType==0){
				cout << "ERROR:  0 experimental reads imported, check formatting of reads" << endl;
				cout << "ERROR: number of columns found is"<< i << endl;
				return(1);
		}
		if(input_slist.size()<=1 & dataType==1) {
				cout << "ERROR:  0 input reads imported, check formatting of readse" << endl;
				cout << "ERROR: number of columns found is"<< i << endl;
				return(1);
		}
		fclose(fh);
		cout << "\tSkipped " << num_skip << " reads" << endl;
		return 0;
	}
}


int import::importTagAlign(const char * signalFile,int extension,int dataType){
	
	cout << "\tImporting tagAlign formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	
	char cChrom[128];
	unsigned long int pos;
	char strand[8];
	//char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char seq[128];int score;
	//int extend = (int)(extension/2);
	char *delim = "\t";
	char *ptr= NULL;  
  char str_buf[MAX_LEN];
	int num_skip = 0;
  int formatFlag = 0 ;
	int m=0; 
	int i;

	//rval = fscanf(fh,"%s%lu%lu%s%i%s",cChrom,&start,&stop,seq,&score,strand);

 	while(fgets(str_buf, MAX_LEN, fh) != NULL){
   i = 0;
   if ((ptr = strtok(str_buf, delim)) != NULL) {
   do {
     i++;
     if(i==1){
       strcpy(cChrom, ptr);
     }else if(i==2){
       start = atol(ptr);
     }else if(i==3){
       stop = atol(ptr);
     }else if(i==4){
      //strcpy(seq, ptr);
     }else if(i==5){
			//score = atoi(ptr);
     }else if(i==6){
			strcpy(strand, ptr);
     }
		//anything greater is ignored 
   } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", signalFile);
  }

	if(i==6){
			if(strcmp(strand,"-") == 0 | strcmp(strand,"-\n") == 0){
				if(stop >= extension)
					pos = stop - extension + 1;
				else
					pos = 1;
			}else{
				if((start + extension-1) <= chr_size[getHashValue(cChrom)])
					pos = start + extension - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				signal_slist.push_back(sig);
			else if(dataType == 1)
				input_slist.push_back(sig);
		}else{
			formatFlag=1;
			num_skip++;
		}
		m++;
	}

	if(i<6){
		cout << "WARNING:  Less than 6 tagAlign columns detected.  ZINBA expects exactly 6 tagAlign columns (see ZINBA wiki for import file format)" << endl;
	}else if(i>6){ 	
			cout << "WARNING:  More than 6 tagAlign columns detected.  ZINBA expects exactly 6 tagAlign columns (see ZINBA wiki for import file format)" << endl;
	}
	
	if(signal_slist.size() <=1 & dataType==0){
		cout << "ERROR: 0 experimental reads imported, check formatting of reads" << endl;
		cout << "ERROR: number of columns found is"<< i << endl;
		return(1);
	}
	if(input_slist.size()<=1 & dataType==1) {
		cout << "ERROR: 0 input reads imported, check formatting of reads" << endl;
		cout << "ERROR: number of columns found is"<< i << endl;
		return(1);
	}

	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}


int import::importBed(const char * signalFile,int extension,int dataType){
	
	cout << "\tImporting bed formatted reads" << endl;
	FILE * fh;
	fh = fopen(signalFile,"r");
	if(fh == NULL){
		cout << "ERROR: Unable to open file containing reads " << signalFile << endl;
		return 1;
	}
	char cChrom[128];
	unsigned long int pos;
	char strand[8];
	//char minus[] = "-";
	unsigned long int start;unsigned long int stop;
	char name[128];int bscore;
	//int extend = (int)(extension/2);
	char *delim = "\t";
	char *ptr= NULL;  
  char str_buf[MAX_LEN];
	int num_skip = 0;
  int formatFlag = 0 ;
	int m=0; 
	int i;
	
 	while(fgets(str_buf, MAX_LEN, fh) != NULL){
   i = 0;
   if ((ptr = strtok(str_buf, delim)) != NULL) {
   do {
     i++;
     if(i==1){
       strcpy(cChrom, ptr);
     }else if(i==2){
       start = atol(ptr);
     }else if(i==3){
       stop = atol(ptr);
     }else if(i==4){
       //strcpy(name, ptr);
     }else if(i==5){
			//bscore = atoi(ptr);
     }else if(i==6){
			strcpy(strand, ptr);
     }
		//anything greater is ignored 
   } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", signalFile);
  }
	if(i==6){
			if(strcmp(strand,"-") == 0 | strcmp(strand,"-\n") == 0){
				if(stop >= extension)
					pos = stop - extension + 1;
				else
					pos = 1;
			}else{
				if((start + extension-1) <= chr_size[getHashValue(cChrom)])
					pos = start + extension - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			}
			unsigned short int chromInt = getHashValue(cChrom);
//			bwRead sig(chromInt,pos,sval);
			bwRead sig(chromInt,pos);
			if(dataType == 0)
				signal_slist.push_back(sig);
			else if(dataType == 1)
				input_slist.push_back(sig);
		}else{
			formatFlag=1;
			num_skip++;
		}
	}

	if(i<6){
		cout << "WARNING:  Less than 6 BED columns detected.  ZINBA expects exactly 6 BED columns (see ZINBA wiki for import file format)" << endl;
	}else if(i>6){ 	
			cout << "WARNING:  More than 6 BED columns detected.  ZINBA expects exactly 6 BED columns (see ZINBA wiki for import file format)" << endl;
	}
	
	if(signal_slist.size() <=1 & dataType==0){
		cout << "ERROR: 0 experimental reads imported, check formatting of reads" << endl;
		cout << "ERROR: number of columns found is"<< i << endl;
		return(1);
	}
	if(input_slist.size()<=1 & dataType==1) {
		cout << "ERROR: 0 input reads imported, check formatting of reads" << endl;
		cout << "ERROR: number of columns found is"<< i << endl;
		return(1);
	}

	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;

}


int import::importCustomBed(const char * signalFile,int extension){
	
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
	//int extend = (int)(extension/2);
	int rval;
	int num_skip = 0;
	int score=0;
	while(!feof(fh)){
		rval = fscanf(fh,"%s%lu%lu%d",cChrom,&start,&stop,&score);
		if(rval == 4){
			/*if(strcmp(strand,minus) == 0){
				if(stop >= extension)
					pos = stop - extension + 1;
				else
					pos = 1;
			}else{*/  //Ignoring strand for now
				if((start + extension-1) <= chr_size[getHashValue(cChrom)])
					pos = start + extension - 1;
				else
					pos = chr_size[getHashValue(cChrom)];
			//}
			unsigned short int chromInt = getHashValue(cChrom);
			bwRead2 sig(chromInt,pos,score);
			custom_slist.push_back(sig);
		}else{
			num_skip++;
		}
	}
	fclose(fh);
	cout << "\tSkipped " << num_skip << " reads" << endl;
	return 0;
}


unsigned short int import::getHashValue(char *currChrom){
        map<const char*, int>::iterator i;
        i = chroms.find(currChrom);
        if(i == chroms.end()){
                char * chromosome = new char[128];
                strcpy(chromosome,currChrom);
                chroms[chromosome] = chromCounter;
                intsToChrom[chromCounter] = chromosome;
                return(chromCounter++);
        }else{
                return i->second;
        }
}

const char * import::getKey(unsigned short int chrom){
        map<int, const char*>::iterator i;
        i = intsToChrom.find(chrom);
        if(i == intsToChrom.end()){
                cout << chrom << endl;
                cout << "REALLY REALLY BAD ERROR!" << endl;
                exit(1);
        }else{
              	return i->second;
        }
}

