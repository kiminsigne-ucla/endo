
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#define MAX_LEN 1025

extern "C" {
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;

extern "C" {
void mapout2alignability(char **Rfilelist,char **Routfile){
/*problems: 1) processes chrY twice at the end for some reason when running
		2) if you're not running the code in the directory the files are in
		the track name will put name=*directory/chr* instead of chr, leading
		to a seqfault in alignadjust.  Works when you cd to b.out directory*/
	const char * filelist = Rfilelist[0];
	const char * outfile = Routfile[0];

	FILE * fout;
	fout = fopen(outfile,"w"); 
	fprintf(fout,"track type=wiggle_0 name=\"%s\" desc=\"%s\" visibility=full\n",filelist,filelist);

	FILE * flist;
	flist = fopen(filelist,"r");
	char mapfile[128];
	char mapChar[1];
	int rval;
	while(!feof(flist)){
		rval = fscanf(flist,"%s",mapfile);
		Rprintf("\nProcessing %s .......\n",mapfile);
		string chrm = string(mapfile);
		size_t bout = chrm.find("b.out");
		chrm.erase(bout);
		Rprintf("\t%s .......\n",chrm.c_str());
		fprintf(fout,"fixedStep chrom=%s start=1 step=1\n",chrm.c_str());
		ifstream mfile (mapfile,ios::in|ios::binary);
		if(mfile.is_open()){
			while(!mfile.eof()){
				mfile.read(mapChar,1);
				unsigned char mchr = mapChar[0];
				fprintf(fout,"%d\n",(int) mchr);
			}
		}
		mfile.close();
	}
	fclose(flist);
}
}
