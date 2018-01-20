#include "pc.h"
#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <functional>
#include <utility>
#include <ext/hash_map>

extern "C" {
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include "Rdefines.h"
}

using namespace std;
using namespace __gnu_cxx; 

class lessAbsoluteValue {
public:
  bool operator()(int a, int b) const {
    return abs(a) < abs(b);
  }
};

extern "C" {
// read in regular eland files, adjusting the negative strand coordinate by sequence length
SEXP read_aligned_tags(SEXP filename,SEXP filetype){
	const char* fname=CHAR(asChar(filename));
	const char* ftype=CHAR(asChar(filetype));
	// main data vector
	// chr - pos
	vector< vector<int> > pos;
	vector< vector<int> > posnm; // number of mismatches

	// chromosome map
	hash_map<string, int, hash<string>,equal_to<string> > cind_map;
	vector<string> cnames;
	const char * bowtie = "bowtie";
	const char * bed = "bed";
	const char * tagAlign = "tagAlign";

	FILE *f=fopen(fname,"r");
	if (!f){
		Rprintf("Can't open input file %s\n",fname);
	}else if(strcmp(ftype,bowtie) != 0 && strcmp(ftype,bed) != 0 && strcmp(ftype,tagAlign) != 0){
		Rprintf("Unrecognized format: %s [bowtie|bed]",ftype);
	}else{
		Rprintf("Opened %s\n",fname);
		// read in bowtie line
		int fcount=0;
		char cChrom[128];
		long unsigned int fpos;
		int nm;
		char strand[1];
		char minus[] = "-";
		char line[512];char seq[128];
		string chr;string mm;
		unsigned long int start;unsigned long int stop;
		char name[128];
		int rval;char qscore[128];int eval;
		while(!feof(f)){
			if (strcmp(ftype,bowtie) == 0){
				rval = fscanf(f,"%s%s%s%lu%s%s%d",name,strand,cChrom,&fpos,seq,qscore,&eval);
				fgets(line,512,f);
				if(rval == 7){
					if(strcmp(strand,minus) == 0)
						fpos = -1 * (fpos + strlen(seq) - 1);
					// determine number of mismatches
					mm = string(line);
					nm=0;
					if(mm.size()>0) {
						nm++;
						string::size_type tp(0);
						while(tp!=string::npos) {
							tp = mm.find(",",tp);
							if(tp!=string::npos) {
								tp++;
								++nm;
							}
						}
					}
				}
			}else if(strcmp(ftype,bed) == 0){
				rval = fscanf(f,"%s%lu%lu%s%d%s",cChrom,&start,&stop,name,&eval,strand);
				if(rval == 6){
					if(strcmp(strand,minus) == 0){
						fpos = -1 * stop;
					}else{
						fpos = start;
					}
					nm = 1;
				}
			}else if(strcmp(ftype,tagAlign) == 0){
				rval = fscanf(f,"%s%lu%lu%s%d%s",cChrom,&start,&stop,name,&eval,strand);
				if(rval == 6){
					if(strcmp(strand,minus) == 0){
						fpos = -1 * stop;
					}else{
						fpos = start;
					}
					nm = 1;
				}
			}
			// determine the chromosome index
			chr = string(cChrom);
			hash_map<string, int, hash<string>,equal_to<string> >::const_iterator li=cind_map.find(chr);
			int cind=-1;
			if(li==cind_map.end()) {
				// register new chromosome
				cind=cnames.size();
				cnames.push_back(chr);
				cind_map[chr]=cind;
				// allocate new pos vector
				pos.push_back(vector<int>());
				posnm.push_back(vector<int>());
			}else{
				cind=li->second;
			}
			fcount++;
			(pos[cind]).push_back(fpos);
			(posnm[cind]).push_back(nm);
			
		}
		fclose(f);
		Rprintf("done. read %d fragments\n",fcount);
	}
	// construct output structures
	SEXP chnames;
	int np=0; // number of protections
	PROTECT(chnames = allocVector(STRSXP, cnames.size()));
	for(vector<string>::const_iterator csi=cnames.begin();csi!=cnames.end();++csi) {
		SET_STRING_ELT(chnames, csi-cnames.begin(), mkChar(csi->c_str()));
	}
	np++;

	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, cnames.size()));   np++;
	vector<vector<int> >::const_iterator nsi;
	vector<vector<string> >::const_iterator ssi;
	for(vector<vector<int> >::const_iterator csi=pos.begin();csi!=pos.end();++csi) {
		nsi=posnm.begin()+(csi-pos.begin());
		SEXP dv,dnames_R;
		PROTECT(dnames_R = allocVector(STRSXP, 2)); np++;
		SET_STRING_ELT(dnames_R, 0, mkChar("t"));
		SET_STRING_ELT(dnames_R, 1, mkChar("n"));
    
		SEXP tv,nv,sv;
		PROTECT(tv=allocVector(INTSXP,csi->size()));   np++;
		PROTECT(nv=allocVector(INTSXP,csi->size()));   np++;
		int* i_tv=INTEGER(tv);
		int* i_nv=INTEGER(nv);
    
		int i=0;
		vector<int>::const_iterator ini=nsi->begin();
		for(vector<int> ::const_iterator pi=csi->begin();pi!=csi->end();++pi) {
			i_tv[i]=*pi;
			i_nv[i]=*ini++;
			i++;
		}
		PROTECT(dv = allocVector(VECSXP, 2));   np++;
		SET_VECTOR_ELT(dv, 0, tv);
		SET_VECTOR_ELT(dv, 1, nv);
		setAttrib(dv, R_NamesSymbol, dnames_R);
		SET_VECTOR_ELT(ans, csi-pos.begin(), dv);
	}
	setAttrib(ans,R_NamesSymbol,chnames);
	UNPROTECT(np);
	return(ans);
}
}
