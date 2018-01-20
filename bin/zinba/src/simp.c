#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define min(a,b)  ((a) < (b) ? (a) : (b))
#define maxx(x,y) ((x) > (y) ? (x) : (y))
#define MAX_LEN 1000000


typedef int elem_type ;
elem_type torben(elem_type *, int, int, int);

SEXP mkans(int *, int);
//int feval(int *, SEXP ,SEXP );

SEXP peakboundc(SEXP f, SEXP bpprofile, SEXP outputfile, SEXP rho){

  int i, j, lbasecount, pstart, pstop, med;
  int basecount[MAX_LEN];
  FILE *FI, *FO;
  char str_buf[MAX_LEN], line[MAX_LEN];
  char chr[127], ID[127], strand[2];
	char *delim = "\t";
	char *ptr= NULL;
  const char *input =CHAR(STRING_ELT(bpprofile,0));
  const char *output=CHAR(STRING_ELT(outputfile,0));


/////////////////
double sig;
double qVal;
/////////////////
	
  for(j=0; j<MAX_LEN;j++){
	basecount[j]=0;
  }

  FI = fopen(input, "r");	
  if(FI == NULL){		
    error("cannot open file %s\n", input);
  }

  FO = fopen(output, "w");	
  if(FO == NULL){
    error("cannot open file %s\n", output);
  }
  

//skip the first line of the coordinate file
    if(fgets(str_buf, MAX_LEN, FI) == NULL){
      error("there are no lines in file\n");
    }

//print headers to file
sprintf(line, "PEAKID\tChrom\tStart\tStop\tStrand\tSig\tMaxloc\tMax\tpStart\tpStop\tMedian\tqValue\n");
fputs(line, FO);


Rprintf("Begin Peak Refinement\n"); 
//read in a line, save the information in each, perform peakbounds, then print out, repeat for each line
int m=0;
 while(fgets(str_buf, MAX_LEN, FI) != NULL){

    i = 0;
     if ((ptr = strtok(str_buf, delim)) != NULL) {
    do {
      i++;
      if(i==1){
        strcpy(ID, ptr);
      }else if(i==2){
        strcpy(chr, ptr);
      }else if(i==3){
        pstart = atol(ptr);
      }else if(i==4){
        pstop = atol(ptr);
      }else if(i==5){
        strcpy(strand, ptr);
      }else if(i==6){
	sig = atof(ptr);
      }else if(i==7){
	qVal= atof(ptr);
      }else if(i>7){
		  basecount[i-8]  = atoi(ptr);
      }
      lbasecount=i-7; 
    } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", input);
  }
//Rprintf("%d\n",lbasecount);
 
  defineVar(install("x"), mkans(basecount, lbasecount), rho);
    int l=round(INTEGER(eval(f, rho))[0]);

  defineVar(install("x"), mkans(basecount, lbasecount), rho);
    med=round(INTEGER(eval(f, rho))[1]);
  int maxvec[l];
  for(i=1;i<=l;i++){
  	defineVar(install("x"), mkans(basecount, lbasecount), rho);
  	maxvec[i-1]=INTEGER(eval(f, rho))[i+1];

  }
////////////////////////////////////////////////////////////////////////

int lmaxvec = l;
int results[lmaxvec*2];
int max;
int searchlength=0;
int  peakstart=0;
int  peakend=0;
int  fit=0;
int  bound,h, n;
int maxvecbounds[lmaxvec*2];
double val,sumxy,sumx, sumx2, sumy, sumy2;

for(j=0; j<lmaxvec; j++){
	max=maxvec[j];
	double hold[lbasecount];
	for(i=0; i<lbasecount; i++){
		hold[i]=0;
	}
	searchlength=3000;		
		for(bound=50; bound<=min(searchlength, max); bound++){
			sumxy=0;
			sumx=0;
			sumy=0;
			sumx2=0;
			sumy2=0;
			n=bound;
			for(h=max-bound; h<=max-1;h++){ 
				sumxy=sumxy+(h+1)*basecount[h];
				sumx=sumx+h+1;
				sumy=sumy+basecount[h];
				sumx2=sumx2+pow(h+1,2);
				sumy2=sumy2+pow(basecount[h],2);						
			}
			hold[max-bound]=pow(sumxy-sumy*sumx/n,2)/((sumy2-sumy*sumy/n)*(sumx2-sumx*sumx/n));
		}
		//find max R2 position
		val=0;
		fit=0;
		for(i=0; i<lbasecount; i++){
        		if(hold[i] > val){
            		val = hold[i];
			fit=i;
        		}
			hold[i]=0;
    		}
		peakstart=fit+1;
	//right bound
	searchlength=3000;
		for(bound=50; bound<=min(searchlength, lbasecount-max); bound++){
			sumxy=0;
			sumx=0;			
			sumy=0;
			sumx2=0;
			sumy2=0;
			n=bound;
			for(h=max-1; h<max+bound-1;h++){ 
				sumxy=sumxy+(h+1)*basecount[h];
				sumx=sumx+h+1;
				sumy=sumy+basecount[h];
				sumx2=sumx2+pow(h+1,2);
				sumy2=sumy2+pow(basecount[h],2);	
					
			}
			hold[max-1+bound]=pow(sumxy-sumy*sumx/n,2)/((sumy2-sumy*sumy/n)*(sumx2-sumx*sumx/n));

			
			
		}
		//find max R2 position
		val=0;
		fit=0;
		for(i=0; i<lbasecount; i++){
        		if(hold[i] > val){
            		val = hold[i];
			fit=i;
        		}
		//	hold[i]=0;
    		}
		peakend=fit;
//save boundaries, adjust back to R index start at 1 
maxvecbounds[2*j]=peakstart;
maxvecbounds[2*j+1]=peakend;
}	

int dist, sep;
int hmaxvec[lmaxvec];

//final outputted vector of bounds after merging, contains zeros
//intialize results vector
for(i=0; i<2*lmaxvec;i++){
	results[i]=0;
}
results[0]=maxvecbounds[0];
results[1]=maxvecbounds[1];
hmaxvec[0]=maxvec[0];

h=0;
i=1;
	 
int dthresh = 100;
int sthresh = 100;
	 
while(i<lmaxvec){
	//calculate metrics per round
	dist=maxvec[i]-maxvec[i-1];
	sep=maxvecbounds[2*i]-results[2*h+1];
	
	//	if(dist<=250){
	if(dist<=dthresh){
	//merge if too close
		results[2*h]=min(results[2*h],maxvecbounds[2*i]);
		results[2*h+1]=maxx(results[2*h+1],maxvecbounds[2*i+1]);			if(basecount[hmaxvec[h]]<basecount[maxvec[i]]){		
			hmaxvec[h]=maxvec[i];
		}
		//	}else if(dist>250 & sep<100){
	}else if(dist>dthresh & sep<sthresh){
		//merge if not close enough but if bound overlap and no valley inbetween 
		//need to normalize with respect to min of whole basecount?
		results[2*h]=min(results[2*h],maxvecbounds[2*i]);
		results[2*h+1]=maxx(results[2*h+1],maxvecbounds[2*i+1])	;
		if(basecount[hmaxvec[h]]<basecount[maxvec[i]]){		
			hmaxvec[h]=maxvec[i];
		}
		//	}else if(dist>250 & sep>=100){
	}else if(dist>dthresh & sep>=sthresh){
		//separate peak here, no boundary overlap and too far 
		h=h+1;
		hmaxvec[h]=maxvec[i];
		results[2*h]=maxvecbounds[2*i];
		results[2*h+1]=maxvecbounds[2*i+1]	;
	}
	i=i+1;
}
/////////////////////////////////////////////////////////////////////////
int pos=0;

for(i=0;i<=h;i++){
	med=torben(basecount, lbasecount,results[2*i],results[2*i+1]);
	sprintf(line, "%s\t%s\t%d\t%d\t%s\t%.14f\t%d\t%d\t%d\t%d\t%d\t%.14f\n",ID, chr, pstart, pstop, strand,sig ,pstart+hmaxvec[i],basecount[hmaxvec[i]],pstart+results[2*i],pstart+results[2*i+1], med, qVal);
	fputs(line, FO);
}



}

 fclose(FO);
 fclose(FI);

  FILE *FI2,*FO2;
  FI2 = fopen(output, "r");	
  if(FI2 == NULL){
    error("cannot open file %s\n", output);
  }
  
  const char *outputbed=CHAR(STRING_ELT(outputfile,0));
  strcat(outputbed,".bed");

  FO2 = fopen(outputbed, "w");	
  if(FO2 == NULL){
    error("cannot open file %s\n", outputbed);
  }


if(fgets(str_buf, MAX_LEN, FI2) == NULL){
      error("there are no lines in file\n");
    }

 while(fgets(str_buf, MAX_LEN, FI2) != NULL){
    i = 0;
     if ((ptr = strtok(str_buf, delim)) != NULL) {
    do {
      i++;
      if(i==1){
        strcpy(ID, ptr);
      }else if(i==2){
        strcpy(chr, ptr);
      }else if(i==9){
        pstart = atol(ptr);
      }else if(i==10){
        pstop = atol(ptr);
      }else if(i==5){
        strcpy(strand, ptr);
      }else if(i==6){
	sig = atof(ptr);
      }else if (i==12){
	qVal=atof(ptr);
      }
    } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", input);
  }
  sprintf(line, "%s\t%d\t%d\t%s\t%.14f\t%s\n", chr, pstart, pstop, ID,sig ,strand);
  fputs(line, FO2);
}

  fclose(FO2);
  fclose(FI2);
  
 Rprintf("Peak Refinement Complete\n"); 
 return(bpprofile);
}



SEXP mkans(int *x, int len)
     {
         SEXP ans;
	 	
	 int l = len;
         PROTECT(ans = allocVector(INTSXP, l));
	 for(int i=0;i<l; i++){
         	INTEGER(ans)[i] = x[i];
	 }
         UNPROTECT(1);
         return ans;
     }
     
   /*  int feval(int *x, SEXP f,SEXP rho)
     {
         defineVar(install("x"), mkans(x), rho);
	 //PROTECT(ans = allocVector(INTSXP, l));
	 return(INTEGER(eval(f, rho))[1]);
     }
     
    SEXP zero(SEXP f, SEXP guesses, SEXP tol, SEXP rho)
     {
         int * x0 = INTEGER(guesses);
         int f0, f1, fc, xc;
     
        f0 = feval(x0, f, rho);
	defineVar(install("x"), mkans(x0), rho);
	int l=round(INTEGER(eval(f, rho))[0]);
	int maxvec[l];
	for(int i=1;i<=l;i++){
		defineVar(install("x"), mkans(x0), rho);
		maxvec[i-1]=INTEGER(eval(f, rho))[i];
	}

		
	return mkans(x0);
     }*/


/*
 * The following code is public domain.
 * Algorithm by Torben Mogensen, implementation by N. Devillard.
 * This code in public domain.
 */

elem_type torben(elem_type m[], int n, int exstart, int exstop)
{
    int         i, less, greater, equal,nnew=0;
    elem_type  min, max, guess, maxltguess, mingtguess;

        min = max = m[0] ;
    for (i=250 ; i<n-250 ; i++) {
         if(i<exstart | i>exstop){
        if (m[i]<min){ min=m[i];}
        if (m[i]>max){ max=m[i];}
           nnew++;
         }
    }
//        Rprintf("%d exstart %d exstop %d\n", nnew, exstart, exstop);
    while (1) {
        guess = (min+max)/2;
        less = 0; greater = 0; equal = 0;
        maxltguess = min ;
        mingtguess = max ;
        for (i=250; i<n-250; i++) {
                if(i<exstart | i>exstop){

            if (m[i]<guess) {
                less++;
                if (m[i]>maxltguess) maxltguess = m[i] ;
            } else if (m[i]>guess) {
                greater++;
                if (m[i]<mingtguess) mingtguess = m[i] ;
            } else equal++;
                }
}
        if (less <= (nnew+1)/2 && greater <= (nnew+1)/2) break ;
        else if (less>greater) max = maxltguess ;
        else min = mingtguess;
    }

    if (less >= (nnew+1)/2) return maxltguess;
    else if (less+equal >= (nnew+1)/2) return guess;
    else return mingtguess;
}


