#include "cmath"
#include "glm_test.h"
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "math.h"
#include "R.h"





extern "C"{
int pglm_fit(int *Rfamily, int *RN, int* RM, const double *y, const double *prior, const double * offset, const double *X, const int *stratum, int *Rinit,int *rank, double *Xb, double *fitted, double *resid, double *weights, double *scale, int *df_resid, double *Rtheta){ 

int N=*RN;
int M=*RM;
int maxit=25;
double conv=0.00001;
int init=*Rinit; 
int failure=0;
int family=*Rfamily;
double theta=*Rtheta;
if(family==0){

failure=glm_nb(N, M, 1,y, prior, offset, X,stratum, maxit, conv,init, 
	    rank,Xb,fitted, resid, weights, 
	    scale, df_resid, Rtheta);
		if (failure == 1) {
//    			Rprintf("Glm.nb : Failure to converge\n");
		} 
} else if(family==2){
failure=glm_fit(POISSON, LOG, N, M, 1,y, prior, offset, X,stratum, maxit, conv,init, 
	    rank,Xb,fitted, resid, weights, 
	    scale, df_resid, theta);
		if (failure == 1) {
//    			Rprintf("Poisson : Failure to converge\n");
		} 
} else if(family==1){
failure=glm_fit(BINOMIAL, LOGIT, N, M, 1,y, prior, offset, X,stratum, maxit, conv,init, 
	    rank,Xb,fitted, resid, weights, 
	    scale, df_resid, theta);
		if (failure == 1) {
//    			Rprintf("Binomial : Failure to converge\n");
		} 
}



}
}
