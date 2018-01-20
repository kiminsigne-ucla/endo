#include "cmath"
#include <cstdlib>
#include <R.h>
extern "C" {
#include <math.h>
#include "glm_test.h"
#include "asa103.h"
#include "asa121.h"
#define max(a,b) (((a) > (b)) ? (a) : (b))

/* Fit negbin GLM 

Input:

family       GLM family (see below)
link         Link function (see below)
N            # units
M            # X variables
S            # strata (0 means no intercept)
y            y-variable (N-vector)
prior        prior weights (if present)
X            If M>0, N*M matrix of X variables
stratum      If S>1, stratum assignments coded 1...S (N-vector)
maxit        Maximum number of iterations of IRLS algorithm
conv         Proportional change in weighted sum of squares residuals to
             declare convergence
init         If true (non-zero), the iteration starts from initial estimates 
             of fitted values (see below). This option has no effect if
	     no iteration is required

Output:

rank         rank of X after regression on strata
Xb           orthogonal basis for X space (N*rank matrix)
fitted       fitted values 
resid        working residuals (on linear predictor scale) (N-vector)
weights      weights (N-vector)
scale        scale factor (scalar)
df_resid     residual degrees of freedom

Return

0            convergence
1            no convergence after maxit iterations

*/

int glm_nb(int N, int M, int S,
	    const double *y, const double *prior, const double * offset, const double *X, 
	    const int *stratum, int maxit, double conv, int init, 
	    int *rank, double *Xb, 
	    double *fitted, double *resid, double *weights, 
	    double *scale, int *df_resid, double *Rtheta){ 


 
int failure=0;
double theta=*Rtheta;
/*calculation of convenience variables*/
double n=0;//sum prior
int i=0;
for(i=0; i<N; i++){n=n+prior[i];}
//If no starting theta it is passed as less than zero.
  if(theta<0){
	//Initial fit of data to poisson 
	failure=glm_fit(POISSON, LOG, N, M, 1, y, prior, offset, X, stratum, maxit, conv, 0,  //input
			rank, Xb, fitted, resid, weights, scale,df_resid, theta);   //output

  }else{
	/*else use existing theta to fit NB*/
	failure=glm_fit(NEGBIN, LOG, N, M, 1, y, prior, offset, X, stratum, maxit, conv, 1,  //input
			rank, Xb, fitted, resid, weights, scale,df_resid, theta);   //output
  }
  /*initial estimate of theta based on poisson fitted values*/ 

  theta=theta_ml( N, n, y, fitted, prior, maxit);
  int iter=0;
  double d1=sqrt(2.0*max(1.0,*df_resid));
  double del=1.0;
  double d2=del;
  double Lm=loglik(N, y, fitted, prior, theta);
  double Lm0=Lm + 2.0*d1;
  double t0=0.0;
  while( (iter<maxit) && (fabs(Lm0-Lm)/d1 +abs(del)/d2) > pow(10.0, -8)){
	//update fitted values in *fitted
	failure=glm_fit(NEGBIN, LOG, N, M, 1, y, prior, offset, X, stratum, maxit, conv, 1,  //input
			rank, Xb, fitted, resid, weights, scale,df_resid, theta);   //output
	//set current value of theta to t0
	t0=theta;
	//update theta
	theta=theta_ml( N, n, y, fitted, prior, maxit);
	del=t0-theta;
	Lm0=Lm;
	Lm=loglik(N, y, fitted, prior, theta);
	iter++;
  }	
	*Rtheta=theta;
  if(iter>maxit){Rprintf("Glm_nb alternation limit reached\n");}
return(failure);
}


/*Log liklihood of NB*/
  double loglik(int N,const double  *y, double *fitted, const double *prior, double theta){
  double ll=0;
  int i=0;
  for(i=0; i<N;i++){
  	ll=ll+prior[i]*((y[i] + theta) * log(fitted[i] + theta) - y[i] * log(fitted[i]) + lgamma(y[i] + 1) - theta * log(theta) + lgamma(theta) - lgamma(theta+y[i]));
  }
  return(ll);
}


/* Maximum Likelihood Estimation of Theta*/
//notes
//check to see if lgamma approximation to digamma is much faster
//look into NBR book's alternative estimation of theta, may be much faster

double theta_ml(int N, double n, const double *y, double *mu, const double *prior, int maxit)
{
  double tol=0.0001220703;
  double t0=0;
  double del=1;
  int it=0;
  int flag;
  double score=0;
  double info=0;
  int i=0;
  //compute initial for theta (t0)
  for(i=0; i<N;i++){
	t0=t0+prior[i]*(y[i]/mu[i]-1)*(y[i]/mu[i]-1);
  }
  t0=n/t0;
  while((it < maxit) && (fabs(del) > tol)){
  	t0= fabs(t0);
	for(i=0; i<N;i++){score=score+ prior[i]*(digama(t0 + y[i], &flag) - digama(t0, &flag) + log(t0) +1.0 - log(t0 + mu[i]) - (y[i] + t0)/(mu[i] + t0));
	info=info + prior[i]*( - trigam(t0 + y[i], &flag) + trigam(t0, &flag) - 1/t0 + 2/(mu[i] + t0) - (y[i] + t0)/((mu[i] + t0)*(mu[i] + t0)));}
	del=score/info;
        t0=t0 + del;
	it++;
	score=0;
	info=0;
  } 
   if(t0 < 0){t0= 0;}
//   if(it==maxit){ Rprintf("iteration Limit Reached\n");}
   return(t0);
}

}

