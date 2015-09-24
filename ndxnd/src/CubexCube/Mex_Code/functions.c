#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "quadpoints.h"
#include <float.h>

#define eps 0.3183099/*DBL_EPSILON */
  
double* Function1(QuadraturePoints q , int d, int k){
   int i, j; 
/* alpha = -2*d+k+epsilon;*/
/* F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^alpha;*/
   double* F1;
   F1=(double*)malloc(q.n*sizeof(double));
   
   double* sum;
   sum=(double*)calloc(q.n, sizeof(double));

   double alpha=-2*d+k+eps;

   for(i=0; i<q.n; i++){
	    for(j=0; j<d; j++){
	        sum[i]+=q.t[2*d+j][i]*q.t[2*d+j][i];
	    }
   }
   for(i=0; i<q.n; i++){
       F1[i]=pow(sum[i], alpha/2.0);
   }
   return F1;
}

double* Function2(QuadraturePoints q, int d, int k){
    int i, j;
/* alpha = -2*d+k+epsilon;*/
/* F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^alpha;*/
   double* F2;
   F2=(double*)malloc(q.n*sizeof(double));
   double alpha=-2*d+k+eps;

   double* sum;
   sum=(double*)calloc(q.n, sizeof(double));

   for(i=0; i<d; i++){
	   for(j=0; j<q.n; j++){
		   sum[j]+=(q.t[i][j]-q.t[d+i][j])*(q.t[i][j]-q.t[d+i][j]);
	   }
   }
   for(i=0; i<q.n; i++){
	 F2[i]=pow(sum[i], alpha/2.0);
   }
   return F2;
}

double* Function3(QuadraturePoints q, int d, int k){
/* F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^-1;*/
   double* F3;
   int i, j; 
   F3=(double*)malloc(q.n*sizeof(double));
   double* sum;
   sum=(double*)calloc(q.n, sizeof(double));      
   for(i=0; i<d; i++){
	   for(j=0; j<q.n; j++){
		   sum[j]+=q.t[2*d+i][j]*q.t[2*d+i][j];
	   }
   }
   for(i=0; i<q.n; i++){
	   F3[i]=pow(sum[i], -0.5);
   }
   return F3;
}

double* Function4(QuadraturePoints q, int d, int k){
/* F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^-1;*/
   double* F4;
   F4=(double*)malloc(q.n*sizeof(double));
   int i, j;
   double* sum;
   sum=(double*)calloc(q.n, sizeof(double));  
   for(i=0; i<d; i++){
	   for(j=0; j<q.n; j++){
		   sum[j]+=(q.t[i][j]-q.t[d+i][j])*(q.t[i][j]-q.t[d+i][j]);
	   }
   }
   for(i=0; i<q.n; i++){
       F4[i]=pow(sqrt(sum[i]),-1.0);
   }
   return F4;
}







