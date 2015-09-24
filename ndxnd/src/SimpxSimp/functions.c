#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "QuadRule.h"
#include "functions.h"

#define eps 0.318309886183791//DBL_EPSILON

double findalpha(int w, int k, int d){
 if(w==1||w==2){
    return -2*d+k+eps;
 }
 else if(w==5){
     return -d+eps;
 }
 else if(w==3){
    return -0.5;
 }
 else if(w==4){
    return -1;
 }
 printf("invalid whichF\n");
 return 0;
}


double* Function1(Quad2 q2 , int d, int k){
   int i, j;
// alpha = -2*d+k+epsilon;
// F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^alpha;
   double* F1;
   F1=(double*)malloc(q2.rzh*sizeof(double));

   double* sum;
   sum=(double*)calloc(q2.rzh, sizeof(double));

   double alpha=-2*d+k+eps;

   for(i=0; i<q2.rzh; i++){
	    for(j=0; j<q2.lzh; j++){
	        sum[i]+=q2.zh[j][i]*q2.zh[j][i];
	    }
   }
   for(i=0; i<q2.rzh; i++){
       F1[i]=pow(sum[i], alpha/2.0);
   }
   free(sum);
   return F1;
}

double* Function2(Quad2 q2, int d, int k){
    int i, j;
// alpha = -2*d+k+epsilon;
// F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^alpha;
   double* F2;
   F2=(double*)malloc(q2.ruv*sizeof(double));
   double alpha=-2*d+k+eps;

   double* sum;
   sum=(double*)calloc(q2.ruv, sizeof(double));

   for(i=0; i<q2.luv; i++){
	   for(j=0; j<q2.ruv; j++){
		   sum[j]+=(q2.u[i][j]-q2.v[i][j])*(q2.u[i][j]-q2.v[i][j]);
	   }
   }
   for(i=0; i<q2.ruv; i++){
	 F2[i]=pow(sum[i], alpha/2.0);
   }
   free(sum);
   return F2;
}

double* Function3(Quad2 q2, int d, int k){
// F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^-1;
   double* F3;
   int i, j;
   F3=(double*)malloc(q2.rzh*sizeof(double));
   double* sum;
   sum=(double*)calloc(q2.rzh, sizeof(double));
   for(i=0; i<q2.lzh; i++){
	   for(j=0; j<q2.rzh; j++){
		   sum[j]+=q2.zh[i][j]*q2.zh[i][j];
	   }
   }
   for(i=0; i<q2.rzh; i++){
	   F3[i]=pow(sum[i], -0.5);
   }
   free(sum);
   return F3;
}

double* Function4(Quad2 q2, int d, int k){
// F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^-1;
   double* F4;
   F4=(double*)malloc(q2.ruv*sizeof(double));
   int i, j;
   double* sum;
   sum=(double*)calloc(q2.ruv, sizeof(double));  
   for(i=0; i<q2.luv; i++){
	   for(j=0; j<q2.ruv; j++){
		   sum[j]+=(q2.u[i][j]-q2.v[i][j])*(q2.u[i][j]-q2.v[i][j]);
	   }
   }
   for(i=0; i<q2.ruv; i++){
       F4[i]=pow(sqrt(sum[i]),-1.0);
   }
   free(sum);
   return F4;
}

double* Function5(Quad2 q2 , int d, int k){
   int i, j;
// alpha = -d+epsilon;
// F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^alpha;
   double* F5;
   F5=(double*)malloc(q2.rzh*sizeof(double));

   double* sum;
   sum=(double*)calloc(q2.rzh, sizeof(double));

   double alpha=-d+eps;

   for(i=0; i<q2.rzh; i++){
	    for(j=0; j<q2.lzh; j++){
	        sum[i]+=q2.zh[j][i]*q2.zh[j][i];
	    }
   }
   for(i=0; i<q2.rzh; i++){
       F5[i]=pow(sum[i], alpha/2.0);
   }
   free(sum);
   return F5;
}
