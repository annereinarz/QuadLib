#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "QuadRule.h"
#include <float.h>

#define eps 0.3183099/*DBL_EPSILON */
#define pi 3.141592653589793

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
double* Function7(QuadRule q , int d, int k){
	int i, j;
	/* alpha = -2*d+k+epsilon;*/
	/* F = @(x,y,z,alpha) sqrt(sum(v.^2+(1+u).^2,2)).^alpha;*/
	double* F1;
	F1=(double*)malloc(q.n*sizeof(double));

	double* sum;
	sum=(double*)calloc(q.n, sizeof(double));

	double alpha=-2*d+k+eps;

	for(i=0; i<q.n; i++){
		for(j=0; j<d-1; j++){
			sum[i]+=pow(q.t[j][i]+q.t[d+j-1][i],2.0);
		}
      sum[i]+=pow(1.0+q.t[2*d-2][i],2.0);
	}

	for(i=0; i<q.n; i++){
		F1[i]=pow(sum[i], alpha/2.0);
	}
   free(sum);
	return F1;
}
double* Function11(QuadRule q , int d, int k){
	int i, j;
	/* alpha = -2*d+k+epsilon;*/
	/* F = @(x,y,z,alpha) sin(hat u)cos(hat z)sin(check u)sin(check v)sqrt(sum(z.^2,2)).^alpha;*/
	double* F1;
	F1=(double*)malloc(q.n*sizeof(double));

	double* sum;
	sum=(double*)calloc(q.n, sizeof(double));
	double* prod;
	prod=(double*)malloc(q.n*sizeof(double));
	for(i=0; i<q.n; i++){
		prod[i]=1;
	}

	double alpha=-2*d+k+eps;

	for(i=0; i<q.n; i++){
		for(j=0; j<d; j++){
			sum[i]+=q.t[2*d+j][i]*q.t[2*d+j][i];
		}
	}
   for(i=0; i<q.n; i++){
      for(j=0; j<k; j++){
         prod[i]*=sin(q.t[j][i])*cos(q.t[j+d][i]-q.t[j][i]);
      }
      for(j=max(k,0); j<d; j++){
         prod[i]*=sin(q.t[j][i])*sin(q.t[j+d][i]);
      }
   }
	for(i=0; i<q.n; i++){
		F1[i]=prod[i]*pow(sum[i], alpha/2.0);
	}
   free(sum); free(prod);
	return F1;
}

double* Function1(QuadRule q , int d, int k){
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
	        sum[i]+=q.t[j][i]*q.t[j][i];
	    }
   }
   for(i=0; i<q.n; i++){
       F1[i]=pow(sum[i], alpha/2.0);
   }
   free(sum);
   return F1;
}

double* Function5(QuadRule q , int d, int k){
   int i, j;
/* alpha = -d+epsilon;*/
/* F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^alpha;*/
   double* F5;
   F5=(double*)malloc(q.n*sizeof(double));

   double* sum;
   sum=(double*)calloc(q.n, sizeof(double));

   double alpha=-d+eps;

   for(i=0; i<q.n; i++){
	    for(j=0; j<d; j++){
	        sum[i]+=q.t[j][i]*q.t[j][i];
	    }
   }
   for(i=0; i<q.n; i++){
       F5[i]=pow(sum[i], alpha/2.0);
   }
   free(sum);
   return F5;
}

double* Function2(QuadRule q, int d, int k){
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
   free(sum);
   return F2;
}

double* Function3(QuadRule q, int d, int k){
/* F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^-1;*/
   double* F3;
   int i, j;
   F3=(double*)malloc(q.n*sizeof(double));
   double* sum;
   sum=(double*)calloc(q.n, sizeof(double));
   for(i=0; i<d; i++){
	   for(j=0; j<q.n; j++){
		   sum[j]+=q.t[i][j]*q.t[i][j];
	   }
   }
   for(i=0; i<q.n; i++){
	   F3[i]=pow(sum[i], -0.5);
   }
   free(sum);
   return F3;
}

double* Function4(QuadRule q, int d, int k){
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
   free(sum);
   return F4;
}







