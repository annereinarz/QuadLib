#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "QuadRule.h"
#include "Quadrature.h"
#include "TensorQuad.h"

void OriginQuad(int* active, int size, QuadRule* qp){
   int i,j,k;
   for(i=1; i<size; i++){
      for(j=0; j<qp->n; j++){
          qp->t[active[i]][j]=qp->t[active[0]][j]*qp->t[active[i]][j];
      }
   }
   for(j=0; j<qp->n; j++){
      qp->wt[j]=qp->wt[j]*pow(qp->t[active[0]][j],size-1);
   }

   int lt=qp->n;
   double th;
   *qp=resize_quadrule_n(*qp, size);
   for(k=1; k<size; k++){
      for(j=0; j<lt; j++){
         for(i=0; i<qp->d; i++){
             qp->t[i][k*lt+j]=qp->t[i][j];
             qp->wt[k*lt+j]=qp->wt[j];
         }
         th=qp->t[active[0]][k*lt+j];
         qp->t[active[0]][k*lt+j]=qp->t[active[k]][k*lt+j];
         qp->t[active[k]][k*lt+j]=th;
      }
   }
}

void C2S(QuadRule* qp,int startidx, int endidx){
   //This routine transforms the quadrature rule from [0,1]^d to S_d
   int i,j;
   for(i=startidx+1; i<=endidx; i++){
       for(j=0; j<qp->n; j++){
           qp->t[i][j]*=qp->t[i-1][j];
       }
   }
   for(i=startidx+1; i<=endidx; i++){
      for(j=0; j<qp->n; j++){
           qp->wt[j]*=qp->t[i-1][j];
      }
   }
   for(j=0; j<qp->n; j++){           //loop over all quadrature points
      for(i=startidx; i<=endidx-1; i++){
         qp->t[i][j]-=qp->t[i+1][j];
      }
   }
}

void BasicSubQuad( int j, int k, int d, QuadRule* quadpoints){
// This routine realizes Steps 8,7,6 of the algorithm
//INPUT: j,k,d,t,wt
//OUTPUT: t,wt
int i,m;
int s1=0;
int s2=d-k;
int s3=2*d-2*k+j;

// Step 8
if(k==0){
   int active[2];
   active[0]=s1;  active[1]=s2;
   OriginQuad(active, 2, quadpoints);
}
else if(k==d){
      double* oldt=quadpoints->t[0];
      quadpoints->t[0]=quadpoints->t[j];
      quadpoints->t[j]=oldt;
}
else{
   int active[3];
   active[0]=s1;  active[1]=s2; active[2]=s3;
   OriginQuad(active,3, quadpoints);
}

// Step 7
if(k<d-1){      //vanishing a-b coordinates for k=d-1, k=d
   //a-s1 coordinates
   C2S(quadpoints ,1, d-k-1);
   for(i=0; i<quadpoints->n; i++){
      quadpoints->wt[i]*=pow(quadpoints->t[s1][i],d-k-1);
   }
   double* sigma=(double*)malloc((quadpoints->n)* sizeof(double));
   for(i=0; i<quadpoints->n; i++){
       sigma[i]=1;
       for(m=1;m<=d-k-1;m++){
           sigma[i]-=quadpoints->t[m][i];
       }
   }
   for(i=0; i<quadpoints->n; i++){
       for(m=1;m<=d-k-1;m++){
           quadpoints->t[m][i]*=quadpoints->t[s1][i];
       }
       quadpoints->t[s1][i]*=sigma[i];
   }
   //b-s2 coordinates
   C2S(quadpoints, d-k+1, 2*d-2*k-1);
   for(i=0; i<quadpoints->n; i++){
      quadpoints->wt[i]*=pow(quadpoints->t[s2][i],d-k-1);
   }
   for(i=0; i<quadpoints->n; i++){
       sigma[i]=1;
       for(m=d-k+1; m<=2*d-2*k-1; m++){
           sigma[i]-=quadpoints->t[m][i];
       }
   }
   for(i=0; i<quadpoints->n; i++){
       for(m=d-k+1; m<=2*d-2*k-1; m++){
           quadpoints->t[m][i]*=quadpoints->t[s2][i];
       }
       quadpoints->t[s2][i]*=sigma[i];
   }
}
if(k>0){    //vanishing p-q-s3 coordinates for k=0
   //p-q-s3 coordinates
   for(i=0; i<quadpoints->n; i++){
      quadpoints->wt[i]*=pow(quadpoints->t[s3][i], k-1);
   }
   if(j>0){
      C2S(quadpoints,2*d-2*k, 2*d-2*k+j-1);
      for(i=0; i<quadpoints->n; i++){
         for(m=2*(d-k);m<2*d-2*k+j;m++){
             quadpoints->t[m][i]*=-quadpoints->t[s3][i];
         }

      }
   }
   if(j<k-1){
     C2S(quadpoints,2*d-2*k+j+1, 2*d-k-1);
     double* sigma=(double*)malloc((quadpoints->n)* sizeof(double));
     for(i=0; i<quadpoints->n; i++){
         sigma[i]=1;
         for(m=2*d-2*k+j+1;m<2*d-k;m++){
             sigma[i]-=quadpoints->t[m][i];
         }
     }
     for(i=0; i<quadpoints->n; i++){
        for(m=2*d-2*k+j+1;m<2*d-k;m++){
           quadpoints->t[m][i]*=quadpoints->t[s3][i];
        }
        quadpoints->t[s3][i]*=sigma[i];
     }
   }
// Step 6
   C2S(quadpoints, 2*d-k, 2*d-1);
   double* sigma=(double*)malloc((quadpoints->n)* sizeof(double));
   for(i=0; i<quadpoints->n; i++){
       sigma[i]=1-quadpoints->t[s3][i];
       for(m=2*(d-k)+j+1;m<2*d-k;m++){
          sigma[i]-=quadpoints->t[m][i];
       }
   }
   for(i=0; i<quadpoints->n; i++){
      quadpoints->wt[i]*=pow(sigma[i],k);
   }
   for(i=0; i<quadpoints->n; i++){
      for(m=2*d-k;m<=2*d-1;m++){
           quadpoints->t[m][i]*=sigma[i];
      }
      if(j>0){
         for(m=2*d-k; m<=2*d-k+j-1; m++){
            quadpoints->t[m][i]-=quadpoints->t[m-k][i] ;
         }
      }
   } 
}

}

void BasicSubQuad_jacobi( int j, int k, int d, QuadRule* qp){
// This routine realizes Steps 7,6 of the algorithm
//INPUT: j,k,d,t,wt
//OUTPUT: t,wt
int i,m;
int s1=0;
int s2=d-k;
int s3=2*d-2*k+j;

// Step 7
if(k<d-1){      //vanishing a-b coordinates for k=d-1, k=d
   //a-s1 coordinates
   double* sigma=(double*)malloc((qp->n)* sizeof(double));
   for(i=0; i<qp->n; i++){
       sigma[i]=1;
       for(m=1;m<=d-k-1;m++){
           sigma[i]-=qp->t[m][i];
       }
   }
   for(i=0; i<qp->n; i++){
       for(m=1;m<=d-k-1;m++){
           qp->t[m][i]*=qp->t[s1][i];
       }
       qp->t[s1][i]*=sigma[i];
   }
   //b-s2 coordinates
   for(i=0; i<qp->n; i++){
       sigma[i]=1;
       for(m=d-k+1; m<=2*d-2*k-1; m++){
           sigma[i]-=qp->t[m][i];
       }
   }
   for(i=0; i<qp->n; i++){
       for(m=d-k+1; m<=2*d-2*k-1; m++){
           qp->t[m][i]*=qp->t[s2][i];
       }
       qp->t[s2][i]*=sigma[i];
   }
   free(sigma);
}
if(k>0){    //vanishing p-q-s3 coordinates for k=0
   //p-q-s3 coordinates
   if(j>0){
      for(i=0; i<qp->n; i++){
         for(m=2*(d-k);m<2*d-2*k+j;m++){
             qp->t[m][i]*=-qp->t[s3][i];
         }
      }
   }
   if(j<k-1){
     double* sigma=(double*)malloc((qp->n)* sizeof(double));
     for(i=0; i<qp->n; i++){
         sigma[i]=1;
         for(m=2*d-2*k+j+1;m<=2*d-k-1;m++){
             sigma[i]-=qp->t[m][i];
         }
     }
     for(i=0; i<qp->n; i++){
        for(m=2*d-2*k+j+1;m<=2*d-k-1;m++){
           qp->t[m][i]*=qp->t[s3][i];
        }
        qp->t[s3][i]*=sigma[i];
     }
     free(sigma);
   }
// Step 6
   double* sigma=(double*)malloc((qp->n)* sizeof(double));
   for(i=0; i<qp->n; i++){
       sigma[i]=1-qp->t[s3][i];
       for(m=2*(d-k)+j+1;m<=2*d-k-1;m++){
          sigma[i]-=qp->t[m][i];
       }
   }
   for(i=0; i<qp->n; i++){
      qp->wt[i]*=pow(sigma[i],k);
   }
   for(i=0; i<qp->n; i++){
      for(m=2*d-k;m<=2*d-1;m++){
           qp->t[m][i]*=sigma[i];
      }
      if(j>0){
         for(m=2*d-k; m<=2*d-k+j-1; m++){
            qp->t[m][i]-=qp->t[m-k][i];
         }
      }
   }
   free(sigma);
}
}




