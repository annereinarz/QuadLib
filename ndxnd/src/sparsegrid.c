#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include "math.h"
#include <string.h>

#include "QuadRule.h"
#include "NextPoint.h"
#include "unique.h"
#include "TensorQuad.h"
#include "sparsegrid.h"
#include "quadrature.h"

int maximum(int *K, int d){
 int i;
 int max;
 max=K[0];
 for(i=1; i<d; i++){
    if(K[i]>max){
       max=K[i];
    }
 }
 return max;
}

bool elemIl(int d, int* K, int l, int p, int T){
   if(p==-1){
      int i;
      int sum=0;
      for(i=0; i<d; i++){
         sum+=pow(2,K[i]+1);
      }
      if(sum>pow(2,l)+2*(d-1)){
         return false;
      }
      return true;
   }
   else if(p==1){
     int i;
     int sum=0;
     for(i=0; i<d; i++){
          sum+=K[i]+1;
     }
     if(sum>l+d-1){
        return false;
     }
     return true;
   }
   else if(p==0){   //full tensor
      int i;
      for(i=0; i<d; i++){
         if(K[i]+1>l){
            return false;
         }
      }
      return true;

   }
   else if(p==-4){
     int i;
     int sum=0;
     for(i=0; i<d;i++){
        sum+=K[i]+1;
     }
     int max;
     max=maximum(K,d);
     if(-sum+T*max<-(d+l-1)+T*l){
        return false;
     }
     return true;
   }
   else{
      printf("invalid choice of p\n");
      return false;
   }
}

void set_up_1dquadrules(int l, QuadRule *QP1D, int regular){
   int i;
   for(i=0; i<l; i++){
      if(regular==CC){ /*Clenshaw-Curtis*/
         char fname[17];
         sprintf(fname, "1DQuadRules/CC%d",i+1);
         import_1d_quadrule(fname,&QP1D[i]);
      }
      else if(regular==KP){  /*Kronrod-Patterson*/
         init_quadrule(&QP1D[i], (int)pow(2,i+1)-1,1);
         KPquad(i+1,&QP1D[i]);
         int j;
         for(j=0; j<QP1D[i].n; j++){
            QP1D[i].wt[j]*=0.5;
            QP1D[i].t[0][j]=(QP1D[i].t[0][j]+1)*0.5;
         }
      }
      else if(regular==GL){   /*Gauss-Legendre*/
         init_quadrule(&QP1D[i], (int)pow(2,i+1)-1,1);
         GLquad(&QP1D[i],(int)pow(2,i+1)-1,0,1);
      }
      else{
         fprintf(stderr,"Error: Invalid choice in --arguments.regular.\n");
         exit(1);
      }
   }
}

int next_sparse(QuadRule *QPtemp, int l, QuadRule* QP1D, int* K, int d, int p, int T){
   int i;
   if(! morePoints(K)){
      return 0;
   }
   /*Remove invalid indices*/
   while(! elemIl(2*d-1,K, l,p, T)){
      if (! NextPoint(K,2*d-1,l)) return 0;
   }
   /*Calculate coefficients and sort out zeros*/
   int coeff=0;  int absz=0;
   int *Kplusz;
   Kplusz=(int*)malloc((2*d-1)*sizeof(int));
   if(Kplusz==NULL){
      fprintf(stderr,"Out of memory\n");
      exit(1);
   }
   memcpy(Kplusz, K, (2*d-1)*sizeof(int));
   char *z;
   z=(char*)calloc(2*d-1,sizeof(char));
   if(z==NULL){
      fprintf(stderr,"Out of memory\n");
      exit(1);
   }
   for(i=0; i<(int)pow(2,2*d-1); i++){
      if(elemIl(2*d-1,Kplusz,l,p,T)){
         if(absz % 2){
            coeff--;
         }
         else{
            coeff++;
         }
      }

      /*Update z, absz, and Kplusz*/
      int j=0;
      while(j<2*d-1){
         z[j] = (z[j] + 1) % 2;
         if(z[j]!=0){
            absz++;
            Kplusz[j]++;
            break;
         }
         absz--;
         Kplusz[j]--;
         j++;
      }
   }
   if(coeff==0){
      if (! NextPoint(K,2*d-1,l)) return 0;
      return next_sparse(QPtemp,l,QP1D,K,d,p,T);
   }
   /*Tensor Product quadrature rules*/
   init_quadrule(QPtemp,QP1D[K[0]].n,1);
   for(i=0; i<QPtemp->n; i++){
      QPtemp->t[0][i]=QP1D[K[0]].t[0][i];
      QPtemp->wt[i]  =QP1D[K[0]].wt[i];
   }
   for(i=1; i<2*d-1; i++){
      TensorQuad_sparse(QPtemp, QP1D[K[i]]);
   }

   for(i=0; i<QPtemp->n; i++){
      QPtemp->wt[i]*=coeff;
   }
   NextPoint(K,2*d-1,l);
   return 1;
}

void append_quadrule(QuadRule *QPreg, QuadRule QPtemp){
      int oldn=QPreg->n;
      int i,j;
      QPreg->n=QPtemp.n+QPreg->n;
      if(QPreg->d!=QPtemp.d){
         fprintf(stderr, "Error: Quadrature rules must be of the same dimension.\n");
         exit(1);
      }
      for(i=0; i<QPreg->d; i++){
         QPreg->t[i]=(double*)realloc(QPreg->t[i], QPreg->n*sizeof(double));
      }
      QPreg->wt=realloc(QPreg->wt, QPreg->n*sizeof(double));
      for(i=0; i<QPtemp.n; i++){
         for(j=0; j<QPreg->d; j++){
            QPreg->t[j][oldn+i]=QPtemp.t[j][i];
         }
         QPreg->wt[oldn+i]=QPtemp.wt[i];
      }
}

QuadRule Sparse(int d, int l, int p, int regular, int T){
   int i;
   QuadRule QPreg, QP1D[l], QPtemp;
   init_quadrule(&QPreg,0,2*d-1);

   set_up_1dquadrules(l,QP1D,regular);

   int *K;
   K=(int*)malloc((2*d-1)*sizeof(int));
   memset(K,0,(2*d-1)*sizeof(int));

   while(1) {
      if (! next_sparse(&QPtemp,l,QP1D, K, d,p, T)){
         break;
      }
      append_quadrule(&QPreg, QPtemp);
      free_quadrule(QPtemp);
   }

   //printf("n=%d Qpreg.d=%d\n d=%d\n",QPreg.n,QPreg.d,d);
   unique(&QPreg);
   //printf("n=%d Qpreg.d=%d\n d=%d\n",QPreg.n,QPreg.d,d);
   for(i=0; i<l; i++){
      free_quadrule(QP1D[i]);
   }
   return QPreg;
}
