#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "quadpoints.h"
#include "Tensorquad.h"
#include "AffineTrafo.h"
#include "BasicSubQuad.h"
#include "functions.h"
#include "squad.h"

#include "clapack.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))

QuadraturePoints QuadratureRule(int k, int d, int nr, int ns){
   //Compute quadrature rule on [0,1]^(2*d)
   int i,lx;
   int powh;
   QuadraturePoints QP, qp;
   init_quadpoints(&qp ,nr ,1);
   GLquad(&qp,nr, 0, 1);
   if(k==-1){
      powh =(int) pow(nr, 2*d),
      init_quadpoints(&QP, powh, 2*d);
      lx=nr;
      GLquad(&QP,nr, 0, 1);     //only regular
   }
   else{
      lx=0.5*ns*(ns+1);
      powh=(int)pow(nr,2*d-1) ;
      init_quadpoints(&QP, lx*powh, 2*d);
      CGLquad(&QP, ns);         //only singular
   }
   for(i=1; i<2*d; i++){
      TensorQuad(i,lx, &QP, qp);
   }

   free_qp(qp);
   return QP;
}

double squad_nonperm(int whichF, AffineTrafo A, int d, int k, int nr, int ns){
      //This routine realises a portionwise computation and evaluation of the quadrature rule
      double Q=0;

     /* printf("k=%d, d=%d, nr=%d, ns=%d, whichF=%d\n",k,d,nr,ns,whichF);
      {
         int i,j;
         for(i=0; i<d; i++){
            for(j=0; j<d; j++){
               printf("A1[%d][%d]=%lf ",i,j,A.A1[i][j]);
            }
            printf("\n");
         }
         for(i=0; i<d; i++){
            for(j=0; j<d; j++){
               printf("A2[%d][%d]=%lf ",i,j,A.A2[i][j]);
            }
            printf("\n");
         }
         for(i=0; i<d; i++){
               printf("v1[%d]=%lf \n",i,A.v10[i]);
         }
         for(i=0; i<d; i++){
               printf("v2[%d]=%lf \n",i,A.v20[i]);
         }
      } */

      //Calculate Quadrature Rule on [0,1]^2*d
      QuadraturePoints QP;
      QP=QuadratureRule(k,d,nr, ns);

      if(k==-1){
         C2S(&QP, 0, d-1);
         C2S(&QP, d, 2*d-1);
      }
      else{
         BasicSubQuad(0,0,d,&QP);
      }
      int* P=0, flag=0;
      PermRefl( k, d, &QP, A, P, &Q, flag, whichF);

      return Q;

}

double squad_perm(int whichF, AffineTrafo A, int d, int k, int nr, int ns, int j, int perm, int process){
      int i;
      double Q=0;
      int flag;
      //Quadrature Rule on [0,1]^2*d
      QuadraturePoints QP;
      //Permutation
      int*  P;

      P=(int*)malloc(2*d*sizeof(int));
      for(i=0; i<2*d; i++){
         P[i]=i+1;
      }
      for(i=0; i<perm; i++){
         int check=NextPerm(j, k,d, P);
         if(check!=0) {
            printf("Wrong permutation number");
         }
      }
      QP=QuadratureRule(k,d,nr, ns);
      BasicSubQuad(j,k,d,&QP);
      //  Permutation
      flag=1;
      PermRefl( k, d, &QP, A, P, &Q, flag, whichF);
      flag=2;
      PermRefl( k, d, &QP, A, P, &Q, flag, whichF);
      free_qp(QP);

      return Q;
}
