#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "QuadRule.h"
#include "BasicSubQuad.h"
#include "TensorQuad.h"
#include "vtx.h"
#include "utils.h"
#include "AffineTrafo.h"
#include "functions.h"
#include "simpquad.h"
#include "Quadrature.h"
#include "quadrature.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))


int main(int argc, char **argv){
   /*Get Input data */
   assert_correct_number_of_inputs (argc, 3);
   // number of regular quadrature points1
   int nr = get_integer_argument(argv[1], "first input argument must be integer!");
   // select dimensions: first number gives dimension in which to integrate and the
   // second gives the dimension of the overlap
   int sim = get_integer_argument(argv[2], "third input argument must be integer!");
   //Give number of subquadrature to calculate
   int j=get_integer_argument(argv[3], "fourth imput argument must be integer!");
   int whichF=1;        // select which function to integrate, 1,2,3 or 4
   int singular = GJ;
   int regular  = GJ;

   double Q;

   Vertexlist vtxlist;
   vtxlist=get_vertexlist(sim);

   // compute space dimension d and dimension of the intersection k
   int d= min(vtxlist.s1, vtxlist.s2);
   int k= 2*d+1 - max(vtxlist.s1, vtxlist.s2);

   //Calculate Affine Transformation
   AffineTrafo A;
   A=determineAffineTrafo(d,k,vtxlist);

   //Calculate Quadrature Rule
   QuadRule QP;
   int ns;           // number of singular quadrature points
   if(singular==GJ){
      ns=nr;
   }
   else if(singular==CGL){
      ns=2*nr;
   }
   double alpha=findalpha(whichF, k, d);
   if(regular==GJ){
      if(k==-1){
         if(regular==GJ){
            QuadRule SP;
            SJacobi(&SP, nr, d);
            init_quadrule(&QP, SP.n*SP.n, SP.d+SP.d);
            TensorQuad_nonit(&QP, SP,SP);
            free_quadrule(SP);
         }
      }
      else{
         if(singular==GJ){
            QP = QuadratureRule_gaujacsing(k,d,ns,j,alpha);
         }
         else if(singular==CGL){
            QP = QuadratureRule_gaulegsing(ns);
         }
         QuadratureRule_gaujacreg2(k,d,j,nr, &QP);
      }

      int ndof=0;
      Q=squadtrafoj_jacobi(whichF, A, &QP,d,k,j, &ndof);
      printf("Q=%3.18lf;\nndof=%d;\n", Q, ndof);
      free_quadrule(QP);
   }
   else if(regular==GL){
      double alpha=findalpha(whichF, k, d);
      QuadRule QPreg, QPsing;
      if(k==-1){
         init_quadrule(&QPsing,nr,1);
         GLquad(&QPsing, nr, 0,1);
      }
      else if(singular==GJ){
         QPsing=QuadratureRule_gaujacsing(k,d,nr, j, alpha);
      }
      else if(singular==CGL){
         QPsing=QuadratureRule_gaulegsing(2*nr);
      }
      QPreg=gauleg_reg(d,nr);

      int i_por;
      int ndof=0;
      QuadRule QPh;
      init_quadrule(&QPh, 1, 1);

      for(i_por=0; i_por<QPsing.n; i_por++){
         init_quadrule(&QP, QPreg.n, QPreg.d+1);
         QPh.t[0][0]=QPsing.t[0][i_por]; QPh.wt[0]=QPsing.wt[i_por];
         TensorQuad_nonit(&QP, QPh, QPreg);
         Q=Q+squadtrafoj(whichF, A, &QP,d,k,j, &ndof);
         free_quadrule(QP);
      }
      printf("Q=%3.18lf;\n",Q);
      free_quadrule(QPreg); free_quadrule(QPsing);
   }

   return 0;
}
