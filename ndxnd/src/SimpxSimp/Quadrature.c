#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "QuadRule.h"
#include "TensorQuad.h"
#include "quadrature.h"

QuadRule QuadratureRule_gaulegsing(int ns){
   QuadRule QP;
   init_quadrule(&QP, (int)(0.5*ns*(ns+1)), 1);
   CGLquad(&QP, ns);
   return QP;
}

QuadRule QuadratureRule_cgl_gl(int k, int d, int nr, int ns){
   //Compute quadrature rule on [0,1]^(2*d)
   int i,lx;
   int powh;
   QuadRule QP, qp;
   init_quadrule(&qp ,nr ,1);
   GLquad(&qp,nr, 0, 1);
   if(k==-1){
      powh =(int) pow(nr, 2*d),
      init_quadrule(&QP, powh, 2*d);
      lx=nr;
      GLquad(&QP,nr, 0, 1);     //only regular
   }
   else{
      lx=0.5*ns*(ns+1);
      powh=(int)pow(nr,2*d-1) ;
      init_quadrule(&QP, lx*powh, 2*d);
      CGLquad(&QP, ns);         //only singular
   }
   for(i=1; i<2*d; i++){
      TensorQuad(i,lx, &QP, qp);
   }

   free_quadrule(qp);
   return QP;
}

QuadRule QuadratureRule(int k, int d,int nr, int ns){
   //Compute quadrature rule on [0,1]^(2*d)
   int i,lx;
   int powh;
   QuadRule QP, qp;
   init_quadrule(&qp ,nr ,1);
   GLquad(&qp,nr, 0, 1);
   if(k==-1){
      powh =(int) pow(nr, 2*d),
      init_quadrule(&QP, powh, 2*d);
      lx=nr;
      GLquad(&QP,nr, 0, 1);     //only regular
   }
   else{
      lx=0.5*ns*(ns+1);
      powh=(int)pow(nr,2*d-1) ;
      init_quadrule(&QP, lx*powh, 2*d);
      CGLquad(&QP, ns);         //only singular
   }
   for(i=1; i<2*d; i++){
      TensorQuad(i,lx, &QP, qp);
   }

   free_quadrule(qp);
   return QP;
}

void SJacobi(QuadRule* QP, int n, int d){
   int i,l,l2;
   QuadRule qp;
   init_quadrule(QP, pow(n,d), d);
   GJquad01(QP, n, 0, d-1);
   for(i=2; i<=d; i++){
      init_quadrule(&qp, n, 1);
      GJquad01(&qp, n, 0, d-i);
      TensorQuad(i-1, n, QP, qp);
      for(l=0; l<QP->n; l++){     //chi_1
          QP->t[i-1][l]=QP->t[i-1][l]*QP->t[i-2][l];
      }
      free_quadrule(qp);
   }
   for(l=0; l<QP->n; l++){      //chi_2
       for(l2=0; l2<QP->d-1; l2++){
           QP->t[l2][l]=QP->t[l2][l]-QP->t[l2+1][l];
       }
   }
}

void QuadratureRule_gaujacreg2(int k, int d, int j, int nr, QuadRule *qp){

int s2=d-k;
int s3=2*d-2*k+j;

 // Step 8
for(int i=0; i<qp->n; i++){
   qp->wt[i]=qp->wt[i]*pow(qp->t[0][i],2*d-k-1);
}
if(k==0){
  QuadRule SP1;
  init_quadrule(&SP1,nr,1);
  GJquad01(&SP1,nr,0,d-1);
  TensorQuad_nonit_ip(qp,SP1);
  *qp=resize_quadrule_n(*qp, 2);
  for(int i=0; i<(int)(qp->n/2); i++){
     qp->wt[i+(int)(qp->n/2)]=qp->wt[i];
     qp->t[1][i]=qp->t[1][i]*qp->t[0][i];
     qp->t[0][i+(int)(qp->n/2)]=qp->t[1][i];
     qp->t[1][i+(int)(qp->n/2)]=qp->t[0][i];
  }
  free_quadrule(SP1);
}
else if(k<d){
  QuadRule SP,SP1,SP2;
  init_quadrule(&SP1, nr,1);
  init_quadrule(&SP2, nr,1);
  GJquad01(&SP1, nr, 0, d-k-1);
  GJquad01(&SP2, nr, 0, k-1);
  TensorQuad_nonit_ip(qp, SP1);
  init_quadrule(&SP, qp->n*SP2.n, qp->d+SP2.d);
  TensorQuad_nonit(&SP, *qp, SP2);
  SP=resize_quadrule_n(SP, 2);
  for(int i=0; i<SP.n/2; i++){
     SP.wt[i+SP.n/2]=SP.wt[i];
     SP.t[1][i]=SP.t[1][i]*SP.t[0][i];
     SP.t[2][i]=SP.t[2][i]*SP.t[0][i];
     SP.t[0][i+SP.n/2]=SP.t[1][i];
     SP.t[1][i+SP.n/2]=SP.t[0][i];
     SP.t[2][i+SP.n/2]=SP.t[2][i];
  }
  TensorQuad_nonit_ip(qp,SP1);
  *qp=resize_quadrule_n(*qp,3);
  for(int i=0; i<qp->n/3; i++){
      qp->t[1][i]=qp->t[1][i]*qp->t[0][i];
      qp->t[2][i]=qp->t[2][i]*qp->t[0][i];
      qp->wt[i+2*qp->n/3]=qp->wt[i];
      qp->t[0][i+2*qp->n/3]=qp->t[1][i];
      qp->t[1][i+2*qp->n/3]=qp->t[2][i];
      qp->t[2][i+2*qp->n/3]=qp->t[0][i];
      qp->wt[i]=SP.wt[i];
      qp->wt[i+qp->n/3]=SP.wt[i+qp->n/3];
      qp->t[0][i]=SP.t[0][i];
      qp->t[1][i]=SP.t[1][i];
      qp->t[2][i]=SP.t[2][i];
      qp->t[0][i+qp->n/3]=SP.t[0][i];
      qp->t[1][i+qp->n/3]=SP.t[1][i];
      qp->t[2][i+qp->n/3]=SP.t[2][i];
  }
  free_quadrule(SP); free_quadrule(SP1); free_quadrule(SP2);
}
if(d-k-1>0){
   QuadRule SP;
   SJacobi(&SP, nr,d-k-1);
   TensorQuad_nonit_ip(qp, SP);
   TensorQuad_nonit_ip(qp, SP);
   free_quadrule(SP);
   if(k==0){
      double* th;
      th=qp->t[s2];
      qp->t[s2]=qp->t[1];
      qp->t[1]=th;
   }
   else{
      if(k<d){
         if(d==k+2){
            double *th1,*th2;
            th1=qp->t[1];
            qp->t[1]=qp->t[3];
            th2=qp->t[2];
            qp->t[2]=th1;
            qp->t[3]=qp->t[4];
            qp->t[4]=th2;
         }
         else{
            double *th1, *th2;
            th1=qp->t[1];
            th2=qp->t[2];
            qp->t[1]=qp->t[d-k];
            qp->t[2]=qp->t[d-k+1];
            qp->t[d-k]=th1;
            qp->t[d-k+1]=th2;

            th1=qp->t[d-k+1];
            qp->t[d-k+1]=qp->t[2*d-2*k];
            qp->t[2*d-2*k]=th1;
         }
      }
   }
}
if(j>0){
  QuadRule SP;
  SJacobi(&SP, nr, j);
  TensorQuad_nonit_ip(qp, SP);
  free_quadrule(SP);
  double* th;
  th=qp->t[s3];
  qp->t[s3]=qp->t[2*d-2*k];
  qp->t[2*d-2*k]=th;
}

if(k-j-1>0){
   QuadRule SP;
   SJacobi(&SP, nr, k-j-1);
   TensorQuad_nonit_ip(qp, SP);
   free_quadrule(SP);
}
if(k>0){
   QuadRule SP;
   SJacobi(&SP, nr, k);
   TensorQuad_nonit_ip(qp, SP);
   free_quadrule(SP);
}

}

QuadRule QuadratureRule_gaujacsing(int k, int d,  int ns, int j, double alpha){
   int i;
   QuadRule QP;
   init_quadrule(&QP, ns, 1);
   GJquad01(&QP, ns, 0, 2*d-k-1+alpha);
   for(i=0; i<QP.n; i++){
       QP.wt[i]=QP.wt[i]/pow(QP.t[0][i],alpha+2*d-k-1);
   }
   return QP;
}

QuadRule QuadratureRule_gaujacreg(int k, int d, int nr, int j, int opt){
   //Compute quadrature rule on [0,1]^(2*d)
   QuadRule QP, SP, qp;
   if(k==-1){
      SJacobi(&QP, nr,d);
      init_quadrule(&SP, QP.n*QP.n, 2*d);
      TensorQuad_nonit(&SP, QP, QP);
      free_quadrule(QP);
   }
   else if(d-k>0){
      if(d-k-1>0){   //a
         SJacobi(&qp, nr, d-k-1);       //nr^d-k-1 x d-k-1
         init_quadrule(&QP, nr, 1);
         GJquad01(&QP, nr, 0, 0);
         init_quadrule(&SP, QP.n*qp.n, QP.d+qp.d);
         TensorQuad_nonit(&SP, qp, QP);
         TensorQuad_nonit_ip(&SP, qp);
         free_quadrule(QP); free_quadrule(qp);
      }
      else{
         init_quadrule(&SP, nr, 1);
         GJquad01(&SP, nr, 0,0);
      }
   }
   if(j>0){  //p
      if(k==d){   //SP has not been defined yet,
         SJacobi(&SP, nr, j);          //nr^j x j
      }
      else{
         SJacobi(&qp, nr, j);
         TensorQuad_nonit_ip(&SP, qp);
         free_quadrule(qp);
      }
   }
   if(k>0){  //s3
     if(k!=d || j!=0){
        if(k!=d){
           init_quadrule(&QP, nr, 1);
           if(opt==1||opt==2){
              GJquad01(&QP, nr,0, 0);
           }
           else{
              GJquad01(&QP, nr, 0,0);
           }
           TensorQuad_nonit_ip(&SP, QP);
           free_quadrule(QP);
        }
     }
   }
   if(k-j-1>0){  //q
      if(k!=d || j!=0){
         SJacobi(&qp, nr, k-j-1);   //nr^k-j-1 x k-j-1
         TensorQuad_nonit_ip(&SP, qp);
         free_quadrule(qp);
      }
      else{
         SJacobi(&SP, nr, k-j-1);
      }
   }
   if(k>0){     //ut
      SJacobi(&qp, nr, k);
      TensorQuad_nonit_ip(&SP, qp);
      free_quadrule(qp);
   }
   return SP;
}

QuadRule QuadratureRule_gaujac(int k, int d, int nr, int ns, int j, double alpha){
   //Compute quadrature rule on [0,1]^(2*d)
   int i;
   QuadRule QP, SP, qp;
   if(k==-1){
      SJacobi(&QP, nr,d);
      init_quadrule(&SP, QP.n*QP.n, 2*d);
      TensorQuad_nonit(&SP, QP, QP);
      free_quadrule(QP);
   }
   else if(d-k>0){
      init_quadrule(&QP, ns, 1);
      GJquad01(&QP, ns, 0, 2*d-k-1+alpha);
      for(i=0; i<QP.n; i++){
          QP.wt[i]=QP.wt[i]/pow(QP.t[0][i],alpha+2*d-k-1);
      }
      if(d-k-1>0){   //a
         SJacobi(&qp, nr, d-k-1);       //nr^d-k-1 x d-k-1
         init_quadrule(&SP, QP.n*qp.n, QP.d+qp.d);
         TensorQuad_nonit(&SP, QP, qp);
         free_quadrule(QP);
         init_quadrule(&QP, nr, 1);
         GLquad(&QP, nr, 0, 1);
         TensorQuad_nonit_ip(&SP, QP);
         TensorQuad_nonit_ip(&SP, qp);
         free_quadrule(QP); free_quadrule(qp);
      }
      else{
         init_quadrule(&qp, nr, 1);
         GLquad(&qp, nr, 0,1);
         init_quadrule(&SP, qp.n*QP.n, qp.d+QP.d);
         TensorQuad_nonit(&SP, QP, qp);
         free_quadrule(QP); free_quadrule(qp);
      }
   }
   if(j>0){  //p
      if(k==d){   //SP has not been defined yet,
         SJacobi(&SP, nr, j);          //nr^j x j
      }
      else{
         SJacobi(&qp, nr, j);
         TensorQuad_nonit_ip(&SP, qp);
         free_quadrule(qp);
      }
   }
   if(k>0){  //s3
     if(k==d && j==0){
          init_quadrule(&QP, ns, 1);
          GJquad01(&QP, ns, 0, 2*d-k-1+alpha);
          for(i=0; i<QP.n; i++){
              QP.wt[i]=QP.wt[i]/pow(QP.t[0][i],alpha+2*d-k-1);
          } 
     }
     else{
        if(k==d){ 
             init_quadrule(&QP, ns, 1);
             GJquad01(&QP, ns, 0, 2*d-k-1+alpha);
             for(i=0; i<QP.n; i++){
                 QP.wt[i]=QP.wt[i]/pow(QP.t[0][i],alpha+2*d-k-1);
              }
        }
        else{
           init_quadrule(&QP, nr, 1);
           GLquad(&QP, nr,0, 1);
        }
        TensorQuad_nonit_ip(&SP, QP);
        free_quadrule(QP);
     }
   }
   if(k-j-1>0){  //q
      SJacobi(&qp, nr, k-j-1);   //nr^k-j-1 x k-j-1
      TensorQuad_nonit_ip(&SP, qp);
      free_quadrule(qp);
   }
   if(k>0){     //ut
      SJacobi(&qp, nr, k);
      TensorQuad_nonit_ip(&SP, qp);
      free_quadrule(qp);
   }
   return SP;
}

QuadRule gauleg_reg(int d, int nr){
   //Compute quadrature rule on [0,1]^(2*d-1)
   int i,lx;
   int powh;
   QuadRule QP, qp;
   init_quadrule(&qp ,nr ,1);
   GLquad(&qp,nr, 0, 1);
   powh =(int) pow(nr, 2*d-1),
   init_quadrule(&QP, powh, 2*d-1);
   lx=nr;
   GLquad(&QP,nr, 0, 1);     //only regular

   for(i=1; i<2*d-1; i++){
      TensorQuad(i,lx, &QP, qp);
   }

   free_quadrule(qp);
   return QP;

}
