#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "quadpoints.h" 

void GLquad(QuadraturePoints* qp, int n, double a, double b){
  /*Computes nodes x and weights w of Gauss_Legendre Quadrature Rule on [a,b]*/
  int i,j;
  double z,z1=0,p1,p2,p3,pp;
  double pi=3.141592653589793;
  double eps=DBL_EPSILON;
  double m=(n+1)/2;

  /*Compute quadrature rule on [-1,1]*/
  for(i=1;i<=m; i++){
    z=cos(pi*(i-0.25)/(n+0.5));
    while(1){
      p1=1.0;
      p2=0.0;
      for(j=1; j<=n; j++){
        p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
      if(fabs(z-z1)<eps){
         break;
      }
    }
    qp->t[0][i-1]=-z;
    qp->t[0][n-i]=z;
    qp->wt[i-1]=2.0/((1.0-z*z)*pp*pp);
    qp->wt[n-i]=qp->wt[i-1];
  }
  /*Convert Quadrature rule to [a,b]*/
  for(i=0; i<n; i++){
     qp->t[0][i]=qp->t[0][i]*(b-a)/2+(b+a)/2;
     qp->wt[i]*=(b-a)/2.0;
  }
  return;
}

void CGLquad(QuadraturePoints* qp, int n){

  int i, j;
  double sigma=0.1;        /*parameter of geometric subdivision*/

  QuadraturePoints qpnew;

  double xl=sigma;
  double xr=1.0;
  int nj;
  int cnt=0;

  for(j=1; j<=n-1; j++){
      nj= n+1-j;
      init_quadpoints(&qpnew, nj, 1);
      GLquad(&qpnew,nj,xl,xr);
      for(i=1; i<=nj; i++){
          qp->t[0][cnt+i-1]=qpnew.t[0][i-1];
          qp->wt[cnt+i-1]=qpnew.wt[i-1];
      }
   
      free_qp(qpnew);
      xr=xl;
      xl=xl*sigma;
      cnt+=nj;
  }
  nj=1;
  init_quadpoints(&qpnew, nj, 1);
  GLquad(&qpnew,nj,0,xr);
  for(i=1; i<=nj; i++){
     qp->t[0][cnt+i-1]=qpnew.t[0][i-1];
     qp->wt[cnt+i-1]=qpnew.wt[i-1];
  }
  free_qp(qpnew);
  return;
}


void TensorQuad(int it, int lx, QuadraturePoints* QP, QuadraturePoints qp){
   /*Computes the tensor product of quadrature rules QP and qp*/
   int i, j, k;
   int ly=qp.n;     /*number of quadrature points in qp*/
   lx=lx*(int)pow(ly,it-1); /*number of quadrature points in QP*/
   int nx=it*qp.d;
   int ny=qp.d;
   double** th;
   th=(double**)malloc(nx*sizeof(double*));
   for(i=0; i<nx; i++){
      th[i]=(double*)malloc(lx*sizeof(double));
   }
   for(i=0; i<nx; i++){
      for(j=0; j<lx;j++){
         th[i][j]=QP->t[i][j];
      }
   }
   double* wh;
   wh = (double*)malloc(lx*sizeof(double));
   for(i=0; i<lx; i++){
      wh[i]=QP->wt[i];
   }
   /*slow changing QP*/
   for(i=0; i<ly; i++){
      for(j=0; j<lx; j++){
         for(k=0; k<nx; k++){
            QP->t[k][i*lx+j] = th[k][j];
         }
         QP->wt[i*lx+j]=wh[j];
      }
   }
   /*fast changing qp*/
   for(i=0; i<lx; i++){
      for(k=0; k<ny; k++){
         for(j=0; j<ly; j++){
              QP->t[k+nx][j*lx+i]=qp.t[k][j];
          }
       }
   }
   for(i=0; i<lx; i++){
      for(j=0; j<ly;j++){
          QP->wt[j*lx+i]*=qp.wt[j];
      }
   }
}
