#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mxMemory.h"
#include "quadpoints.h"

void init_vtxlist(Vertexlist *vtx, int s1, int s2){
    int i;
   vtx->s1=s1;
   vtx->s2=s2;
   vtx->vtxlist=(double **)mxMalloc(vtx->s1*sizeof(double*));
   for(i=0; i<vtx->s1; i++){
      vtx->vtxlist[i] = (double*)mxMalloc(vtx->s2*sizeof(double));
   }
return;
}

void free_vtxlist(Vertexlist *vtx){
    int i;
   for(i=0; i<vtx->s1; i++){
      mxFree(vtx->vtxlist[i]);
   }
   mxFree(vtx->vtxlist);
return;
}


void resizequadpointsn(QuadraturePoints* quadpoints, int size){
   int i;
   quadpoints->n=size*quadpoints->n;

   for(i=0; i<quadpoints->d; i++){
      quadpoints->t[i]=(double*)mxRealloc(quadpoints->t[i], quadpoints->n*sizeof(double));
   }
   quadpoints->wt=(double*)mxRealloc(quadpoints->wt, quadpoints->n*sizeof(double));
return;
}

void resizequadpointsd(QuadraturePoints* quadpoints, int size){
   int i;
   int oldd=quadpoints->d;
   quadpoints->d=quadpoints->d+size;
   quadpoints->t=(double**)mxRealloc(quadpoints->t, quadpoints->d*sizeof(double*));
   for(i=oldd; i<quadpoints->d; i++){
       quadpoints->t[i]=(double*)mxCalloc(quadpoints->n,sizeof(double));
   }
return;
}

void init_quadpoints(QuadraturePoints* qp, int n, int d) {
   int i;
   qp->n = n;
   qp->d = d;
   qp->t = (double**)mxMalloc(qp->d*sizeof(double*));
   if(qp->t == NULL){
       /*mexPrintf("out of memory in init_quadpoints\n");*/
   }
   for(i=0; i<qp->d; i++){
      qp->t[i] = (double*)mxMalloc(qp->n*sizeof(double));
      if(qp->t[i] == NULL){
         /*mexPrintf("out of memory in init_quadpoints\n");*/
      }
   }
   qp->wt = (double*)mxMalloc(qp->n*sizeof(double));
return;
}

void free_qp(QuadraturePoints qp){
   int i;
   for(i=0; i<qp.d; i++){
      mxFree(qp.t[i]);
   }
   mxFree(qp.t);
   mxFree(qp.wt);
return;
}

AffineTrafo determineAffineTrafo(int d, int k, Vertexlist vtxlist){
   /* parameters of the affine transformation to the physical coordinates*/
   int i,j;
   int* idx2;

   AffineTrafo AfTr;
   AfTr.lvtx=vtxlist.s1;     AfTr.d=d;

   idx2=(int*)mxMalloc((d+1)*sizeof(int));

   for(i=0; i<k+1; i++){
      idx2[i]=i+1;
   }
   for(i=0; i<d-k; i++){
      idx2[k+1+i]=d+2+i;
   }

   AfTr.v10 = (double*)mxMalloc(vtxlist.s1*sizeof(double));
   AfTr.v20 = (double*)mxMalloc(vtxlist.s1*sizeof(double));

   for(i=0; i<vtxlist.s1; i++){
      AfTr.v10[i]=vtxlist.vtxlist[i][0];
      AfTr.v20[i]=vtxlist.vtxlist[i][idx2[0]-1];
   }

   AfTr.A1 = (double**)mxMalloc(vtxlist.s1*sizeof(double*));
   AfTr.A2 = (double**)mxMalloc(vtxlist.s1*sizeof(double*));

   for(i=0; i<vtxlist.s1; i++){
      AfTr.A1[i] = (double*)mxMalloc(d*sizeof(double));
      AfTr.A2[i] = (double*)mxMalloc(d*sizeof(double));
   }
   for(i=0; i<vtxlist.s1; i++){
      for(j=0; j<d; j++){
         AfTr.A1[i][j]=vtxlist.vtxlist[i][j+1]- AfTr.v10[i];
         AfTr.A2[i][j]=vtxlist.vtxlist[i][idx2[j+1]-1]-AfTr.v20[i];
      }
   }

   return AfTr;
}



void GLquad(QuadraturePoints* qp, int n, double a, double b){
/* Computes nodes x and weights w of Gauss_Legendre Quadrature Rule on [a,b]*/
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


void TensorQuad(QuadraturePoints* SP, QuadraturePoints QP, QuadraturePoints qp){
   /*Computes the tensor product of quadrature rules QP and qp*/
   int i, j, k;
   int lx,nx,ly,ny;

   lx=QP.n;   nx=QP.d;   ly=qp.n;   ny=qp.d;

   if(SP->n!=lx*ly || SP->d!=nx+ny){
       printf("SP has the wrong size.\n");
       return;
   }

   for(j=0; j<SP->n; j++){
         for(k=0; k<SP->d; k++){
            SP->t[k][j] = 0;
         }
         SP->wt[j]=0;
    }

   /*slow changing QP*/
   for(i=0; i<nx; i++){
      for(j=0; j<SP->n; j++){
            SP->t[i][j] = QP.t[i][j%lx];
      }
   }
   for(j=0; j<SP->n; j++){
      SP->wt[j]=QP.wt[j%lx];
   }
     
   /*fast changing qp*/
   for(i=0; i<lx; i++){
      for(k=0; k<ny; k++){
         for(j=0; j<ly; j++){
              SP->t[k+nx][j*lx+i]=qp.t[k][j];
          }
       }
   }
   
   for(i=0; i<lx; i++){
      for(j=0; j<ly;j++){
          SP->wt[j*lx+i]*=qp.wt[j];
      }
   }
return;
}

void QuadRule_cgl_gl(QuadraturePoints *SP, int nr, int ns, int s2){

   /*Compute quadrature rule on [0,1]^2*/
   QuadraturePoints QP, qp, SP2, SP3;

   init_quadpoints(&qp ,nr ,1);
   GLquad(&qp,nr, 0, 1);       /* regular*/

   if(s2==6){
       init_quadpoints(&QP ,nr ,1);
       GLquad(&QP,nr, 0, 1);       /* regular*/
   }
   else{
       init_quadpoints(&QP, ns*(ns+1)/2,1);
       CGLquad(&QP, ns);           /* singular*/
   }

   init_quadpoints(&SP2, qp.n*QP.n, 2);
   TensorQuad(&SP2, QP, qp);
   free_qp(QP);
   init_quadpoints(&SP3, qp.n*qp.n, 2);
   TensorQuad(&SP3, qp, qp);
   free_qp(qp);
   init_quadpoints(SP, SP3.n*SP2.n, 4);
   TensorQuad(SP, SP2, SP3);
   free_qp(SP2);
   free_qp(SP3);

}
