#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "clapack.h"

#include "AffineTrafo.h"
#include "quadpoints.h"
#include "Tensorquad.h"
#include "functions.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))

int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

void PermRefl(int k, int d, QuadraturePoints* QP, AffineTrafo A, int* P, double* Q, int flag, int whichF){
   int i, j, r;
   Quad2 qp2;
   if(flag==0){
      // Step 1 of the algorithm,
      // transference to the physical simplex
      init_Quad2(&qp2, d, QP->n, A.lvtx, QP->n);
      for(i=0; i< QP->n; i++){
         for(j=0; j<A.lvtx; j++){
            qp2.zh[j][i]=0;
            for(r=0; r<d; r++){
               qp2.zh[j][i]+=QP->t[d+r][i] * A.A2[j][r]  - QP->t[r][i]*A.A1[j][r];
            }
         }
      }
      if(k==-1){
         for(i=0; i< qp2.rzh; i++){
            for(j=0; j<A.d; j++){
               qp2.zh[j][i]+=A.v20[j]-A.v10[j];
            }
         }
      }
      for(i=0; i< QP->n; i++){
         for(j=0; j<A.d; j++){
            qp2.u[j][i] = 0;
            qp2.v[j][i] = 0;
            for(r=0; r<d; r++){
               qp2.u[j][i] += QP->t[r][i] * A.A1[j][r];
               qp2.v[j][i] += QP->t[d+r][i] * A.A2[j][r];
            }
            qp2.u[j][i] += A.v10[j];
            qp2.v[j][i] += A.v20[j];
         }
      }
      double detdet=determinant(A.A1, d);
      detdet*=determinant(A.A2, d);
      detdet=fabs(detdet);
      for(i=0; i< qp2.ruv; i++){
         qp2.w[i]=detdet*QP->wt[i];
      }
      Quadrature(d,k, whichF, Q, qp2);
      free_Quad2(qp2);
   }

   else{
      if(flag==1){     //permutations
         Quad2 q2;
         init_Quad2(&q2, d, QP->n,k, QP->n );
         for(i=0; i<k; i++){
            for(r=0; r<q2.ruv; r++){
               q2.u[i][r]=QP->t[P[2*d-k+i]-1][r];
            }
         }
         for(i=0; i<d-k; i++){
            for(r=0; r<q2.ruv; r++){
               q2.u[k+i][r]=QP->t[P[i]-1][r];
            }
         }
         for(i=0; i<k; i++){
            for(r=0; r<q2.rzh; r++){
               q2.zh[i][r]=QP->t[P[2*(d-k)+i]-1][r];
            }
         }
         for(i=0; i<k; i++){
            for(r=0; r<q2.ruv; r++){
               q2.v[i][r]=q2.u[i][r]+q2.zh[i][r];
            }
         }
         for(i=0; i<d-k; i++){
            for(r=0; r<q2.ruv; r++){
               q2.v[k+i][r]=QP->t[P[d-k+i]-1][r];
            }
         }
         for(i=0; i<q2.ruv; i++){
            q2.w[i]=QP->wt[i];
         }
         Quad2RefS(k, d, &q2);
         Quad2 qp2 = Quad2PhyS(k,d,q2, A);
         free_Quad2(q2);
         Quadrature(d,k, whichF, Q, qp2);
         free_Quad2(qp2);
      }
      if(flag==2){         //reflections
         Quad2 q3;
         init_Quad2(&q3,d, QP->n,k, QP->n);
         for(i=0; i<k; i++){
            for(r=0; r<q3.ruv; r++){
               q3.u[i][r]=QP->t[P[2*d-k+i]-1][r]+QP->t[P[2*(d-k)+i]-1][r];
            }
         }
         for(i=0; i<d-k; i++){
            for(r=0; r<q3.ruv; r++){
               q3.u[k+i][r]=QP->t[P[i]-1][r];
            }
         }
         for(i=0; i<k; i++){
            for(r=0; r<q3.rzh; r++){
               q3.zh[i][r]=-QP->t[P[2*(d-k)+i]-1][r];
            }
         }
         for(i=0; i<k; i++){
            for(r=0; r<q3.ruv; r++){
               q3.v[i][r]=QP->t[P[2*d-k+i]-1][r];
            }
         }
         for(i=0; i<d-k; i++){
            for(r=0; r<q3.ruv; r++){
               q3.v[k+i][r]=QP->t[P[d-k+i]-1][r];
            }
         }
         for(i=0; i<q3.ruv; i++){
            q3.w[i]=QP->wt[i];
         }
         Quad2RefS(k,d,&q3);
         Quad2 qp2=Quad2PhyS(k,d,q3, A);
         free_Quad2(q3);
         Quadrature(d,k, whichF, Q, qp2);
         free_Quad2(qp2);
      }
   }
}


void free_afftrafo(AffineTrafo AfTr){
   int i;
   free(AfTr.v10);
   free(AfTr.v20);
   
   for(i=0; i<AfTr.lvtx; i++){
      free(AfTr.A1[i]);
      free(AfTr.A2[i]);
   }
   free(AfTr.A1);
   free(AfTr.A2);
}

double determinant(double**A, int m){
   int i, j;
   //if(m==3){
     //double det=A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])-A[1][0]*(A[0][1]*A[2][2]-A[2][1]*A[0][2])+A[2][0]*(A[0][1]*A[1][2]-A[1][1]*A[0][2]);
       /*for(i=0; i<3; i++){
          for(j=0; j<3; j++){
               printf("%lf ", A[i][j]);
           }
           printf("\n");
       }*/
       //printf("det=%lf\n", det);
     //return det;
   //}
   //else{
   double Ah[m*m];
   for(i=0; i<m; i++){
      for(j=0; j<m; j++){
         Ah[i+j*m]=A[i][j];
      }
   }
   int IPEV[m*m];
   int INFO;
   dgetrf_(&m ,&m ,Ah ,&m ,IPEV, &INFO);
   if(INFO==0){
      double det=1;
      for(i=0; i<m; i++){
         det*=Ah[i+i*m];
      }
      return det;
   }
   else{
      if(INFO<0){
         printf("A(%d, %d) is illegal in determinant function\n", INFO,INFO);
      }
      else{
         printf("A(%d, %d) is exactly zero in determinant function, the determinant is zero\n", INFO,INFO);
      }
         return 0;
      }//end else
   //}
}




void Quadrature(int d, int k,int whichF, double* ptrQ, Quad2 q){
   int i;
   double* F;
   int lengthF;
   //  Quadrature Q=dot(F(u, v, z, alpha),w)
   if(whichF==1){
      F=Function1(q, d, k);
      lengthF=q.rzh;
   }
   else if(whichF==2){
      F=Function2(q, d, k);
      lengthF=q.ruv;
   }
   else if(whichF==3){
      F=Function3(q, d, k);
      lengthF=q.rzh;
   }
   else if(whichF==4){
      F=Function4(q, d, k);
      lengthF=q.ruv;
   }
   else if(whichF==5){
      F=Function5(q,d,k);
      lengthF=q.rzh;
   }
   else{
      printf("invalid whichF, must be 1, 2, 3 or 4\n");
   }
   /*for(i=0; i<lengthF; i++){
      printf("F[%d]=%lf\n",i, F[i]);
   }*/
   for(i=0; i<lengthF; i++){
      *ptrQ+=F[i]*q.w[i];
   }
   free(F);
}

AffineTrafo determineAffineTrafo(int d, int k, Vertexlist vtxlist){
   // parameters of the affine transformation to the physical coordinates
   int i,j;
   AffineTrafo AfTr;
   AfTr.lvtx=vtxlist.s1;     AfTr.d=d;
   int* idx2;
   idx2=(int*)malloc((d+1)*sizeof(int));
   for(i=0; i<k+1; i++){
      idx2[i]=i+1;
   }
   for(i=0; i<d-k; i++){
      idx2[k+1+i]=d+2+i;
   }

   AfTr.v10 = (double*)malloc(vtxlist.s1*sizeof(double));
   AfTr.v20 = (double*)malloc(vtxlist.s1*sizeof(double));
   for(i=0; i<vtxlist.s1; i++){
      AfTr.v10[i]=vtxlist.vtxlist[i][0];
      AfTr.v20[i]=vtxlist.vtxlist[i][idx2[0]-1];
   }
   AfTr.A1 = (double**)malloc(vtxlist.s1*sizeof(double*));
   AfTr.A2 = (double**)malloc(vtxlist.s1*sizeof(double*));
   for(i=0; i<vtxlist.s1; i++){
      AfTr.A1[i] = (double*)malloc(d*sizeof(double));
      AfTr.A2[i] = (double*)malloc(d*sizeof(double));
   }
   for(i=0; i<vtxlist.s1; i++){
      for(j=0; j<d; j++){
         AfTr.A1[i][j]=vtxlist.vtxlist[i][j+1]- AfTr.v10[i];
         AfTr.A2[i][j]=vtxlist.vtxlist[i][idx2[j+1]-1]-AfTr.v20[i];
      }
   }
   free(idx2);
   return AfTr;
}


 int NextPerm(int j, int k, int d, int* U){
   // This routine generates the next permutation of U in the lexiographic order
   int i,m;
   int* P = (int*)malloc((2*d)*sizeof(int));
   for(i=0; i<2*d; i++){
      P[i]=i+1;
   }
   for(i=0; i<2*d; i++){
      P[U[i]-1]=i+1;    //back permutation
   }
   for(i=0; i<j; i++){
      U[i]=P[2*(d-k)+i] - 2*(d-k);
   }

   // compute the subpermutation u
   int flag=0;
   for(i=j-1; i>=0; i--){
      if(U[i]<k-j+i+1){
         U[i] += 1;
         for(m=i+1; m<j; m++){
            U[m]=U[m-1]+1;
         }
         flag=1;
         break;
      }
   }
   if(flag==0){
      free(P);
      return 3;
   }
   //compute full permutation P
   for(i=0; i<k; i++){
      P[i]=i+1;
   }
   for(i=0; i<j; i++){
      P[U[i]-1]=0;
   }

   int cnt=0;
   int* w=(int*)malloc((k-j)*sizeof(int));
   for(i=0; i<k; i++){
      if(P[i]!=0){
         w[cnt]=P[i];
         cnt++;
      }
   }
   for(i=0; i<2*(d-k); i++){
      P[i]=i+1;
   }

   for(i=2*(d-k); i<2*d-k; i++){
      if(i<2*(d-k)+j){
         P[i]=2*(d-k)+U[i-2*(d-k)];
      }
      else{
         P[i]=2*(d-k)+w[i-j-2*(d-k)];
      }
   }
   for(i=2*d-k; i<2*d; i++){
      if(i<2*d-k+j){
         P[i]=2*d-k+U[i-2*d+k];
      }
      else{
         P[i]=2*d-k+w[i-2*d+k-j];
      }
   }
   free(w);
   for(i=0; i<2*d; i++){
      U[i]=i+1;
   }
   for(i=0; i<2*d; i++){
      U[P[i]-1]=i+1;      //back permutation
   }
   free(P);
   return 0;
}




void Quad2RefS(int k, int d,  Quad2* qp){
   // This function realises Step 2 of the algorithm,
   // transference to the reference simplex
   // zh=\hat z in the Quad2 structure is necessary to avoid subtractive cancellation
   int i, m;
   double *vdelta, *udelta;
   vdelta = (double*)calloc(qp->ruv, sizeof(double));
   udelta = (double*)calloc(qp->ruv, sizeof(double));

   for(i=0; i<qp->ruv; i++){
      for(m=max(k,0); m<d; m++){
         udelta[i]+=qp->u[m][i];
         vdelta[i]+=qp->v[m][i];
      }
   }
   for(i=0; i<qp->ruv; i++){
      for(m=0; m<k; m++){
         qp->zh[m][i]=qp->zh[m][i]*(1-vdelta[i])+qp->u[m][i]*(udelta[i]-vdelta[i]);
      }
   }

   for(i=0; i<qp->ruv; i++){
      udelta[i]=1-udelta[i];
      vdelta[i]=1-vdelta[i];
   }
   for(i=0; i< qp->ruv; i++){
      for(m=0; m<k; m++){
         qp->u[m][i]*=udelta[i];
         qp->v[m][i]*=vdelta[i];
      }
   }

   for(i=0; i<qp->ruv; i++){
      qp->w[i]*= pow(udelta[i]*vdelta[i], k);
   }
   free(udelta);
   free(vdelta);
}

Quad2 Quad2PhyS(int k, int d, Quad2 qp, AffineTrafo A ){
   // This routine realises step 1 of the algorithm,
   // transference to the physical simplex
   // zh in the Quad2 structure is necessary to avoid subtractive cancellation
   int i, j, r;
   Quad2 qp2;
   init_Quad2(&qp2, qp.luv, qp.ruv, d, qp.ruv);

   if(k<=0){
      for(i=0; i< qp.ruv; i++){
         for(j=0; j<d; j++){
            qp2.zh[j][i]=0;
            for(r=0; r<d; r++){
               qp2.zh[j][i]+=qp.v[r][i] * A.A2[j][r]  - qp.u[r][i]*A.A1[j][r];
            }
         }
      }
      if(k==-1){
         for(i=0; i< qp2.rzh; i++){
            for(j=0; j<d; j++){
               qp2.zh[j][i]+=A.v20[j]-A.v10[j];
            }
         }
      }
   }
   else if(k<d){
      for(i=0; i< qp.rzh; i++){
         for(j=0; j<d; j++){
            qp2.zh[j][i]=0;
            for(r=0; r<k; r++){
               qp2.zh[j][i]+=qp.zh[r][i] * A.A1[j][r];
            }
         }
      }
      for(i=0; i< qp.ruv; i++){
         for(j=0; j<d; j++){
            for(r=max(k,0); r<d; r++){
               // v * A2' - u * A1', wobei z, u, v mit Spalten zuerst gespeichert werden
               qp2.zh[j][i]+=qp.v[r][i]*A.A2[j][r]  - qp.u[r][i]*A.A1[j][r];
            }
         }
      }
   }
   else if(k==d){
      for(i=0; i< qp.ruv; i++){
         for(j=0; j<d; j++){
            qp2.zh[j][i]=0;
            for(r=0; r<k; r++){
               qp2.zh[j][i]+=qp.zh[r][i]*A.A1[j][r];
            }
         }
      }
   }
   for(i=0; i< qp.ruv; i++){
      for(j=0; j<A.d; j++){
         qp2.u[j][i] = 0;
         qp2.v[j][i] = 0;
         for(r=0; r<d; r++){
            qp2.u[j][i] += qp.u[r][i] * A.A1[j][r];
            qp2.v[j][i] += qp.v[r][i] * A.A2[j][r];
         }
         qp2.u[j][i] += A.v10[j];
         qp2.v[j][i] += A.v20[j];
      }
   }
   double detdet=determinant(A.A1, d);
   detdet*=determinant(A.A2, d);
   detdet=fabs(detdet);
   for(i=0; i< qp2.ruv; i++){
      qp2.w[i]=detdet*qp.w[i];
   }

   return qp2;
}
