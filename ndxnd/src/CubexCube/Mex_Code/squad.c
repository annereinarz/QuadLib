#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "lapack.h"

#include "quadpoints.h"
#include "Tensorquad.h"
#include "functions.h"
#include "squad.h"

#include "mex.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))

double determinant(double**A, int m){
   int i, j;
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
   }/*end else*/
}

void Quadrature(int d, int k,int whichF, double* ptrQ, QuadraturePoints q){
   int i;
   double* F;
   int lengthF;
   /*  Quadrature Q=dot(F(u, v, z, alpha),w)  */
   if(whichF==1){
      F=Function1(q, d, k);
      lengthF=q.n;
   }
   else if(whichF==2){
      F=Function2(q, d, k);
      lengthF=q.n;
   }
   else if(whichF==3){
      F=Function3(q, d, k);
      lengthF=q.n;
   }
   else if(whichF==4){
      F=Function4(q, d, k);
      lengthF=q.n;
   }
   else{
      printf("invalid whichF, must be 1, 2, 3 or 4\n");
   }
   for(i=0; i<lengthF; i++){
      *ptrQ+=F[i]*q.wt[i];
   }
}

AffineTrafo determineAffineTrafo(int d, int k, Vertexlist vtxlist){
   /* parameters of the affine transformation to the physical coordinates   */
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
   return AfTr;
}

QuadraturePoints QuadratureRule(int k, int d, int nr, int ns){
   /*Compute quadrature rule on [0,1]^(2*d)*/
   int i,lx;
   int powh;
   QuadraturePoints QP, qp;
   init_quadpoints(&qp ,nr ,1);
   GLquad(&qp,nr, 0, 1);
   if(k==-1){
      powh =(int) pow(nr, 2*d),
      init_quadpoints(&QP, powh, 2*d);
      lx=nr;
      GLquad(&QP,nr, 0, 1);     /*only regular*/
   }
   else{
      lx=0.5*ns*(ns+1);
      powh=(int)pow(nr,2*d-1) ;
      init_quadpoints(&QP, lx*powh, 2*d);
      CGLquad(&QP, ns);         /*only singular*/
   }
   for(i=1; i<2*d; i++){
      TensorQuad(i,lx, &QP, qp);
   }

   free_qp(qp);
   return QP;
}  

int* FindSubset(int k, int numN){
   /* This routine generates the subset for the next Reflection */
   int remain;
   int* U;
   
   int i=0;
   U=(int*)calloc(k,sizeof(int)); 

   while(numN>0){
      remain=numN%2;
      numN=numN/2;
      U[i]=remain;
      i++;
   }
   return U;
}


QuadraturePoints Quad2PhyP(int k, int d, QuadraturePoints* QP, AffineTrafo A ){
   /* This routine realises step 1 of the algorithm,
    * transference to the physical parallelotope
    * zh in the Quad2 structure is necessary to avoid subtractive cancellation   */
   int i, j, r;
   QuadraturePoints sp;
   init_quadpoints(&sp, QP->n,3*d);

   if(k<=0){
      for(i=0; i< QP->n; i++){
         for(j=0; j<A.lvtx; j++){
            sp.t[j+2*d][i]=0;
            for(r=0; r<d; r++){
               sp.t[j+2*d][i]+=QP->t[r+d][i] * A.A2[j][r]  - QP->t[r][i]*A.A1[j][r];
            }
         }
      }
      if(k==-1){
         for(i=0; i< sp.n; i++){
            for(j=0; j<A.d; j++){
               sp.t[j+2*d][i]+=A.v20[j]-A.v10[j];
            }
         }
      }
   }
   else if(k<d){
      for(i=0; i< QP->n; i++){
         for(j=0; j<min(A.lvtx, QP->d); j++){
            sp.t[j+2*d][i]=0;
            for(r=0; r<k; r++){
               sp.t[j+2*d][i]+=QP->t[r+2*d][i] * A.A1[j][r];
            }
         }
      }
      for(i=0; i< QP->n; i++){
         for(j=0; j<A.lvtx; j++){
            for(r=max(k,0); r<d; r++){
               /*v*A2' - u*A1', wobei z, u, v mit Spalten zuerst gespeichert werden*/
               sp.t[j+2*d][i]+=QP->t[r+d][i] * A.A2[j][r]  - QP->t[r][i]*A.A1[j][r];
            }
         }
      }
   }
   else if(k==d){
      for(i=0; i< QP->n; i++){
         for(j=0; j<A.lvtx; j++){
            sp.t[j+2*d][i]=0;
            for(r=0; r<k; r++){
               sp.t[j+2*d][i]+=QP->t[r+2*d][i] * A.A1[j][r];
            }
         }
      }
   }
   for(i=0; i< QP->n; i++){
      for(j=0; j<A.d; j++){
         sp.t[j][i] = 0;
         sp.t[j+d][i] = 0;
         for(r=0; r<d; r++){
            sp.t[j][i]   += QP->t[r][i]   * A.A1[j][r];
            sp.t[j+d][i] += QP->t[r+d][i] * A.A2[j][r];
         }
         sp.t[j][i]   += A.v10[j];
         sp.t[d+j][i] += A.v20[j];
      }
   }
   double detdet=determinant(A.A1, d);
   detdet*=determinant(A.A2, d);
   detdet=fabs(detdet);
   for(i=0; i< QP->n; i++){
      sp.wt[i]=detdet*QP->wt[i];
   }
   return sp;
}


QuadraturePoints Step5(int j, int k, int d, QuadraturePoints* QP){
    int i,l,m;
    QuadraturePoints sp;
    init_quadpoints(&sp, QP->n, QP->d);
    for(i=0; i<sp.n; i++){
        sp.wt[i]=pow(QP->t[0][i], 2*d-k-1)*QP->wt[i];
    }
    /*Split into 2d-k pyramids*/
    for(i=0; i<2*d-k; i++){
        if(i<j){
            for(m=0; m<sp.n; m++){
              sp.t[i][m]=QP->t[0][m]*QP->t[i+1][m];
            }
        }
        else{
            for(m=0; m<sp.n; m++){
              sp.t[i][m]=QP->t[0][m]*QP->t[i][m];
            }
        }
        for(m=0; m<sp.n; m++){
          sp.t[j][m]=QP->t[0][m];
          for(l=0; l<k; l++){
              sp.t[2*d-k+l][m]=QP->t[2*d-k+l][m];
          }
        }
    }
    return sp;
}

void Step4(int numN, int k, QuadraturePoints* QP){
     int i,m, cnt=0;
     int* U;
     U=FindSubset(k, numN);
     /*mexPrintf("U=%d\n", U[0]);*/
     for(i=0; i<k; i++){
         if(U[i]!=0){
             for(m=0;m<QP->n; m++){
                 QP->t[cnt][m]=-QP->t[cnt][m];
             }
         }
         cnt++;
     }
     free(U);
}

void Step3(int numN, int k, int d, QuadraturePoints* QP){
    int i, m, cnt=0;
    int* U;
    U=FindSubset(k, numN);
    for(i=0; i<k; i++){
         if(U[i]!=0){
             for(m=0; m<QP->n; m++){
                 QP->t[2*d-k+i][m]=-QP->t[cnt][m]+(1-fabs(QP->t[i][m]))*QP->t[2*d-k+i][m];
             }
         }
         else{
             for(m=0; m<QP->n; m++){
                 QP->t[2*d-k+i][m]=                (1-fabs(QP->t[i][m]))*QP->t[2*d-k+i][m];
             }
         }
         for(m=0; m<QP->n; m++){
             QP->wt[m]=(1-fabs(QP->t[i][m]))*QP->wt[m];
         }
         cnt++;
     }
    free(U);
}

void Step2(int k, int d, QuadraturePoints* QP){
    int i,j;
    double ** s;

    if (QP->d != 2*d+k) {
        mexPrintf("Error in Step2: Wrong dimension.\n");
        exit(1);
    }

    for(i=0; i<QP->n; i++){
        for(j=0; j<k; j++){
           QP->t[2*d+j][i]=QP->t[j][i];
        }
    }

    for(i=0; i<QP->n; i++){
        for(j=0; j<k; j++){
           QP->t[j][i]=QP->t[j][i]+QP->t[2*d-k+j][i];
        }
    }

    /*t(:, d+1 : 2*d)=t(:,[(2*d-k+1 : 2*d)  (d+1 : 2*d-k)]);
     does the same as:
     s=t(:,d+1:2*d);
     t(:,d+1:d+k)=s(:,d-k+1:d);
     t(:, d+k+1:2*d)=s(:,1:d-k);*/

   s=(double**)malloc(d*sizeof(double*));

    for(j=0; j<d; j++){
        s[j]=QP->t[d+j];
    }
    for(j=0; j<k; j++){
        QP->t[d+j]=s[d-k+j];
    }
    for(j=0; j<d-k; j++){
        QP->t[d+k+j]=s[j];
    }
    free(s);
}

double squad(int nr, int ns, int whichvtxlist, int whichF){
/*QuadraturePoints squad(int nr, int ns, int whichvtxlist, int whichF){*/
    int j;
    int numN;
    double Q=0;

    Vertexlist vtxlist;
    vtxlist=get_vertexlist(whichvtxlist);

    /*compute space dimension d and dimension of the intersection k*/
    int d= vtxlist.s1;
    int k= 2*(d+1)-vtxlist.s2-1;
    
    /*Calculate Affine Transformation*/
    AffineTrafo A;
    A=determineAffineTrafo(d, k, vtxlist);

    /*Quadrature Rule on [0,1]^2*d*/
    QuadraturePoints QP, qp, sp;
    
    /*Initialise Quadrature Rule*/
    QP=QuadratureRule(k, d, nr, ns);

    if(k==-1){
        /*Transform to physical parallelotope */
        qp=Quad2PhyP(k, d, &QP, A );
        Quadrature(d, k, whichF, &Q, qp);
        free_qp(qp);
    }
    else{
        if(k==0){
            for(j=0; j<2*d; j++){
                /*Transform [0,1]^2d to D^N_j*/
                qp=Step5(j, k, d, &QP);
                /*Transform to physical parallelotope*/
                sp=Quad2PhyP(k, d, &qp, A );
                Quadrature(d, k, whichF, &Q, sp);
                free_qp(sp);
            }
        }
        else{
            for(numN=0; numN<pow(2,k); numN++){
                for(j=0; j<2*d-k; j++){
                   /*Transform [0,1]^2d to D^N_j*/
                   qp=Step5(j, k, d, &QP);
                   /*Reflections*/
                   Step4(numN, k, &qp);
                   /*utilde in [0,1]^k to uhat in F_N,M\N(zhat)*/
                   Step3(numN, k, d, &qp);
                   /*[-1,1]^k x F_k(zhat) to J^kxJ^k*/
                   qp=resizequadpointsd(qp,k);
                   Step2(k,d,&qp);
                   /*Transform to physical parallelotope*/
                   sp=Quad2PhyP(k, d, &qp, A);
                   Quadrature(d, k, whichF, &Q, sp);
                   free_qp(sp);
                }
            }
        }
    }      
    free_qp(QP);
    free_vertexlist(vtxlist);
    return Q;
}

Vertexlist get_vertexlist(int sim){
 int i;
 Vertexlist vtxlist;

 switch(sim){

      case 19:
         vtxlist.s1=2;
         vtxlist.s2=6;
         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }
         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         vtxlist.vtxlist[0][3]=-0.5;        vtxlist.vtxlist[1][3]=-0.5;
         vtxlist.vtxlist[0][4]=-1.5;         vtxlist.vtxlist[1][4]=-0.5;
         vtxlist.vtxlist[0][5]=-0.5;        vtxlist.vtxlist[1][5]=-1.5;
         break;

      case 20:
         vtxlist.s1=2;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         vtxlist.vtxlist[0][3]=-1;        vtxlist.vtxlist[1][3]=0;
         vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=-1;
         break;

      case 21:
         vtxlist.s1=2;
         vtxlist.s2=4;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=0;         vtxlist.vtxlist[1][1]=1;
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         vtxlist.vtxlist[0][3]=-1;         vtxlist.vtxlist[1][3]=0;
         break;

      case 22:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=0;         vtxlist.vtxlist[1][1]=1;
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         break;

     case 29:
         vtxlist.s1=3;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=0;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=1;
         vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=0;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
         vtxlist.vtxlist[0][4]=-1;   vtxlist.vtxlist[1][4]=-1;     vtxlist.vtxlist[2][4]=-1;
         vtxlist.vtxlist[0][5]=-2;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=-1;
         vtxlist.vtxlist[0][6]=-1;   vtxlist.vtxlist[1][6]=-2;      vtxlist.vtxlist[2][6]=-1;
         vtxlist.vtxlist[0][7]=-1;   vtxlist.vtxlist[1][7]=-1;     vtxlist.vtxlist[2][7]=-2;
         break;

       case 30:
         vtxlist.s1=3;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }
         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=0;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=1;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=0;
         vtxlist.vtxlist[0][4]=0;   vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=-1;
         vtxlist.vtxlist[0][5]=-1;    vtxlist.vtxlist[1][5]=0;     vtxlist.vtxlist[2][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;      vtxlist.vtxlist[2][6]=0;
         break;

       case 31:
         vtxlist.s1=3;
         vtxlist.s2=6;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=-1;     vtxlist.vtxlist[2][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=-1;
         break;

      case 32:
         vtxlist.s1=3;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=0;   vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=-1;
         break;

      case 33:
         vtxlist.s1=3;
         vtxlist.s2=4;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         break;

       /*case 39:
         vtxlist.s1=4;
         vtxlist.s2=10;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=-1;   vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=-1;    vtxlist.vtxlist[3][5]=-1;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;     vtxlist.vtxlist[2][6]=-1;    vtxlist.vtxlist[3][6]=-1;
         vtxlist.vtxlist[0][7]=-1;   vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=-1;     vtxlist.vtxlist[3][7]=-1;
         vtxlist.vtxlist[0][8]=-1;   vtxlist.vtxlist[1][8]=-1;     vtxlist.vtxlist[2][8]=0;     vtxlist.vtxlist[3][8]=-1;
         vtxlist.vtxlist[0][9]=-1;   vtxlist.vtxlist[1][9]=-1;     vtxlist.vtxlist[2][9]=-1;     vtxlist.vtxlist[3][9]=0;

          break;

     case 40:
         vtxlist.s1=4;
         vtxlist.s2=9;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=-1;   vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;     vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=-1;     vtxlist.vtxlist[3][7]=0;
         vtxlist.vtxlist[0][8]=0;    vtxlist.vtxlist[1][8]=0;      vtxlist.vtxlist[2][8]=0;     vtxlist.vtxlist[3][8]=-1;
         break;

      case 41:
         vtxlist.s1=4;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;    vtxlist.vtxlist[3][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=0;     vtxlist.vtxlist[3][7]=-1;
         break;

      case 42:
         vtxlist.s1=4;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=-1;    vtxlist.vtxlist[3][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=-1;
         break;

      case 43:
         vtxlist.s1=4;
         vtxlist.s2=6;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=-1;
         break;

      case 44:
         vtxlist.s1=4;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         break;   */
     default: mexPrintf("Invalid vertexlist chosen\n");
              exit(1);

   }
 return vtxlist;

}

 void free_vertexlist(Vertexlist vtxlist){
   int i;
   for(i=0; i<vtxlist.s1; i++){
       free(vtxlist.vtxlist[i]);
   }
   free(vtxlist.vtxlist);
}
