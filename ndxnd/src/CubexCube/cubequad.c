#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "QuadRule.h"
#include "TensorQuad.h"
#include "functions.h"
#include "cubequad.h"
#include "sobol.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))


double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
}


void Quadrature(int d, int k,int whichF, double* ptrQ, QuadRule q, int whichoutput){
   int i;
   double* F;
   int lengthF;
   /*  Quadrature Q=dot(F(u, v, z, alpha),w)  */
   if(whichF==1){
	   if(whichoutput!=1){
		   printf("For this integrand whichoutput should be 1\n");
	   }
	   else {
           F=Function1(q, d, k);
           lengthF=q.n;
	   }
   }
   else if(whichF==11){
	   if(whichoutput!=7){
		   printf("For this integrand whichoutput should be 7\n");
	   }
	   else {
           F=Function11(q, d, k);
           lengthF=q.n;
	   }
   }
   else if(whichF==5){
	   if(whichoutput!=1){
		   printf("For this integrand whichoutput should be 1\n");
	   }
	   else {
           F=Function5(q, d, k);
           lengthF=q.n;
	   }
   }
   else if(whichF==2){
	   if(whichoutput!=6){
		   printf("For this integrand whichoutput should be 6\n");
	   }
	   else {
		  F=Function2(q, d, k);
		  lengthF=q.n;
	   }
   }
   else if(whichF==3){
	   if(whichoutput!=1){
		   printf("For this integrand whichoutput should be 1\n");
	   }
	   else {
		  F=Function3(q, d, k);
		  lengthF=q.n;
	   }
   }
   else if(whichF==4){
	   if(whichoutput!=6){
		   printf("For this integrand whichoutput should be 6\n");
	   }
	   else {
		  F=Function4(q, d, k);
		  lengthF=q.n;
	   }
   }
   else if(whichF==7){
	   if(whichoutput!=7){
		   printf("For this integrand whichoutput should be 7\n");
	   }
	   else {
		  F=Function7(q, d, k);
		  lengthF=q.n;
	   }
   }
   else{
      printf("invalid whichF, must be 1, 2, 3 or 4 or 5\n");
   }
   for(i=0; i<lengthF; i++){
      *ptrQ+=F[i]*q.wt[i];
   }
   free(F);
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


QuadRule Quad2PhyP(int k, int d, QuadRule* QP, AffineTrafo A , int whichoutput){
   /* This routine realises step 1 of the algorithm,
    * transference to the physical parallelotope
    * zh in the Quad2 structure is necessary to avoid subtractive cancellation   */
   int i, j, r;
   QuadRule sp;
   switch (whichoutput) {
      case 1:
      case 2:
      case 4:
         init_quadrule(&sp, QP->n, d);
         break;
      case 3:
      case 5:
      case 6:
         init_quadrule(&sp, QP->n, 2*d);
         break;
      case 7:
         init_quadrule(&sp, QP->n, 3*d);
         break;
      default:
         printf("invalid whichoutput\n");
         break;
   }
   int cnt=0;

   if (whichoutput!=1){
      for(i=0; i< QP->n; i++){
         for(j=0; j<A.d; j++){
            if(whichoutput & 2){
                sp.t[j][i] = 0;
                cnt=1;
            }
            if(whichoutput & 4){
                sp.t[j+cnt*d][i] = 0;
            }
            for(r=0; r<d; r++){
               if (whichoutput & 2){
                  sp.t[j][i]   += QP->t[r][i]   * A.A1[j][r];    // u
                  sp.t[j][i]   += A.v10[j];
               }
               if (whichoutput & 4){
                  sp.t[j+cnt*d][i] += QP->t[r+d][i] * A.A2[j][r];    // v
                  sp.t[cnt*d+j][i] += A.v20[j];
               }
            }
         }
      }
   }
   if(whichoutput==3||whichoutput==6||whichoutput==7){
      cnt=2;
   }

   if(k<=0){
      if (whichoutput & 1) {
         for(i=0; i< QP->n; i++){
            for(j=0; j<A.lvtx; j++){
               sp.t[j+cnt*d][i]=0;
               for(r=0; r<d; r++){
                  sp.t[j+cnt*d][i]+=QP->t[r+d][i] * A.A2[j][r]  - QP->t[r][i]*A.A1[j][r];
               }
            }
         }
         if(k==-1){
            for(i=0; i< sp.n; i++){
               for(j=0; j<A.d; j++){
                  sp.t[j+cnt*d][i]+=A.v20[j]-A.v10[j];
               }
            }
         }
      }
   }
   else {
      if (whichoutput & 1) {
         if(k<d){
            for(i=0; i< QP->n; i++){
               for(j=0; j<min(A.lvtx, QP->d); j++){
                  sp.t[j+cnt*d][i]=0;
                  for(r=0; r<k; r++){
                     sp.t[j+cnt*d][i]+=QP->t[r+2*d][i] * A.A1[j][r];   //z
                  }
               }
            }
            for(i=0; i< QP->n; i++){
               for(j=0; j<A.lvtx; j++){
                  for(r=max(k,0); r<d; r++){
                     /*v*A2' - u*A1', wobei z, u, v mit Spalten zuerst gespeichert werden*/
                     sp.t[j+cnt*d][i]+=QP->t[r+d][i] * A.A2[j][r]  - QP->t[r][i]*A.A1[j][r];   //z
                  }
               }
            }
         }
         else if(k==d){
            for(i=0; i< QP->n; i++){
               for(j=0; j<A.lvtx; j++){
                  sp.t[j+cnt*d][i]=0;
                  for(r=0; r<k; r++){
                     sp.t[j+cnt*d][i]+=QP->t[r+2*d][i] * A.A1[j][r];   //z
                  }
               }
            }
         }
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


 void Step5(int j, int k, int d, QuadRule* QP,QuadRule* sp){
    int i,l,m;
      for(i=0; i<sp->n; i++){
        sp->wt[i]=pow(QP->t[0][i], 2*d-k-1)*QP->wt[i];
      }
    /*Split into 2d-k pyramids*/
    for(i=0; i<2*d-k; i++){
        if(i<j){
            for(m=0; m<sp->n; m++){
              sp->t[i][m]=QP->t[0][m]*QP->t[i+1][m];
            }
        }
        else{
            for(m=0; m<sp->n; m++){
              sp->t[i][m]=QP->t[0][m]*QP->t[i][m];
            }
        }
        for(m=0; m<sp->n; m++){
          sp->t[j][m]=QP->t[0][m];
          for(l=0; l<k; l++){
              sp->t[2*d-k+l][m]=QP->t[2*d-k+l][m];
          }
        }
    }
}

void Step4(int numN, int k, QuadRule* QP){
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

void Step3(int numN, int k, int d, QuadRule* QP){
    int i, m, cnt=0;
    int* U;
    U=FindSubset(k, numN);
    for(i=0; i<k; i++){
         if(U[i]!=0){
             for(m=0; m<QP->n; m++){
                 QP->t[2*d-k+i][m]= QP->t[cnt][m]+(1-fabs(QP->t[i][m]))*QP->t[2*d-k+i][m];
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


void Step2(int k, int d, QuadRule* QP){
    int i,j;
    double ** s;

    if (QP->d != 2*d+k) {
        printf("Error in Step2: Wrong dimension.\n");
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

void cubetransform(int k, int d, FILE *f, QuadRule *QP){
   int j;
   int numN;

   /*Quadrature Rule on [0,1]^2*d*/
   QuadRule  qp;

   if(k==-1){
      // save QP in a file
      write_quadrule(f,QP);
   }
   else{
      if(k==0){
         for(j=0; j<2*d; j++){
            /*Transform [0,1]^2d to D^N_j*/
            init_quadrule(&qp, QP->n, QP->d);
            Step5(j, k, d, QP, &qp);
            // save qp in file
            write_quadrule(f,&qp);
            free_quadrule(qp);
         }
      }
      else{
         for(numN=0; numN<pow(2,k); numN++){
            for(j=0; j<2*d-k; j++){
               qp=BasicSubQuad(j,k,d,numN,QP);
               /*utilde in [0,1]^k to uhat in F_N,M\N(zhat)*/
               Step3(numN, k, d, &qp);
               /*[-1,1]^k x F_k(zhat) to J^kxJ^k*/
               qp=resize_quadrule_d(qp,k);
               Step2(k,d,&qp);
               // save qp in file
               write_quadrule(f,&qp);
               free_quadrule(qp);
            }
         }
      }
   }
}

void cubeaffine(int k, int d, AffineTrafo A, int whichoutput, QuadRule *QP, FILE *f){
   int j;
   int numN;

   /*Quadrature Rule on [0,1]^2*d*/
   QuadRule  qp, sp;

   if(k==-1){
      /*Transform to physical parallelotope */
      qp=Quad2PhyP(k, d, QP, A , whichoutput);
      // save qp in a file
      write_quadrule(f,&qp);
      free_quadrule(qp);
   }
   else{
      if(k==0){
         for(j=0; j<2*d; j++){
            /*Transform [0,1]^2d to D^N_j*/
            init_quadrule(&qp, QP->n, QP->d);
            Step5(j, k, d, QP, &qp);
            /*Transform to physical parallelotope*/
            sp=Quad2PhyP(k, d, &qp, A , whichoutput);
            free_quadrule(qp);
            // save sp in file
            write_quadrule(f,&sp);
            free_quadrule(sp);
         }
      }
      else{
         for(numN=0; numN<pow(2,k); numN++){
            for(j=0; j<2*d-k; j++){
               qp=BasicSubQuad(j,k,d,numN,QP);
               /*utilde in [0,1]^k to uhat in F_N,M\N(zhat)*/
               Step3(numN, k, d, &qp);
               /*[-1,1]^k x F_k(zhat) to J^kxJ^k*/
               qp=resize_quadrule_d(qp,k);
               Step2(k,d,&qp);
               /*Transform to physical parallelotope*/
               sp=Quad2PhyP(k, d, &qp, A, whichoutput);
               free_quadrule(qp);
               // save sp in file
               write_quadrule(f,&sp);
               free_quadrule(sp);
            }
         }
      }
   }
}

double cubequad(int k, int d, AffineTrafo A, int whichF, int whichoutput, int* ndof, QuadRule *QP){
   int j;
   int numN;
   double Q=0;

   /*Quadrature Rule on [0,1]^2*d*/
   QuadRule  qp, sp;

   if(k==-1){
      /*Transform to physical parallelotope */
      qp=Quad2PhyP(k, d, QP, A , whichoutput);
      *ndof+=QP->n;
      Quadrature(d, k, whichF, &Q, qp, whichoutput);
      free_quadrule(qp);
   }
   else{
      if(k==0){
         for(j=0; j<2*d; j++){
            /*Transform [0,1]^2d to D^N_j*/
            init_quadrule(&qp, QP->n, QP->d);
            Step5(j, k, d, QP, &qp);
            /*Transform to physical parallelotope*/
            sp=Quad2PhyP(k, d, &qp, A , whichoutput);
            free_quadrule(qp);
            *ndof+=sp.n;
            Quadrature(d, k, whichF, &Q, sp, whichoutput);
            free_quadrule(sp);
         }
      }
      else{
         for(numN=0; numN<pow(2,k); numN++){
            for(j=0; j<2*d-k; j++){
               qp=BasicSubQuad(j,k,d,numN,QP);
               /*utilde in [0,1]^k to uhat in F_N,M\N(zhat)*/
               Step3(numN, k, d, &qp);
               /*[-1,1]^k x F_k(zhat) to J^kxJ^k*/
               qp=resize_quadrule_d(qp,k);
               Step2(k,d,&qp);
               /*Transform to physical parallelotope*/
               sp=Quad2PhyP(k, d, &qp, A, whichoutput);
               free_quadrule(qp);
               *ndof+=sp.n;
               Quadrature(d, k, whichF, &Q, sp, whichoutput);
               free_quadrule(sp);
            }
         }
      }
   }
   return Q;
}


QuadRule BasicSubQuad(int j,int k,int d,int numN, QuadRule* QP){
  QuadRule qp;
  /*Transform [0,1]^2d to D^N_j*/
  init_quadrule(&qp, QP->n, QP->d);
  Step5(j, k, d, QP, &qp);
  /*Reflections*/
  Step4(numN, k, &qp);
  return qp;
}


