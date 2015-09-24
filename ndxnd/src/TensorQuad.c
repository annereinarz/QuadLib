#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stddef.h>

#include "QuadRule.h"
#include "TensorQuad.h"

#define eps DBL_EPSILON
#define NR_END 1
#define MAXIT 50
#define pi 3.141592653589793


void TensorQuad_new(QuadRule *SP,QuadRule QP, QuadRule qp){
	/*Computes the tensor product of quadrature rules QP and qp*/
	int i, j, k;
	int lx,nx,ly,ny;

	lx=QP.n;   nx=QP.d;   ly=qp.n;   ny=qp.d;

	init_quadrule(SP, lx*ly, nx+ny);

	if(SP->n!=lx*ly || SP->d!=nx+ny){
		printf("SP has the wrong size in ip.\n");
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
}

void TensorQuad_sparse(QuadRule *QP, QuadRule qp){
	/*Computes the tensor product of quadrature rules QP and qp*/
	int i, j, k;
	int lx,nx,ly,ny;

	lx=QP->n;   nx=QP->d;   ly=qp.n;   ny=qp.d;

	QuadRule SP;
	init_quadrule(&SP, lx*ly, nx+ny);

	if(SP.n!=lx*ly || SP.d!=nx+ny){
		printf("SP has the wrong size in ip.\n");
		return;
	}

	for(j=0; j<SP.n; j++){
		for(k=0; k<SP.d; k++){
		    SP.t[k][j] = 0;
		}
		SP.wt[j]=0;
	}

	/*slow changing QP*/
	for(i=0; i<nx; i++){
		for(j=0; j<SP.n; j++){
			SP.t[i][j] = QP->t[i][j%lx];
		}
	}
	for(j=0; j<SP.n; j++){
		SP.wt[j]=QP->wt[j%lx];
	}

	/*fast changing qp*/
	for(i=0; i<lx; i++){
		for(k=0; k<ny; k++){
			for(j=0; j<ly; j++){
				SP.t[k+nx][j*lx+i]=qp.t[k][j];
			}
		}
	}

	for(i=0; i<lx; i++){
		for(j=0; j<ly;j++){
			SP.wt[j*lx+i]*=qp.wt[j];
		}
	}

	*QP = resize_quadrule(*QP, lx*ly, nx+ny);
	for(i=0; i<QP->n; i++){
		for(j=0; j<QP->d; j++){
			QP->t[j][i]=SP.t[j][i];
		}
		QP->wt[i]=SP.wt[i];
	}
	free_quadrule(SP);
	return;
}

void TensorQuad(int it, int lx, QuadRule* QP, QuadRule qp){
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
   free(wh);
   for(i=0; i<nx; i++){
      free(th[i]);
   }
   free(th);
}

void TensorQuad_nonit(QuadRule* SP, QuadRule QP, QuadRule qp){
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


void TensorQuad_nonit_ip(QuadRule *QP, QuadRule qp){
	/*Computes the tensor product of quadrature rules QP and qp*/
	int i, j, k;
	int lx,nx,ly,ny;

	lx=QP->n;   nx=QP->d;   ly=qp.n;   ny=qp.d;

	QuadRule SP;
	init_quadrule(&SP, lx*ly, nx+ny);

	if(SP.n!=lx*ly || SP.d!=nx+ny){
		printf("SP has the wrong size in ip.\n");
		return;
	}

	for(j=0; j<SP.n; j++){
		for(k=0; k<SP.d; k++){
			SP.t[k][j] = 0;
		}
		SP.wt[j]=0;
	}

	/*slow changing QP*/
	for(i=0; i<nx; i++){
		for(j=0; j<SP.n; j++){
			SP.t[i][j] = QP->t[i][j%lx];
		}
	}
	for(j=0; j<SP.n; j++){
		SP.wt[j]=QP->wt[j%lx];
	}

	/*fast changing qp*/
	for(i=0; i<lx; i++){
		for(k=0; k<ny; k++){
			for(j=0; j<ly; j++){
				SP.t[k+nx][j*lx+i]=qp.t[k][j];
			}
		}
	}

	for(i=0; i<lx; i++){
		for(j=0; j<ly;j++){
			SP.wt[j*lx+i]*=qp.wt[j];
		}
	}

	*QP = resize_quadrule(*QP, lx*ly, nx+ny);
	for(i=0; i<QP->n; i++){
		for(j=0; j<QP->d; j++){
			QP->t[j][i]=SP.t[j][i];
		}
		QP->wt[i]=SP.wt[i];
	}
	free_quadrule(SP);
	return;
}
