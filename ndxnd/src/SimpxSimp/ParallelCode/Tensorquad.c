#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "quadpoints.h"
#include "Tensorquad.h"

#include "clapack.h"

#define eps DBL_EPSILON
#define pi 3.141592653589793

double tgamma(double x);
int dsyev_ (char *jobz, char *uplo, int *n, double *fa, int *lda, double *w, double *work, int *lwork, int *info);

void GJquad01(QuadraturePoints* qp, int n, double alpha, double beta){
	int i;
	GJquad(qp,n,alpha,beta);
	for(i=0; i<n; i++){
		qp->t[0][i]=(qp->t[0][i]+1)/2;
	}
}

void GJquad(QuadraturePoints* qp, int n, double alpha, double beta){
	/*Given alpha and beta, the parameters of the Jacobi polynomials, this routine returns
	  arrays x and w containing the abscissas and weights of the n-point Gauss-Jacobi
	  quadrature formula.*/
	int i;

	if (n==1){
		qp->t[0][0]=(alpha-beta)/(alpha+beta+2);
		qp->wt[0] = 2;
		return;
	}

	/*Compute J*/
	double * J;
	J=(double*)calloc(n*n, sizeof(double));

	double *h1;
	h1=(double*)malloc(n*sizeof(double));
	for(i=0; i<=n-1; i++){
		h1[i]=2*i+alpha+beta;
	}

	for(i=0; i<n; i++){
		J[i+i*n]=-(pow(alpha,2.0)-pow(beta,2.0))/(h1[i]+2)/h1[i];
	}
	for(i=0; i<n-1; i++){
		J[i+(i+1)*n]=2.0/(h1[i]+2)*sqrt( (i+1)*((i+1)+alpha+beta)*
				((i+1)+alpha)*((i+1)+beta)/(h1[i]+1)/(h1[i]+3) );
	}
	if (fabs(alpha+beta)<10*eps){
		J[0]=0.0;
	}
   free(h1);
	/*{int j;
	  for(i=0; i<n; i++){
	  for(j=0; j<n; j++){
	  printf("%lf ",J[i+j*n]);
	  }printf("\n");
	  }
	  }*/
	/*Compute quadrature nodes and weights*/
	int info;

	double dbl_lwork;
	int lwork=-1;

	dsyev_("V", "U", &n, J, &n, qp->t[0], &dbl_lwork, &lwork, &info);
	if (0 != info) {
		printf("Error using dsyev_");
		exit(1);
	}

	lwork=(int)dbl_lwork;
	double *WORK;
	WORK=(double*)malloc(lwork*sizeof(double));

	dsyev_("V", "U", &n, J, &n, qp->t[0], WORK, &lwork, &info);

	free(WORK);

	/*{int j;
	  for(i=0; i<n; i++){
	  for(j=0; j<n; j++){
	  printf("%lf ",J[i+j*n]);
	  }printf("\n");
	  }}*/
	int j;
	double sum;
	for(i=0; i<n; i++){
		sum=0;
		for(j=0; j<n; j++){  /*used to normalise the eigenvector*/
			sum+=pow(J[i*n+j],2.0);
		}
		qp->wt[i]=pow(J[i*n]/sqrt(sum),2.0)/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);
	}

	free(J);
}

void GLquad(QuadraturePoints* qp, int n, double a, double b){
	// Computes nodes x and weights w of Gauss_Legendre Quadrature Rule on [a,b]
	int i,j;
	double z,z1=0,p1,p2,p3,pp;
	double m=(n+1)/2;

	//Compute quadrature rule on [-1,1]
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
	//Convert Quadrature rule to [a,b]
	for(i=0; i<n; i++){
		qp->t[0][i]=qp->t[0][i]*(b-a)/2+(b+a)/2;
		qp->wt[i]*=(b-a)/2.0;
	}
}

void CGLquad(QuadraturePoints* qp, int n){

	int i, j;
	double sigma=0.1;        //parameter of geometric subdivision

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
			qp->t[0][cnt+(i-1)]=qpnew.t[0][i-1];
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
}




void TensorQuad(int it, int lx, QuadraturePoints* QP, QuadraturePoints qp){
	// Computes the tensor product of quadrature rules QP and qp
	int i, j, k;
	int ly=qp.n;     //number of quadrature points in qp
	lx=lx*(int)pow(ly,it-1); //number of quadrature points in QP
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
	//slow changing QP
	for(i=0; i<ly; i++){
		for(j=0; j<lx; j++){
			for(k=0; k<nx; k++){
				QP->t[k][i*lx+j] = th[k][j];
			}
			QP->wt[i*lx+j]=wh[j];
		}
	}
	for(i=0; i<nx; i++){
      free(th[i]);
	}
   free(th);
   free(wh);
	//fast changing qp
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

void TensorQuad_nonit(QuadraturePoints* SP, QuadraturePoints QP, QuadraturePoints qp){
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


void TensorQuad_nonit_ip(QuadraturePoints *QP, QuadraturePoints qp){
	/*Computes the tensor product of quadrature rules QP and qp*/
	int i, j, k;
	int lx,nx,ly,ny;

	lx=QP->n;   nx=QP->d;   ly=qp.n;   ny=qp.d;

	QuadraturePoints SP;
	init_quadpoints(&SP, lx*ly, nx+ny);

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

	*QP = resizequadpoints(*QP, lx*ly, nx+ny);
	for(i=0; i<QP->n; i++){
		for(j=0; j<QP->d; j++){
			QP->t[j][i]=SP.t[j][i];
		}
		QP->wt[i]=SP.wt[i];
	}
	free_qp(SP);
	return;
}




