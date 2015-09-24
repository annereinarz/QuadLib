#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "quadpoints.h"


QuadraturePoints resizequadpointsn(QuadraturePoints quadpoints, int size){
   int i;
   quadpoints.n=size*quadpoints.n;
   for(i=0; i<quadpoints.d; i++){
      quadpoints.t[i]=(double*)realloc(quadpoints.t[i], quadpoints.n*sizeof(double));
   }
   quadpoints.wt=realloc(quadpoints.wt, quadpoints.n*sizeof(double));
   return quadpoints;
}

QuadraturePoints resizequadpointsd(QuadraturePoints quadpoints, int size){
   int i;
   int oldd=quadpoints.d;
   quadpoints.d=size+quadpoints.d;
   quadpoints.t=(double**)realloc(quadpoints.t, quadpoints.d*sizeof(double*));
   for(i=oldd; i<quadpoints.d; i++){
       quadpoints.t[i]=(double*)malloc(quadpoints.n*sizeof(double));
   }
   return quadpoints;
}

void init_quadpoints(QuadraturePoints* qp, int n, int d) {
   int i;
   qp->n = n;
   qp->d = d;
   qp->t = (double**)malloc(qp->d*sizeof(double*));
   for(i=0; i<qp->d; i++){
      qp->t[i] = (double*)malloc(qp->n*sizeof(double));
   }
   qp->wt = (double*)malloc(qp->n*sizeof(double));
}

void free_qp(QuadraturePoints qp){
   int i;
   for(i=0; i<qp.d; i++){
      free(qp.t[i]);
   }
   free(qp.t);
   free(qp.wt);
}


