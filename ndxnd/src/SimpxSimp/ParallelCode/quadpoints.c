#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "quadpoints.h"

QuadraturePoints resizequadpoints(QuadraturePoints quadpoints, int size1, int size2){
   int i;
   quadpoints.n=size1;
   for(i=0; i<quadpoints.d; i++){
      quadpoints.t[i]=(double*)realloc(quadpoints.t[i], quadpoints.n*sizeof(double));
   }
   quadpoints.wt=realloc(quadpoints.wt, quadpoints.n*sizeof(double));
   int oldd=quadpoints.d;
   quadpoints.d=size2;
   quadpoints.t=(double**)realloc(quadpoints.t, quadpoints.d*sizeof(double*));
   for(i=oldd; i<quadpoints.d; i++){
       quadpoints.t[i]=(double*)malloc(quadpoints.n*sizeof(double));
   }
   return quadpoints;
}

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
   quadpoints.d=size*quadpoints.d;
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

void init_Quad2(Quad2* qp, int luv, int ruv, int lzh, int rzh) {
   int i;
   qp->luv = luv;
   qp->ruv = ruv;
   qp->lzh = lzh;
   qp->rzh = rzh;
   qp->u = (double**)malloc(qp->luv*sizeof(double*));
   for(i=0; i<qp->luv; i++){
      qp->u[i] = (double*)malloc(qp->ruv*sizeof(double));
   }
   qp->v = (double**)malloc(qp->luv*sizeof(double*));
   for(i=0; i<qp->luv; i++){
      qp->v[i] = (double*)malloc(qp->ruv*sizeof(double));
   }
   qp->zh = (double**)malloc(qp->lzh*sizeof(double*));
   for(i=0; i<qp->lzh; i++){
      qp->zh[i] = (double*)malloc(qp->rzh*sizeof(double));
   }
   qp->w = (double*)malloc(qp->ruv*sizeof(double));
}

void import_quadpoints(char* fname, QuadraturePoints* qp) {
    int i,n, d;
    FILE* f = fopen(fname, "rb");
    fread(&n, sizeof(int), 1, f);
    fread(&d, sizeof(int), 1, f);
    init_quadpoints(qp, n,d);
    for (i=0; i<d; i++) {
        fread(qp->t[i], sizeof(double), n, f);
    }
    fread(qp->wt, sizeof(double), n, f);
    fclose(f);
}

void export_quadpoints(char* fname, QuadraturePoints* qp) {
    int i;
    FILE* f = fopen(fname, "wb");
    fwrite(&qp->n, sizeof(int), 1, f);
    fwrite(&qp->d, sizeof(int), 1, f);
    for (i=0; i<qp->d; i++) {
        fwrite(qp->t[i], sizeof(double), qp->n, f);
    }
    fwrite(qp->wt, sizeof(double), qp->n, f);
    fclose(f);
}

void export_quad2(char* fname, Quad2* qp) {
    int i;
    FILE* f = fopen(fname, "wb");

    fwrite(&qp->luv, sizeof(int), 1, f);
    fwrite(&qp->ruv, sizeof(int), 1, f);

    for (i=0; i<qp->luv; i++) {
        fwrite(qp->u[i], sizeof(double), qp->ruv, f);
    }

    for (i=0; i<qp->luv; i++) {
        fwrite(qp->v[i], sizeof(double), qp->ruv, f);
    }

    fwrite(&qp->lzh, sizeof(int), 1, f);
    fwrite(&qp->rzh, sizeof(int), 1, f);

    for (i=0; i<qp->lzh; i++) {
        fwrite(qp->zh[i], sizeof(double), qp->rzh, f);
    }

    fwrite(qp->w, sizeof(double), qp->ruv, f);

    fclose(f);
}

void free_qp(QuadraturePoints qp){
   int i;
   for(i=0; i<qp.d; i++){
      free(qp.t[i]);
   }
   free(qp.t);
   free(qp.wt);
}

void free_Quad2(Quad2 qp){
   int i;
   for(i=0; i<qp.luv; i++){
      free(qp.u[i]);
      free(qp.v[i]);
   }
   free(qp.u);
   free(qp.v);
   for(i=0; i< qp.lzh; i++){
       free(qp.zh[i]);
   }
   free(qp.zh);
   free(qp.w);
}
