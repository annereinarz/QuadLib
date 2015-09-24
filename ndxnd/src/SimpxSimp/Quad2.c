#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Quad2.h"


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
