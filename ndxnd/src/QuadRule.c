#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "QuadRule.h"

void write_quadrule(FILE* f, QuadRule *qp){
   int i;
   fwrite(&qp->n, sizeof(int), 1, f);
   fwrite(&qp->d, sizeof(int), 1, f);
   for (i=0; i<qp->d; i++) {
       fwrite(qp->t[i], sizeof(double), qp->n, f);
   }
   fwrite(qp->wt, sizeof(double), qp->n, f);
   fflush(f);
}

void write_1d_quadrule(FILE* f, QuadRule *qp){
   fwrite(&qp->n, sizeof(int), 1, f);
   fwrite(qp->t,  sizeof(double), qp->n, f);
   fwrite(qp->wt, sizeof(double), qp->n, f);
   fflush(f);
 }

void import_1d_quadrule(char* fname, QuadRule* qp) {
    int n;
    FILE* f = fopen(fname, "rb");
    fread(&n, sizeof(int), 1, f);
    init_quadrule(qp, n,1);
    fread(qp->t[0], sizeof(double), n, f);
    fread(qp->wt, sizeof(double), n, f);
    fclose(f);
}

void import_quadrule(char* fname, QuadRule* qp) {
    int i,n, d;
    FILE* f = fopen(fname, "rb");
    fread(&n, sizeof(int), 1, f);
    fread(&d, sizeof(int), 1, f);
    init_quadrule(qp, n,d);
    for (i=0; i<d; i++) {
        fread(qp->t[i], sizeof(double), n, f);
    }
    fread(qp->wt, sizeof(double), n, f);
    fclose(f);
}

QuadRule resize_quadrule(QuadRule quadpoints, int size1, int size2){
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
       if(!quadpoints.t[i]){
          printf("Out of memory in resizequadrule\n");
          exit(1);
       }
   }
   return quadpoints;
}

QuadRule resize_quadrule_n(QuadRule quadpoints, int size){
   int i;
   quadpoints.n=size*quadpoints.n;
   for(i=0; i<quadpoints.d; i++){
      quadpoints.t[i]=(double*)realloc(quadpoints.t[i], quadpoints.n*sizeof(double));
   }
   quadpoints.wt=realloc(quadpoints.wt, quadpoints.n*sizeof(double));
   return quadpoints;
}

QuadRule resize_quadrule_d(QuadRule quadpoints, int size){
   int i;
   int oldd=quadpoints.d;
   quadpoints.d=size+quadpoints.d;
   quadpoints.t=(double**)realloc(quadpoints.t, quadpoints.d*sizeof(double*));
   for(i=oldd; i<quadpoints.d; i++){
       quadpoints.t[i]=(double*)malloc(quadpoints.n*sizeof(double));
   }
   return quadpoints;
}

void init_quadrule(QuadRule* qp, int n, int d) {
   int i;
   qp->n = n;
   qp->d = d;
   qp->t = (double**)malloc(qp->d*sizeof(double*));
   if(!qp->t){
      printf("Out of Memory in init_quadrule, d too large\n");
      exit(1);
   }
   for(i=0; i<qp->d; i++){
      qp->t[i] = (double*)malloc(qp->n*sizeof(double));
      if(!qp->t[i]){
         printf("Out of memory in init_quadrule, n too large\n");
         exit(1);
      }
   }
   qp->wt = (double*)malloc(qp->n*sizeof(double));
   if(!qp->wt){
       printf("Out of memory in init_quadrule, n too large\n");
       exit(1);
   }
}

void print_quadrule(QuadRule qr) {
   int i,j;
   printf("\n");
   for (i=0; i<qr.n; i++) {
      for (j=0; j<qr.d; j++) {
         printf("%lf ", qr.t[j][i]);
      }
      printf("\n");
   }
   printf("\n");

   for (i=0; i<qr.n; i++) {
      printf("%lf\n", qr.wt[i]);
   }
}

void free_quadrule(QuadRule qp){
   int i;
   for(i=0; i<qp.d; i++){
      free(qp.t[i]);
   }
   free(qp.t);
   free(qp.wt);
}


