#include "AffineTrafo.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

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

void freeAffineTrafo(AffineTrafo AfTr){
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

double determinant(double**a, int n){
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1){
       fprintf(stderr,"Error: Can't calculate determinant. Invalid spatial dimension!\n");
       exit(1);
   }
   else if (n == 1){
      det = a[0][0];
   }
   else if (n == 2){
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   }
   else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return det;
}
