#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "QuadRule.h"
#include "unique.h"

const double MARKER = -1;

int compare(const void* a, const void* b){
   double* tmpa=*(double**)a;
   double* tmpb=*(double**)b;

   int i=0;
   while(tmpa[i] != MARKER && tmpb[i] != MARKER) {
      if (fabs(tmpa[i] - tmpb[i]) > DBL_EPSILON) {
         if (tmpa[i] > tmpb[i]) {
            return 1;
         }
         return -1;
      }
      i++;
   }
   return 0;
}

/*int main(){

   QuadRule QP;
   import_quadrule("quadpoints",&QP);
   print_quadpoints(QP);

   unique(&QP);

   print_quadpoints(QP);
   free_quadrule(QP);
   return 0;
} */

void unique(QuadRule* QP){
      /* Uniq */
   int n=QP->n; int d=QP->d;
   int i,j;
   double ** tempQP;
   tempQP=(double**)malloc(n*sizeof(double*));
   for(i=0; i<n; i++){
      tempQP[i]=(double*)malloc((d+2)*sizeof(double));
   }

   for(i=0; i<n; i++){
      for(j=0; j<d;j++){
         tempQP[i][j]=QP->t[j][i];
      }
      tempQP[i][d]=MARKER;
      tempQP[i][d+1]=QP->wt[i];
   }

   free_quadrule(*QP);

   qsort(tempQP, n, sizeof(double*),compare);

   int* remove;
   remove=(int*)malloc(n*sizeof(int));
   memset(remove,0,n*sizeof(int));
   int cnt=0;
   for(i=0; i<n-1; i++){
      if(! compare(&tempQP[i],&tempQP[i+1])){
         remove[i]=1;
         tempQP[i+1][d+1]+=tempQP[i][d+1];
      }
      else{
         cnt++;
      }
   }
   cnt++;

   init_quadrule(QP, cnt, d);

   int k=0;
   for(i=0; i<n; i++){
      if(! remove[i]){
         for(j=0; j<d;j++){
            QP->t[j][k]=tempQP[i][j];
         }
         QP->wt[k]=tempQP[i][d+1];
         k++;
      }
   }

   for(i=0; i<n; i++){
      free(tempQP[i]);
   }
   free(tempQP);
   free(remove);

   /* End of Uniq */
}
