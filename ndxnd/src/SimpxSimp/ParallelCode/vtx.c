#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vtx.h"

void import_vtx(Vertexlist *vtx, FILE* f){

    int i;

    if (f == NULL) {
       printf("Error: Can't open file.");
       exit(1);
    }
    fread(&vtx->s1, sizeof(int), 1, f);
    fread(&vtx->s2, sizeof(int), 1, f);
    vtx->vtxlist = (double**)malloc(vtx->s1*sizeof(double*));
    for(i=0; i<vtx->s1; i++){
       vtx->vtxlist[i] = (double*)malloc(vtx->s2*sizeof(double));
    }
    for (i=0; i<vtx->s1; i++) {
        fread(vtx->vtxlist[i], sizeof(double), vtx->s2, f);
    }

}

Vertexlist validate_results(int i, int sim){
   Vertexlist vtxlist;

   switch(i){
       case 1:   // identical elements
       vtxlist.s1=3;
       vtxlist.s2=4;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=1;
       break;

       case 2:   // elements with a common face
       vtxlist.s1=3;
       vtxlist.s2=5;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=1;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
       vtxlist.vtxlist[0][4]=1;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=1;
       break;

       case 3:   // elements with a common edge 1
       vtxlist.s1=3;
       vtxlist.s2=6;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=1;      vtxlist.vtxlist[2][1]=1;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=0;      vtxlist.vtxlist[2][2]=0;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
       vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=1;
       vtxlist.vtxlist[0][5]=1;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=1;
       break;

       case 4:   // elements with a common edge 2
       vtxlist.s1=3;
       vtxlist.s2=6;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=1;      vtxlist.vtxlist[2][1]=1;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=0;      vtxlist.vtxlist[2][2]=0;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
       vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=1;
       vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=1;      vtxlist.vtxlist[2][5]=1;
       break;
   }

   return vtxlist;
}


Vertexlist get_vertexlist(int sim){
 int i;
 Vertexlist vtxlist;

 switch(sim){
    case 1:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=-1;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=sqrt(3);

         break;
    case 2:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0.5;       vtxlist.vtxlist[1][2]=0.5*sqrt(3);

         break;

    case 4:
         vtxlist.s1=2;
         vtxlist.s2=4;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=0.5;       vtxlist.vtxlist[1][1]=0.5*sqrt(3);
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         vtxlist.vtxlist[0][3]=-0.5;      vtxlist.vtxlist[1][3]=0.5*sqrt(3);

         break;

    case 5:
         vtxlist.s1=2;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0.5;       vtxlist.vtxlist[1][0]=0.5*sqrt(3);
         vtxlist.vtxlist[0][1]=0;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         vtxlist.vtxlist[0][3]=-0.5;      vtxlist.vtxlist[1][3]=0.5*sqrt(3);
         vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=sqrt(3);

         break;

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
         vtxlist.vtxlist[0][3]=-1;        vtxlist.vtxlist[1][3]=-1;
         vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=-1;
         vtxlist.vtxlist[0][5]=-1;        vtxlist.vtxlist[1][5]=0;

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
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         vtxlist.vtxlist[0][3]=0;         vtxlist.vtxlist[1][3]=-1;
         break;

      case 22:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         break;

      case 29:
         vtxlist.s1=3;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=-1;   vtxlist.vtxlist[1][4]=-1;     vtxlist.vtxlist[2][4]=-1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=-1;
         vtxlist.vtxlist[0][6]=-1;   vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;
         vtxlist.vtxlist[0][7]=-1;   vtxlist.vtxlist[1][7]=-1;     vtxlist.vtxlist[2][7]=0;
         break;

       case 30:
         vtxlist.s1=3;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=-1;   vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;
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

      case 39:
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
         break;

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
