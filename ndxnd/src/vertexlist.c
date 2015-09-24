#include "vertexlist.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

Vertexlist import_vertexlist(char* fname){
    int i,j;
    Vertexlist vtxlist;
    FILE* f = fopen(fname, "r");
    if(! f) {
       fprintf(stderr, "Error: Unable to read file \"%s\".\n", fname);
       exit(1);
    }
    fscanf(f, "%d", &vtxlist.s1);
    fscanf(f, "%d", &vtxlist.s2);
    vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
    for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
    }
    for (i=0; i<vtxlist.s1; i++) {
       for(j=0; j<vtxlist.s2; j++) {
           fscanf(f, "%lf", &vtxlist.vtxlist[i][j]);
       }
    }
    fclose(f);
    return vtxlist;
}

void export_vertexlist(Vertexlist vtxlist, char *fname){
    int i,j;
    FILE* f = fopen(fname, "w");
    fprintf(f,"%d ", vtxlist.s1);
    fprintf(f,"%d\n",vtxlist.s2);
    for (i=0; i<vtxlist.s1; i++) {
       for(j=0; j<vtxlist.s2; j++) {
        fprintf(f, "%lf ",vtxlist.vtxlist[i][j]);
       }
       fprintf(f, "\n");
    }
    fclose(f);
}

void print_vertexlist(Vertexlist vl){
    for (int i=0; i<vl.s1; i++) {
      for (int j=0; j<vl.s2; j++) {
        printf("%lf ", vl.vtxlist[i][j]);
      }
      printf("\n");
    }
}

void free_vertexlist(Vertexlist vtxlist){
   int i;
   for(i=0; i<vtxlist.s1; i++){
       free(vtxlist.vtxlist[i]);
   }
   free(vtxlist.vtxlist);
}

 /*
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
       vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
       vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
       vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
       break;
     case 2:
       vtxlist.s1=2;
       vtxlist.s2=3;
       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
       vtxlist.vtxlist[0][1]=0.5;       vtxlist.vtxlist[1][1]=0;
       vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=0.5;
       break;
     case 3:
       vtxlist.s1=2;
       vtxlist.s2=4;
       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0.5;       vtxlist.vtxlist[1][0]=0;
       vtxlist.vtxlist[0][1]=0.5;       vtxlist.vtxlist[1][1]=0.5;
       vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=0;
       vtxlist.vtxlist[0][3]=1;         vtxlist.vtxlist[1][3]=0;
       break;
     case 4:
       vtxlist.s1=2;
       vtxlist.s2=5;
       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0.5;       vtxlist.vtxlist[1][0]=0.5;
       vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0.5;
       vtxlist.vtxlist[0][2]=0.5;       vtxlist.vtxlist[1][2]=1;
       vtxlist.vtxlist[0][3]=0.5;       vtxlist.vtxlist[1][3]=0;
       vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=0.5;
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
         vtxlist.vtxlist[0][3]=-0.5;        vtxlist.vtxlist[1][3]=-0.5;
         vtxlist.vtxlist[0][4]=-1.5;         vtxlist.vtxlist[1][4]=-0.5;
         vtxlist.vtxlist[0][5]=-0.5;        vtxlist.vtxlist[1][5]=-1.5;
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
         vtxlist.vtxlist[0][1]=0;         vtxlist.vtxlist[1][1]=1;
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         vtxlist.vtxlist[0][3]=-1;         vtxlist.vtxlist[1][3]=0;
         break;

      case 22:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=0;         vtxlist.vtxlist[1][1]=1;
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         break;

     case 29:
         vtxlist.s1=3;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=0;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=1;
         vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=0;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
         vtxlist.vtxlist[0][4]=-1;   vtxlist.vtxlist[1][4]=-1;     vtxlist.vtxlist[2][4]=-1;
         vtxlist.vtxlist[0][5]=-2;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=-1;
         vtxlist.vtxlist[0][6]=-1;   vtxlist.vtxlist[1][6]=-2;      vtxlist.vtxlist[2][6]=-1;
         vtxlist.vtxlist[0][7]=-1;   vtxlist.vtxlist[1][7]=-1;     vtxlist.vtxlist[2][7]=-2;
         break;

       case 30:
         vtxlist.s1=3;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }
         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=0;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=1;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=0;
         vtxlist.vtxlist[0][4]=0;   vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=-1;
         vtxlist.vtxlist[0][5]=-1;    vtxlist.vtxlist[1][5]=0;     vtxlist.vtxlist[2][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;      vtxlist.vtxlist[2][6]=0;
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
         vtxlist.vtxlist[0][5]=-2;   vtxlist.vtxlist[1][5]=-2;     vtxlist.vtxlist[2][5]=-2;    vtxlist.vtxlist[3][5]=-2;
         vtxlist.vtxlist[0][6]=-1;    vtxlist.vtxlist[1][6]=-2;     vtxlist.vtxlist[2][6]=-2;    vtxlist.vtxlist[3][6]=-2;
         vtxlist.vtxlist[0][7]=-2;   vtxlist.vtxlist[1][7]=-1;      vtxlist.vtxlist[2][7]=-2;     vtxlist.vtxlist[3][7]=-2;
         vtxlist.vtxlist[0][8]=-2;   vtxlist.vtxlist[1][8]=-2;     vtxlist.vtxlist[2][8]=-1;     vtxlist.vtxlist[3][8]=-2;
         vtxlist.vtxlist[0][9]=-2;   vtxlist.vtxlist[1][9]=-2;     vtxlist.vtxlist[2][9]=-2;     vtxlist.vtxlist[3][9]=-1;

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

      case 49:
         vtxlist.s1=5;
         vtxlist.s2=12;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         vtxlist.vtxlist[0][6]=-2;    vtxlist.vtxlist[1][6]=-2;      vtxlist.vtxlist[2][6]=-2;     vtxlist.vtxlist[3][6]=-2;    vtxlist.vtxlist[4][6]=-2;
         vtxlist.vtxlist[0][7]=-1;    vtxlist.vtxlist[1][7]=-2;      vtxlist.vtxlist[2][7]=-2;     vtxlist.vtxlist[3][7]=-2;    vtxlist.vtxlist[4][7]=-2;
         vtxlist.vtxlist[0][8]=-2;    vtxlist.vtxlist[1][8]=-1;      vtxlist.vtxlist[2][8]=-2;     vtxlist.vtxlist[3][8]=-2;    vtxlist.vtxlist[4][8]=-2;
         vtxlist.vtxlist[0][9]=-2;    vtxlist.vtxlist[1][9]=-2;      vtxlist.vtxlist[2][9]=-1;     vtxlist.vtxlist[3][9]=-2;    vtxlist.vtxlist[4][9]=-2;
         vtxlist.vtxlist[0][10]=-2;    vtxlist.vtxlist[1][10]=-2;      vtxlist.vtxlist[2][10]=-2;     vtxlist.vtxlist[3][10]=-1;    vtxlist.vtxlist[4][10]=-2;
         vtxlist.vtxlist[0][11]=-2;    vtxlist.vtxlist[1][11]=-2;      vtxlist.vtxlist[2][11]=-2;     vtxlist.vtxlist[3][11]=-2;    vtxlist.vtxlist[4][11]=-1;
         break;

      case 50:
         vtxlist.s1=5;
         vtxlist.s2=11;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         vtxlist.vtxlist[0][6]=-1;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=0;    vtxlist.vtxlist[4][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=-1;      vtxlist.vtxlist[2][7]=0;     vtxlist.vtxlist[3][7]=0;    vtxlist.vtxlist[4][7]=0;
         vtxlist.vtxlist[0][8]=0;    vtxlist.vtxlist[1][8]=0;      vtxlist.vtxlist[2][8]=-1;     vtxlist.vtxlist[3][8]=0;    vtxlist.vtxlist[4][8]=0;
         vtxlist.vtxlist[0][9]=0;    vtxlist.vtxlist[1][9]=0;      vtxlist.vtxlist[2][9]=0;     vtxlist.vtxlist[3][9]=-1;    vtxlist.vtxlist[4][9]=0;
         vtxlist.vtxlist[0][10]=0;    vtxlist.vtxlist[1][10]=0;      vtxlist.vtxlist[2][10]=0;     vtxlist.vtxlist[3][10]=0;    vtxlist.vtxlist[4][10]=-1;
         break;

      case 51:
         vtxlist.s1=5;
         vtxlist.s2=10;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;      vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=0;    vtxlist.vtxlist[4][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=-1;     vtxlist.vtxlist[3][7]=0;    vtxlist.vtxlist[4][7]=0;
         vtxlist.vtxlist[0][8]=0;    vtxlist.vtxlist[1][8]=0;      vtxlist.vtxlist[2][8]=0;     vtxlist.vtxlist[3][8]=-1;    vtxlist.vtxlist[4][8]=0;
         vtxlist.vtxlist[0][9]=0;    vtxlist.vtxlist[1][9]=0;      vtxlist.vtxlist[2][9]=0;     vtxlist.vtxlist[3][9]=0;    vtxlist.vtxlist[4][9]=-1;
         break;

      case 52:
         vtxlist.s1=5;
         vtxlist.s2=9;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;     vtxlist.vtxlist[3][6]=0;    vtxlist.vtxlist[4][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=0;     vtxlist.vtxlist[3][7]=-1;    vtxlist.vtxlist[4][7]=0;
         vtxlist.vtxlist[0][8]=0;    vtxlist.vtxlist[1][8]=0;      vtxlist.vtxlist[2][8]=0;     vtxlist.vtxlist[3][8]=0;    vtxlist.vtxlist[4][8]=-1;
         break;

      case 53:
         vtxlist.s1=5;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=-1;    vtxlist.vtxlist[4][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=0;     vtxlist.vtxlist[3][7]=0;    vtxlist.vtxlist[4][7]=-1;
         break;

      case 54:
         vtxlist.s1=5;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=0;    vtxlist.vtxlist[4][6]=-1;
         break;

      case 55:
         vtxlist.s1=5;
         vtxlist.s2=6;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;    vtxlist.vtxlist[4][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;    vtxlist.vtxlist[4][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;    vtxlist.vtxlist[4][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;    vtxlist.vtxlist[4][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;    vtxlist.vtxlist[4][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;    vtxlist.vtxlist[4][5]=1;
         break;
     default: printf("Invalid vertexlist chosen\n");
              exit(1);

   }
 return vtxlist;

}
   */
