#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "../mxMemory.h"
#include "Ref2PhyS2d_spherical.h"

void Ref2PhyS2d_spherical(QuadraturePoints *SP, Vertexlist vtx, int k){

     int i;
     double th1, th2;
     double detF;
     double *Fu1, *Fu2, *Fu3;
     
     switch(vtx.s2){
         case 3 : resizequadpointsd(SP, 2);
                  Fu1=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu2=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu3=(double*)mxMalloc(SP->n*sizeof(double));
                  for(i=0; i<SP->n; i++){
                      th1=SP->t[2][i];
                      th2=SP->t[3][i];
                      
                      SP->t[2][i]=1-SP->t[0][i]-SP->t[1][i];
                      SP->t[3][i]=th1;
                      SP->t[4][i]=th2;
                      SP->t[5][i]=1-th1-th2;
                  }  
                  for(i=0; i<SP->n; i++){
                      SP->t[2][i]=1-SP->t[0][i]-SP->t[1][i];
                  }
                  detF=vtx.vtxlist[0][0]*(vtx.vtxlist[1][1]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][1])
                         -vtx.vtxlist[0][1]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][0])
                         +vtx.vtxlist[0][2]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][1]-vtx.vtxlist[1][1]*vtx.vtxlist[2][0]);
                  
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][0]*SP->t[0][i]+vtx.vtxlist[0][1]*SP->t[1][i]+vtx.vtxlist[0][2]*SP->t[2][i];
                      Fu2[i]=vtx.vtxlist[1][0]*SP->t[0][i]+vtx.vtxlist[1][1]*SP->t[1][i]+vtx.vtxlist[1][2]*SP->t[2][i];
                      Fu3[i]=vtx.vtxlist[2][0]*SP->t[0][i]+vtx.vtxlist[2][1]*SP->t[1][i]+vtx.vtxlist[2][2]*SP->t[2][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[0][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[1][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[2][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                      SP->t[3][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[4][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[5][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                  
         break;
         
         case 4 : resizequadpointsd(SP, 2);
                  Fu1=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu2=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu3=(double*)mxMalloc(SP->n*sizeof(double));
                  for(i=0; i<SP->n; i++){
                      th1=SP->t[2][i];
                      th2=SP->t[3][i];
                      
                      SP->t[2][i]=1-SP->t[0][i]-SP->t[1][i];
                      SP->t[3][i]=th1;
                      SP->t[4][i]=th2;
                      SP->t[5][i]=1-th1-th2;
                  }
                  detF=vtx.vtxlist[0][0]*(vtx.vtxlist[1][1]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][1])
                         -vtx.vtxlist[0][1]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][0])
                         +vtx.vtxlist[0][2]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][1]-vtx.vtxlist[1][1]*vtx.vtxlist[2][0]);
                  
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][0]*SP->t[0][i]+vtx.vtxlist[0][1]*SP->t[1][i]+vtx.vtxlist[0][2]*SP->t[2][i];
                      Fu2[i]=vtx.vtxlist[1][0]*SP->t[0][i]+vtx.vtxlist[1][1]*SP->t[1][i]+vtx.vtxlist[1][2]*SP->t[2][i];
                      Fu3[i]=vtx.vtxlist[2][0]*SP->t[0][i]+vtx.vtxlist[2][1]*SP->t[1][i]+vtx.vtxlist[2][2]*SP->t[2][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[0][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[1][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[2][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                  detF=vtx.vtxlist[0][0]*(vtx.vtxlist[1][1]*vtx.vtxlist[2][3]-vtx.vtxlist[1][3]*vtx.vtxlist[2][1])
                         -vtx.vtxlist[0][1]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][3]-vtx.vtxlist[1][3]*vtx.vtxlist[2][0])
                         +vtx.vtxlist[0][3]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][1]-vtx.vtxlist[1][1]*vtx.vtxlist[2][0]);
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][0]*SP->t[3][i]+vtx.vtxlist[0][1]*SP->t[4][i]+vtx.vtxlist[0][3]*SP->t[5][i];
                      Fu2[i]=vtx.vtxlist[1][0]*SP->t[3][i]+vtx.vtxlist[1][1]*SP->t[4][i]+vtx.vtxlist[1][3]*SP->t[5][i];
                      Fu3[i]=vtx.vtxlist[2][0]*SP->t[3][i]+vtx.vtxlist[2][1]*SP->t[4][i]+vtx.vtxlist[2][3]*SP->t[5][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[3][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[4][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[5][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                  
         break;
         
         case 5 : resizequadpointsd(SP, 2);
                  Fu1=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu2=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu3=(double*)mxMalloc(SP->n*sizeof(double));
                  for(i=0; i<SP->n; i++){
                      th1=SP->t[2][i];
                      th2=SP->t[3][i];
                      
                      SP->t[2][i]=1-SP->t[0][i]-SP->t[1][i];
                      SP->t[3][i]=th1;
                      SP->t[4][i]=th2;
                      SP->t[5][i]=1-th1-th2;
                  }
                  detF=vtx.vtxlist[0][0]*(vtx.vtxlist[1][1]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][1])
                         -vtx.vtxlist[0][1]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][0])
                         +vtx.vtxlist[0][2]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][1]-vtx.vtxlist[1][1]*vtx.vtxlist[2][0]);
                  
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][0]*SP->t[0][i]+vtx.vtxlist[0][1]*SP->t[1][i]+vtx.vtxlist[0][2]*SP->t[2][i];
                      Fu2[i]=vtx.vtxlist[1][0]*SP->t[0][i]+vtx.vtxlist[1][1]*SP->t[1][i]+vtx.vtxlist[1][2]*SP->t[2][i];
                      Fu3[i]=vtx.vtxlist[2][0]*SP->t[0][i]+vtx.vtxlist[2][1]*SP->t[1][i]+vtx.vtxlist[2][2]*SP->t[2][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[0][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[1][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[2][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                  detF=vtx.vtxlist[0][0]*(vtx.vtxlist[1][3]*vtx.vtxlist[2][4]-vtx.vtxlist[1][4]*vtx.vtxlist[2][3])
                         -vtx.vtxlist[0][3]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][4]-vtx.vtxlist[1][4]*vtx.vtxlist[2][0])
                         +vtx.vtxlist[0][4]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][3]-vtx.vtxlist[1][3]*vtx.vtxlist[2][0]);
                  
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][0]*SP->t[3][i]+vtx.vtxlist[0][3]*SP->t[4][i]+vtx.vtxlist[0][4]*SP->t[5][i];
                      Fu2[i]=vtx.vtxlist[1][0]*SP->t[3][i]+vtx.vtxlist[1][3]*SP->t[4][i]+vtx.vtxlist[1][4]*SP->t[5][i];
                      Fu3[i]=vtx.vtxlist[2][0]*SP->t[3][i]+vtx.vtxlist[2][3]*SP->t[4][i]+vtx.vtxlist[2][4]*SP->t[5][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[3][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[4][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[5][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                     

         break;
         
         case 6 : Fu1=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu2=(double*)mxMalloc(SP->n*sizeof(double));
                  Fu3=(double*)mxMalloc(SP->n*sizeof(double));
                  resizequadpointsd(SP, 2);
                  for(i=0; i<SP->n; i++){
                      th1=SP->t[2][i];
                      th2=SP->t[3][i];
                      
                      SP->t[2][i]=1-SP->t[0][i]-SP->t[1][i];
                      SP->t[3][i]=th1;
                      SP->t[4][i]=th2;
                      SP->t[5][i]=1-th1-th2;
                  }
                  
                  detF=vtx.vtxlist[0][0]*(vtx.vtxlist[1][1]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][1])
                         -vtx.vtxlist[0][1]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][2]-vtx.vtxlist[1][2]*vtx.vtxlist[2][0])
                         +vtx.vtxlist[0][2]*(vtx.vtxlist[1][0]*vtx.vtxlist[2][1]-vtx.vtxlist[1][1]*vtx.vtxlist[2][0]);
                  
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][0]*SP->t[0][i]+vtx.vtxlist[0][1]*SP->t[1][i]+vtx.vtxlist[0][2]*SP->t[2][i];
                      Fu2[i]=vtx.vtxlist[1][0]*SP->t[0][i]+vtx.vtxlist[1][1]*SP->t[1][i]+vtx.vtxlist[1][2]*SP->t[2][i];
                      Fu3[i]=vtx.vtxlist[2][0]*SP->t[0][i]+vtx.vtxlist[2][1]*SP->t[1][i]+vtx.vtxlist[2][2]*SP->t[2][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[0][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[1][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[2][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                  
                  
                  detF=vtx.vtxlist[0][3]*(vtx.vtxlist[1][4]*vtx.vtxlist[2][5]-vtx.vtxlist[1][5]*vtx.vtxlist[2][4])
                         -vtx.vtxlist[0][4]*(vtx.vtxlist[1][3]*vtx.vtxlist[2][5]-vtx.vtxlist[1][5]*vtx.vtxlist[2][3])
                         +vtx.vtxlist[0][5]*(vtx.vtxlist[1][3]*vtx.vtxlist[2][4]-vtx.vtxlist[1][4]*vtx.vtxlist[2][3]);
                  
                  for(i=0; i<SP->n; i++){
                      Fu1[i]=vtx.vtxlist[0][3]*SP->t[3][i]+vtx.vtxlist[0][4]*SP->t[4][i]+vtx.vtxlist[0][5]*SP->t[5][i];
                      Fu2[i]=vtx.vtxlist[1][3]*SP->t[3][i]+vtx.vtxlist[1][4]*SP->t[4][i]+vtx.vtxlist[1][5]*SP->t[5][i];
                      Fu3[i]=vtx.vtxlist[2][3]*SP->t[3][i]+vtx.vtxlist[2][4]*SP->t[4][i]+vtx.vtxlist[2][5]*SP->t[5][i];
                  }
                  for(i=0; i<SP->n; i++){
                      SP->t[3][i]=Fu1[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[4][i]=Fu2[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->t[5][i]=Fu3[i]/(sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2)));
                      SP->wt[i]=SP->wt[i]*(detF/pow((sqrt(pow(Fu1[i],2)+pow(Fu2[i],2)+pow(Fu3[i],2))),3));
                  }
                  
         break;
         
         default : printf("Too many vertices given\n");
         exit(1);
     }
    return;
}
