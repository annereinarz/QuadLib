#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "sing_commonvertex2dtria.h"
#include "../mxMemory.h"

void sing_commonvertex2dtria(QuadraturePoints* SP){

     int i;
     double th;
     double *pth;
     
     pth=SP->t[0];
     SP->t[0]=SP->t[3];
     SP->t[3]=pth;

     for(i=0; i<SP->n; i++){
        SP->wt[i]=SP->wt[i]*SP->t[3][i];
        SP->t[2][i]=SP->t[2][i]*SP->t[3][i];
     }

     resizequadpointsn(SP, 2);

     for(i=0; i<SP->n/2; i++){
        SP->wt[SP->n/2+i]=SP->wt[i];
        SP->t[2][SP->n/2+i]=SP->t[3][i];
        SP->t[3][SP->n/2+i]=SP->t[2][i];
        SP->t[0][SP->n/2+i]=SP->t[0][i];
        SP->t[1][SP->n/2+i]=SP->t[1][i];
     }
     
     for(i=0; i<SP->n; i++){
        th=SP->t[0][i];
        SP->t[0][i]=SP->t[2][i];
        SP->t[2][i]=th;
     }

     for(i=0; i<SP->n; i++){
        SP->wt[i]=SP->wt[i]*SP->t[0][i]*SP->t[3][i];
     }

     for(i=0; i<SP->n; i++){
        SP->t[1][i]=SP->t[1][i]*SP->t[0][i];
        SP->t[2][i]=SP->t[2][i]*SP->t[3][i];
     }

     for(i=0; i<SP->n; i++){
        SP->t[0][i]=SP->t[0][i]-SP->t[1][i];
        SP->t[3][i]=SP->t[3][i]-SP->t[2][i];
     }
     return;

}

