#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "../mxMemory.h"
#include "sing_commonedge2dtria.h"

void sing_commonedge2dtria(QuadraturePoints* SP){

   int i;
   double th;
   double *pth;
   
   pth=SP->t[0];
   SP->t[0]=SP->t[2];
   SP->t[2]=pth;
   
   for(i=0; i<SP->n; i++){
       SP->wt[i]=SP->t[2][i]*SP->t[2][i]*SP->wt[i];
   }
   for(i=0; i<SP->n; i++){
       SP->t[0][i]=SP->t[0][i]*SP->t[2][i];
       SP->t[1][i]=SP->t[1][i]*SP->t[2][i];
   }

   resizequadpointsn(SP, 3);            /* triple length */

   for(i=0; i<SP->n/3;  i++){
       SP->t[0][SP->n/3+i]=SP->t[2][i];
       SP->t[1][SP->n/3+i]=SP->t[0][i];
       SP->t[2][SP->n/3+i]=SP->t[1][i];
       
       SP->t[0][2*SP->n/3+i]=SP->t[1][i];
       SP->t[1][2*SP->n/3+i]=SP->t[2][i];
       SP->t[2][2*SP->n/3+i]=SP->t[0][i];
       
       SP->wt[SP->n/3+i]=SP->wt[i];
       SP->wt[2*SP->n/3+i]=SP->wt[i];
       
       SP->t[3][SP->n/3+i]=SP->t[3][i];
       SP->t[3][2*SP->n/3+i]=SP->t[3][i];
   }

   for(i=0; i<SP->n; i++){
       SP->wt[i]=SP->wt[i]*(1-SP->t[2][i]);
       SP->t[3][i]=SP->t[3][i]*(1-SP->t[2][i]);
   }

   resizequadpointsn(SP,2);                  /* double length */

   for(i=0; i<SP->n/2; i++){
      SP->wt[SP->n/2+i]=SP->wt[i];
      SP->t[0][SP->n/2+i]=  SP->t[0][i];
      SP->t[1][SP->n/2+i]=  SP->t[1][i];
      SP->t[2][SP->n/2+i]= -SP->t[2][i];
      SP->t[3][SP->n/2+i]=1-SP->t[3][i];
   }

   resizequadpointsd(SP,1);      /* add 1 row */

   for(i=0; i<SP->n; i++){
       SP->t[4][i]=SP->t[2][i];
       SP->t[2][i]=SP->t[2][i]+SP->t[3][i];
       th=SP->t[0][i];
       SP->t[0][i]=SP->t[3][i];
       SP->t[3][i]=SP->t[1][i];
       SP->t[1][i]=th;
   }

   for(i=0; i<SP->n; i++){
      SP->t[4][i]=(1-SP->t[3][i])*SP->t[4][i]+SP->t[0][i]*(SP->t[1][i]-SP->t[3][i]);
      SP->wt[i]=SP->wt[i]*(1-SP->t[1][i])*(1-SP->t[3][i]);
      SP->t[0][i]=SP->t[0][i]*(1-SP->t[1][i]);
      SP->t[2][i]=SP->t[2][i]*(1-SP->t[3][i]);
   }
   return;
}
