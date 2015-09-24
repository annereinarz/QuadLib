#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "../mxMemory.h"
#include "sing_commonedge2dquad.h"

void sing_commonedge2dquad(QuadraturePoints* SP){

   int i;
   
   /*Step 5: Transform to pyramids,
              split [0,1]x[0,1]x[0,1] into 3 pyramids*/
   resizequadpointsn(SP, 3);
   for(i=0; i<SP->n/3; i++){
       SP->wt[i]=SP->wt[i]*SP->t[0][i]*SP->t[0][i];
       SP->wt[i+SP->n/3]=SP->wt[i];     SP->wt[i+2*SP->n/3]=SP->wt[i];
       
       SP->t[3][i+SP->n/3]=  SP->t[3][i];  /*hat u variable is untouched*/
       SP->t[3][i+2*SP->n/3]=SP->t[3][i];
       
       SP->t[0][i+SP->n/3]=SP->t[0][i]*SP->t[1][i];
       SP->t[1][i+SP->n/3]=SP->t[0][i];
       SP->t[2][i+SP->n/3]=SP->t[0][i]*SP->t[2][i];
       
       SP->t[0][i+2*SP->n/3]=SP->t[0][i]*SP->t[1][i];
       SP->t[1][i+2*SP->n/3]=SP->t[0][i]*SP->t[2][i];
       SP->t[2][i+2*SP->n/3]=SP->t[0][i];
       
       SP->t[1][i]=SP->t[0][i]*SP->t[1][i];
       SP->t[2][i]=SP->t[0][i]*SP->t[2][i];
   }   

   /*Step 4 and 3: Transform utilde in [0,1] to uhat in F_1 and Reflections*/
   resizequadpointsn(SP, 2);
   for(i=0; i<SP->n/2; i++){
       SP->t[0][i+SP->n/2]=-SP->t[0][i];
   }
   
   for(i=0; i<SP->n/2; i++){
       SP->wt[i]=(1-SP->t[0][i])*SP->wt[i];
       SP->wt[i+SP->n/2]=SP->wt[i];
       
       SP->t[1][i+SP->n/2]= SP->t[1][i];
       SP->t[2][i+SP->n/2]= SP->t[2][i];
       
       SP->t[3][i+SP->n/2]=-SP->t[0][i+SP->n/2]+(1-fabs(SP->t[0][i+SP->n/2]))*SP->t[3][i];
       SP->t[3][i]        =                     (1-fabs(SP->t[0][i]))*        SP->t[3][i];
   }
          
   /*Step 2: (zhat, uhat, ucheck, vcheck) to (u,v,zhat)*/
   double th1,th2;
   resizequadpointsd(SP, 1);
   for(i=0; i<SP->n; i++){
      SP->t[4][i]=SP->t[0][i];   /*zhat*/
      SP->t[0][i]=SP->t[3][i]+SP->t[4][i];  /*vhat=zhat+uhat*/
   }
   
   for(i=0; i<SP->n; i++){
      th1=SP->t[2][i];
      SP->t[2][i]=SP->t[0][i];
      th2=SP->t[3][i];
      SP->t[3][i]=th1;
      SP->t[0][i]=th2;
   }

return;
}
