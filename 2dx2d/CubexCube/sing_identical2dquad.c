#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "sing_identical2dquad.h"


void sing_identical2dquad(QuadraturePoints* SP){

     int i;

     /*Step 5: Transform to pyramids */
     resizequadpointsn(SP, 2);
     for(i=0; i<SP->n/2; i++){
         SP->wt[i]        =SP->t[0][i]*SP->wt[i];
         SP->wt[i+SP->n/2]=SP->wt[i];

         SP->t[0][i+SP->n/2]=     SP->t[0][i];
         SP->t[1][i+SP->n/2]=     SP->t[0][i]*SP->t[1][i];
         
         SP->t[1][i]=SP->t[0][i+SP->n/2];
         SP->t[0][i]=SP->t[1][i+SP->n/2];
         
         SP->t[2][i+SP->n/2]=SP->t[2][i];
         SP->t[3][i+SP->n/2]=SP->t[3][i];
     }
     /*Step 4: reflections*/
     resizequadpointsn(SP, 4);
     
     for(i=0; i<SP->n/4; i++){
         SP->wt[i]=SP->wt[i]*(1-fabs(SP->t[0][i]))*(1-fabs(SP->t[1][i]));
         SP->wt[i+SP->n/4]=SP->wt[i];
         SP->wt[i+2*SP->n/4]=SP->wt[i];
         SP->wt[i+3*SP->n/4]=SP->wt[i];
     }
     
     for(i=0; i<SP->n/4; i++){
         SP->t[0][i+SP->n/4]= SP->t[0][i];
         SP->t[1][i+SP->n/4]=-SP->t[1][i];
         
         SP->t[0][i+2*SP->n/4]=-SP->t[0][i];
         SP->t[1][i+2*SP->n/4]= SP->t[1][i];
         
         SP->t[0][i+3*SP->n/4]= -SP->t[0][i];
         SP->t[1][i+3*SP->n/4]= -SP->t[1][i];
         
         SP->t[0][i]= SP->t[0][i];
         SP->t[1][i]= SP->t[1][i];
     }

     /*Step 3: Transform utilde  in J^2 to uhat in F_k*/     
     for(i=0; i<SP->n/4; i++){
         SP->t[2][i+SP->n/4]=                     (1-fabs(SP->t[0][i+SP->n/4]))*SP->t[2][i];
         SP->t[3][i+SP->n/4]=-SP->t[1][i+SP->n/4]+(1-fabs(SP->t[1][i+SP->n/4]))*SP->t[3][i];
         
         SP->t[2][i+2*SP->n/4]=-SP->t[0][i+2*SP->n/4]+(1-fabs(SP->t[0][i+2*SP->n/4]))*SP->t[2][i];
         SP->t[3][i+2*SP->n/4]=                       (1-fabs(SP->t[1][i+2*SP->n/4]))*SP->t[3][i];
         
         SP->t[2][i+3*SP->n/4]=-SP->t[0][i+3*SP->n/4]+(1-fabs(SP->t[0][i+3*SP->n/4]))*SP->t[2][i];
         SP->t[3][i+3*SP->n/4]=-SP->t[1][i+3*SP->n/4]+(1-fabs(SP->t[1][i+3*SP->n/4]))*SP->t[3][i];
         
         SP->t[2][i]=(1-fabs(SP->t[0][i]))*SP->t[2][i];
         SP->t[3][i]=(1-fabs(SP->t[1][i]))*SP->t[3][i];
     }
     
     /*Step 2: (u,v) to (z,u)*/
     resizequadpointsd(SP, 2);
     for(i=0; i<SP->n; i++){
         SP->t[4][i]=SP->t[0][i];    SP->t[5][i]=SP->t[1][i];
         SP->t[0][i]=SP->t[2][i]+SP->t[0][i];
         SP->t[1][i]=SP->t[3][i]+SP->t[1][i];
     }
     
     return;

}

