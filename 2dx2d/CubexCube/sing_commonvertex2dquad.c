#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "sing_commonvertex2dquad.h"

void sing_commonvertex2dquad(QuadraturePoints* SP){

     int i;
          
     resizequadpointsn(SP, 4);
     for(i=0; i<SP->n/4; i++){
         SP->wt[i]=SP->wt[i]*pow(SP->t[0][i], 3);
         SP->wt[i+SP->n/4]=SP->wt[i];
         SP->wt[i+2*SP->n/4]=SP->wt[i];
         SP->wt[i+3*SP->n/4]=SP->wt[i];
         
         SP->t[0][i+SP->n/4]=SP->t[0][i]*SP->t[1][i];         
         SP->t[1][i+SP->n/4]=SP->t[0][i];
         SP->t[2][i+SP->n/4]=SP->t[0][i]*SP->t[2][i];
         SP->t[3][i+SP->n/4]=SP->t[0][i]*SP->t[3][i];
         
         SP->t[0][i+2*SP->n/4]=SP->t[0][i]*SP->t[1][i];
         SP->t[1][i+2*SP->n/4]=SP->t[0][i]*SP->t[2][i];         
         SP->t[2][i+2*SP->n/4]=SP->t[0][i];
         SP->t[3][i+2*SP->n/4]=SP->t[0][i]*SP->t[3][i];

         SP->t[0][i+3*SP->n/4]=SP->t[0][i]*SP->t[1][i];
         SP->t[1][i+3*SP->n/4]=SP->t[0][i]*SP->t[2][i];
         SP->t[2][i+3*SP->n/4]=SP->t[0][i]*SP->t[3][i];
         SP->t[3][i+3*SP->n/4]=SP->t[0][i];      
          
         SP->t[0][i]=SP->t[0][i];
         SP->t[1][i]=SP->t[0][i]*SP->t[1][i];
         SP->t[2][i]=SP->t[0][i]*SP->t[2][i];
         SP->t[3][i]=SP->t[0][i]*SP->t[3][i];
     }     
     
return;

}

