#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "quadpoints.h"
#include "sing_identical1d.h"

void sing_identical1d(QuadraturePoints* SP){

     int i;
     
     for(i=0; i<SP->n; i++){
        SP->wt[i]=SP->wt[i]*(1-SP->t[1][i]);
        SP->t[0][i]=SP->t[0][i]*(1-SP->t[1][i]);
        SP->t[1][i]=SP->t[1][i]+SP->t[0][i];
     }

     resizequadpointsn(SP, 2);
     
     for(i=0; i<SP->n/2; i++){
        SP->t[0][SP->n/2+i]=SP->t[1][i];
        SP->t[1][SP->n/2+i]=SP->t[0][i];
        SP->wt[SP->n/2+i]=SP->wt[i];
     }

    
}
