#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "../mxMemory.h"
#include "sing_distant2dtria.h"

void sing_distant2dtria(QuadraturePoints* SP){
   int i;
   for(i=0; i<SP->n; i++){
       SP->t[1][i]=SP->t[1][i]*(1-SP->t[0][i]);
       SP->t[3][i]=SP->t[3][i]*(1-SP->t[2][i]);
       SP->wt[i]=SP->wt[i]*(1-SP->t[0][i])*(1-SP->t[2][i]);
   }
   
   return;
}
