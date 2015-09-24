#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "squad2dquad.h"
#include "sing_identical2dquad.h"
#include "sing_commonvertex2dquad.h"
#include "sing_commonedge2dquad.h"
#include "Ref2PhyJ2d.h"
#include "Ref2PhyJ2d.h"

QuadraturePoints squad2dquad(Vertexlist vtx, int nr, int ns){
   int k;
   QuadraturePoints SP;

   QuadRule_cgl_gl(&SP, nr, ns, vtx.s2);

   switch(vtx.s2){
         case 3 : sing_identical2dquad(&SP);
                  k=2;
         break;

         case 4 : sing_commonedge2dquad(&SP);
                  k=1;
         break;

         case 5 : sing_commonvertex2dquad(&SP);
                  k=0;
         break;

         case 6 :k=-1;
         break;

         default : printf("Too many vertices given\n");
           exit(1);
     }

   /* Map to Physical Simplex */
   Ref2PhyJ2d(&SP, vtx, k);
   return SP;
}
