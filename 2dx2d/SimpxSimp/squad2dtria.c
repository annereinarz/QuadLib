#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "squad2dtria.h"
#include "sing_identical2dtria.h"
#include "sing_commonvertex2dtria.h"
#include "sing_distant2dtria.h"
#include "sing_commonedge2dtria.h"
#include "Ref2PhyS2d_affine.h"
#include "Ref2PhyS2d_spherical.h"

QuadraturePoints squad2dtria(Vertexlist vtx, int nr, int ns){
   int k;
   QuadraturePoints SP;

   QuadRule_cgl_gl(&SP, nr,ns, vtx.s2);


   switch(vtx.s2){
         case 3 : sing_identical2dtria(&SP);
                  k=2;
         break;

         case 5 : sing_commonvertex2dtria(&SP);
                  k=0;
         break;

         case 4 : sing_commonedge2dtria(&SP);
                  k=1;
         break;

         case 6 : sing_distant2dtria(&SP);
                  k=-1;
         break;

         default : printf("Too many vertices given\n");
           exit(1);
     }     

   /* Map to Physical Simplex */
   if(vtx.s1==2){
       Ref2PhyS2d_affine(&SP, vtx, k);
   }
   if(vtx.s1==3){
       Ref2PhyS2d_spherical(&SP, vtx, k);
   }
   
   return SP;
}

