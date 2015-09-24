#include <stdio.h>
#include <stdlib.h>

#include "quadpoints.h"
#include "squad1d.h"
#include "nomatlab_utils.h"


/*
Calculates a quadrature rule on intervals
INPUT:
    flag  : Can be:
 *             0   :  gauss legendre/composite gauss legendre rules in [0,1]^2 are transformed, 
 *                     number of points to use given by the second and third argument with affine transformation
 *             1      evaluation only on one point in [0,1]^2, point given by second and third argument
 *             2      block wise evaluation of the quadrature rules, NOT IMPLEMENTED YET!
 *             3      gauss legendre/composite gauss legendre with transformation on a circle
    n_gl  : the number of gauss legendre points
    n_cgl : number for the composite gauss legendre quadrature
               number of quadrature points is 0.5*n_cgl*(n_cgl+1)
    vtx    : list of vertices of the two intervals. The common vertices should be listed first,
             then the remaining vertices of the first interval and then the remaining vertices
             of the second interval. Can be one or two dimensional. In the case of two dimensions
             mapping onto the unit circle is used, otherwise an affine mapping to the physical intervals.

OUTPUT:
    [t,wt]  : quadrature rule on the two input intervals, with singularities as described in the sing_*1d routines.

*/  

int main(int argc, char **argv){

    int flag;
    Vertexlist vtx;
    QuadraturePoints SP;
    
   /*Get Input data */
   if(argc != 4){
      printf("Requires three input arguments (and the vertex list on stdin)!\n");
      return 1;
   }
   if (1 != sscanf(argv[1], "%d", &flag)) {
      printf("first input argument must be integer!\n");
      return 1;
   }
    
   read_vertexlist(&vtx);
   
   if(flag==0){
       int  nr, ns;
        if (1 != sscanf(argv[2], "%d", &nr)) {
            printf("second input argument must be integer!\n");
            return 1;
        }
        if (1 != sscanf(argv[3], "%d", &ns)) {
            printf("third input argument must be integer!\n");
            return 1;
        }
        /*Get Quadrature rule */
        SP=squad1d(vtx, nr, ns,0);
   }
   if(flag==3){
       int  nr, ns;
        if (1 != sscanf(argv[2], "%d", &nr)) {
            printf("second input argument must be integer!\n");
            return 1;
        }
        if (1 != sscanf(argv[3], "%d", &ns)) {
            printf("third input argument must be integer!\n");
            return 1;
        }
        /*Get Quadrature rule */
        SP=squad1d(vtx, nr, ns,1);
   }

  print_quadpoints(SP);
    
  free_vtxlist(&vtx); 
  free_qp(SP);
  return 0;
}
