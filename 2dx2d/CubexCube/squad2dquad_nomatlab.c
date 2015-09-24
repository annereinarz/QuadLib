#include <stdio.h>

#include "../quadpoints.h"
#include "../nomatlab_utils.h"
#include "squad2dquad.h"


/* Calculates Quadrature Rules on simplices
 * INPUT:
 * nr, ns  :  number of quadrature points in regular and singular direction
 * vtx     :  list of vertices of the two cubes. The common vertices should be listed first,
              then the remaining vertices of the first cube and then the remaining vertices 
              of the second cube. Can have either 2 or 3 dimension. In the case of two dimension+
              an affine mapping to the simplices is used, in the case of three dimensions, the quadrature
              rule is mapped to its position on a sphere of radius one.
 */


int main(int argc, char **argv){

    int  nr, ns;
    Vertexlist vtx;
    QuadraturePoints SP;
    
   /*Get Input data */
   
   assert_correct_number_of_inputs (argc, 2);
   nr = get_integer_argument(argv[1], "first input argument must be integer!");
   ns = get_integer_argument(argv[2], "second input argument must be integer!");
    
   read_vertexlist(&vtx);
      
   SP=squad2dquad(vtx, nr, ns);
   
  print_quadpoints(SP);
  
  free_vtxlist(&vtx); 
  free_qp(SP);
  return 0;
}
