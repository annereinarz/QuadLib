#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "quadpoints.h"
#include "sing_identical1d.h"
#include "sing_commonvertex1d.h"
#include "Ref2PhyS1d_affine.h"
#include "Ref2PhyS1d_circle.h"
#include "squad1d.h"

QuadraturePoints squad1d(Vertexlist vtx, int nr, int ns, int flag){
    QuadraturePoints SP;
    int k;
    
    /* Get Quadrature Rule*/
    QuadraturePoints qp1, qp2;
    if(vtx.s2==4){   /*distant elements*/
       init_quadpoints(&qp1 ,nr ,1);
       GLquad(&qp1,nr, 0, 1);       /* regular*/
    }
    else{ 
       init_quadpoints(&qp1, ns*(ns+1)/2,1); 
       CGLquad(&qp1, ns);         /* singular*/
    }
    init_quadpoints(&qp2 ,nr ,1);
    GLquad(&qp2,nr, 0, 1);       /* regular*/
    init_quadpoints(&SP, qp1.n*qp2.n, 2);
    TensorQuad(&SP, qp2, qp1);
    free_qp(qp1);  free_qp(qp2);
    
    switch(vtx.s2){
            case 2 : sing_identical1d(&SP);
                     k=1;
            break;
            
            case 3 : sing_commonvertex1d(&SP);
                     k=0;
            break;
            
            case 4 :
                     k=-1;
            break;
            
        default : printf("Wrong number of vertices! Expected: 2, 3 or 4. Got: %d\n", vtx.s2);
            exit(1);
    }

    
    if(flag==0){
      Ref2PhyS1d_affine(&SP, vtx, k);  
    }
    if(flag==1){
      Ref2PhyS1d_circle(&SP, vtx, k);  
    }
 
    return SP;
}

