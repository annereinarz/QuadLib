#include "matrix.h"
#include <mex.h>
#include <memory.h>

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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    int i;
    int  nr, ns;
    Vertexlist vtx;
    double* vtxh;
    QuadraturePoints SP;
    
   /*Get Input data */
   if(nrhs != 3){
      mexErrMsgTxt("Requires three input arguments");
   }
    
   init_vtxlist(&vtx, mxGetM(prhs[2]), mxGetN(prhs[2]));
      
   vtxh=mxGetPr(prhs[2]);
   for(i=0; i<vtx.s1*vtx.s2; i++){
      vtx.vtxlist[i%vtx.s1][i/vtx.s1]=vtxh[i];
   }
   nr=(int)(mxGetScalar(prhs[0]));
   ns=(int)(mxGetScalar(prhs[1]));
   
   SP=squad2dquad(vtx, nr, ns);
   
  /* Create Output Arguments */
  if(nlhs != 2){
    mexErrMsgTxt("Requires two output argument\n");
  }
  plhs[0] = (mxArray *) mxCreateDoubleMatrix(SP.n,SP.d,mxREAL);     /*vector of the quadrature points*/
  plhs[1] = (mxArray *) mxCreateDoubleMatrix(SP.n,1,mxREAL);     /*vector of the quadrature weights*/
  mxAssert(plhs[0] != NULL, "Out of memory");
  mxAssert(plhs[1] != NULL, "Out of memory");
  
  /* Copy the output */
  {
  int i,j;
  double* r = mxGetPr(plhs[0]);
  for (j=0; j<SP.d; j++) {
      memcpy(&r[j*SP.n], SP.t[j], SP.n*sizeof(double));
  }
  }
  memcpy(mxGetPr(plhs[1]), SP.wt, SP.n*sizeof(double));
  
  free_qp(SP);
  return;
}
