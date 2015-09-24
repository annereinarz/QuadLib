#include "matrix.h"
#include <mex.h>
#include "memory.h"

#include "quadpoints.h"
#include "squad1d.h"

/*
Calculates a quadrature rule on intervals
INPUT:
    flag  : Can be:
 *             0   :  gauss legendre/composite gauss legendre rules in [0,1]^2 are transformed, 
 *                     number of points to use given by the second and third argument with affine transformation
 *             1   :  evaluation only on one point in [0,1]^2, point given by second and third argument
 *             2   :  block wise evaluation of the quadrature rules, NOT IMPLEMENTED YET!
 *             3   :  gauss legendre/composite gauss legendre with transformation on the circle   
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    int i;
    int  nr, ns;
    double x1, x2;
    int flag;
    Vertexlist vtx;
    QuadraturePoints SP;
    double* vtxh;
    
   /*Get Input data */
   if(nrhs != 4){
      mexErrMsgTxt("Requires four input arguments");
   }
   flag=(int)(mxGetScalar(prhs[0])); 
    
   init_vtxlist(&vtx, mxGetM(prhs[3]), mxGetN(prhs[3]));
 
   vtxh=mxGetPr(prhs[3]);
   for(i=0; i<vtx.s1*vtx.s2; i++){
      vtx.vtxlist[i%vtx.s1][i/vtx.s1]=vtxh[i];
   }
   
   if(flag==0){
       nr=(int)(mxGetScalar(prhs[1]));
       ns=(int)(mxGetScalar(prhs[2]));
        /*Get Quadrature rule */
        SP=squad1d(vtx, nr, ns,0);
   }
   if(flag==3){
       nr=(int)(mxGetScalar(prhs[1]));
       ns=(int)(mxGetScalar(prhs[2]));
        /*Get Quadrature rule */
        SP=squad1d(vtx, nr, ns,1);
   }

  

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
