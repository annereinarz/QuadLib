#include "matrix.h"
#include <mex.h>
#include <memory.h>

#include "quadpoints.h"
#include "squad.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    int  nr, ns;
    int whichvtxlist;
    int whichF=1;
    double Q;
    /* QuadraturePoints SP;*/

   /*Get Input data */
   if(nrhs != 4){
      mexErrMsgTxt("Requires four input arguments");
   }
   nr=(int)(mxGetScalar(prhs[0]));
   ns=(int)(mxGetScalar(prhs[1]));
   whichvtxlist=(int)(mxGetScalar(prhs[2]));
   whichF=(int)(mxGetScalar(prhs[3]));

   /*Get Quadrature rule */
   Q=squad(nr, ns, whichvtxlist, whichF);

  /* Create Output Arguments */
  if(nlhs != 1){
    mexErrMsgTxt("Requires one output argument\n");
  }
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
  *mxGetPr(plhs[0]) = Q;

  /*if(nlhs != 2){
    mexErrMsgTxt("Requires two output argument\n");
  }
  plhs[0] = (mxArray *) mxCreateDoubleMatrix(SP.n,SP.d,mxREAL);     
  plhs[1] = (mxArray *) mxCreateDoubleMatrix(SP.n,1,mxREAL);    
  mxAssert(plhs[0] != NULL, "Out of memory");
  mxAssert(plhs[1] != NULL, "Out of memory");   
  {
  int i,j;
  double* r = mxGetPr(plhs[0]);
  for (j=0; j<SP.d; j++) {
      memcpy(&r[j*SP.n], SP.t[j], SP.n*sizeof(double));
  }
  }
  memcpy(mxGetPr(plhs[1]), SP.wt, SP.n*sizeof(double));
        
free_qp(SP);*/
  return;
}
