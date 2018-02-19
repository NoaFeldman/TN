
/* compilation:
   mex my_inverse.cc -largeArrayDims -lmwblas -lmwlapack
   a=my_inverse(rand(3)); => crashes MatLab

   if I comment out the dgetri(...) line, MatLab survives.
   
   compilation with -g option and starting in debug mode writes:
   Program received signal SIGSEGV, Segmentation fault.
#0 0x00007fffc314b821 in LSteps1gas_1 () from ...
   /amnt/dist64/DIR/matlab-R2011a/bin/glnxa64/../../bin/glnxa64/../../bin/glnxa64/mkl.so

   is this an MKL library !?

   Wb,Oct04,12
*/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <string.h>

#include <mex.h>
#include <mat.h>
#include <lapack.h>

#define pINT ptrdiff_t

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   if (nargin!=1 || !mxIsDouble(argin[0])) {
   printf("\n  ERR invalid usage\n"); exit(0); }

   const mwSize *S=mxGetDimensions(argin[0]);
   const double *M0=mxGetPr(argin[0]);
   pINT lda=S[0], n=S[1], lwork=lda*n, i=0, ipiv[n];
   size_t s=lda*n*sizeof(double);

   double M[lda*n], work[lwork];
   memcpy(M,M0,s);

   printf("\n   TST: %dx%d @ 0x%lX [0x%lX]\n",lda,n,M,M0);
   if (i) {
      if (i<0)
           printf("\n  ERR dgetri() got invalid argument #%d\n",-i);
      else printf("\n  ERR dgetri() got singular matrix (i=%d)\n",i);
   }

   if (nargout) {
      argout[0]=mxCreateDoubleMatrix(lda,n,mxREAL);
      memcpy(mxGetPr(argout[0]),M,s);
   }
};

