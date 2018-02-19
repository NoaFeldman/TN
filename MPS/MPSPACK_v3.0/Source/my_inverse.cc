
/* compilation: mex my_inverse.cc -lmwlapack
   NB! -lmwlapack is ESSENTIAL
   NB! also setup of IPIV[] using dgetrf(!!)
   otherwise MatLab happily crashes )!(*$)!(*P*I
   Wb,Oct05,12
*/

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <string>

#include <mex.h>
#include <mat.h>

#define pINT ptrdiff_t

#define dgetri dgetri_ 
#define dgetrf dgetrf_ 

extern "C" {
void dgetri(
   const pINT   &,
         double [],
   const pINT   &,
         pINT   [],
         double [],
   const pINT   &,
         pINT   &
);

void dgetrf(
   const pINT   &,
   const pINT   &,
         double [],
   const pINT   &,
         pINT   [],
         pINT   &
);
}

void doflush() { fflush(0); mexEvalString("pause(0);"); };

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   if (nargin!=1 || !mxIsDouble(argin[0])) {
   printf("\n  ERR invalid usage\n"); exit(0); }

   const mwSize *S=mxGetDimensions(argin[0]);
   const double *A0=mxGetPr(argin[0]);
   pINT i=0, m=S[0], n=S[1], N=m*n, lwork=N, ipiv[n]; 
   size_t s=N*sizeof(double);

   double *A=NULL, *work = new double[lwork];
   if (nargout) {
      argout[0]=mxCreateDoubleMatrix(m,n,mxREAL); if (argout[0]) {
      A=mxGetPr(argout[0]); }
   }
   else {
      printf("\n  ERR no output argument required - return.\n");
      mexErrMsgTxt("");
   }

   if (!A || !work) {
      printf("\n  ERR failed to allocate data (%lX, %lX)\n",A,work);
      mexErrMsgTxt("");
   }
   work[0]=0; memcpy(A,A0,s);


   dgetrf(m,n,A,m,ipiv,i);
   if (i) { 
      printf("\n  ERR dgetrf() returned e=%d\n",-i);
      mexErrMsgTxt("");
   }

   dgetri(n,A,m,ipiv,work,lwork,i);
   if (i) {
      if (i<0)
           printf("\n  ERR dgetri() got invalid argument #%d\n",-i);
      else printf("\n  ERR dgetri() got singular matrix (i=%d)\n",i);
      mexErrMsgTxt("");
   }

   delete [] work;
};

