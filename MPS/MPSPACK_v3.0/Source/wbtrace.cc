char USAGE[] =
/* ==================================================================
 * contractTrace.cc */

"   Usage: i = wbtrace(A, [i1 i2; i3 i4; ...] [,rank]);          \n\
                                                                 \n\
       generalized trace for arbitrary dimensional object        \n\
       contracts i1 with i2, i3 with i4, ... (1-based)           \n\
                                                                 \n\
       optional 3rd argument to specify rank if singleton        \n\
       dimensions exist.                                         \n\
                                                                 \n\
   AWb C Jul 2006                                                \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

template<class T>
mxArray* MEX_FUNCTION(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[],
    wbarray<T> &A,
    wbarray<T> &Aout
);


void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin<2 || nargin>3 || nargout>1)
   usage(FL,"ERR Invalid number of I/O arguments.");

   if (mxMatSafeGuard(0,1,argin[0],-1,'c')) wblog(FL,
      "ERR invalid usage (double input array expected)\n%s", str);

   if (!mxIsComplex(argin[0])) {
      wbarray<double> A(argin[0]), Aout;
      argout[0]=MEX_FUNCTION(nargout,argout,nargin,argin,A,Aout);
   }
   else {
      wbarray<wbcomplex> A(argin[0],0), Aout;
      argout[0]=MEX_FUNCTION(nargout,argout,nargin,argin,A,Aout);
   }
}


template<class T>
mxArray* MEX_FUNCTION(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   wbarray<T> &A,
   wbarray<T> &Aout
){
   unsigned i,j,m,n; double k;

   wbMatrix<double> D12;
   wbMatrix<unsigned> I12;
   wbMatrix<unsigned> mark;

   if (!A.isEmpty()) {
      D12.init(argin[1]); I12.initDim(D12);

      m=I12.dim1; n=A.SIZE.len;
      mark.init(m,2);

      if (nargin>2) {
         double r;

         if (!mxIsDblScalar(argin[2])) wblog(FL,
            "ERR invalid arg #3 (rank)");
         if (mxGetNumber(argin[2],r)) wblog(FL,"ERR %s",str); 

         if (r<n) wblog(FL,"ERR invalid rank (%g/%d)",r,n);
         else if (r>n) {
            A.addSingletons((unsigned)r);
            n=A.SIZE.len;
         }
      }

      for (i=0; i<m; i++)
      for (j=0; j<2; j++) { k=D12(i,j);
         if (k<1 || k>n) wblog(FL,
            "ERR trace index out of bounds: %g/%d (%d,%d; %d)\n"
            "Hint: use 3rd argument to specify rank.",
             k,n, i,j, mark(i,j)
         );
         else if (k!=double(unsigned(k)) || mark(i,j)++) wblog(FL,
            "ERR Invalid trace index set (%d,%d: %g/%d; %d).",
             i,j,D12(i,j),n,mark(i,j));
         else

         I12(i,j)=unsigned(D12(i,j))-1;
      }

      A.contract(I12,Aout);
   }
   else {
      Aout.init(1,1); Aout.data[0]=0;
   }

   return Aout.toMx();
};


