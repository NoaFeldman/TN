char USAGE[] =
/* ==================================================================
 * mpsGetDim.cc */

"   Usage: size = traceQS(A, i1, i2)                              \n\
                                                                  \n\
       trace out index set of given QSpace object                 \n\
                                                                  \n\
   AWb © Nov 2006                                                 \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

template<class TQ, class TD>
mxArray* traceQS(
    const QSpace<TQ,TD> &A,
    unsigned i1,
    unsigned i2,
    unsigned narg
);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i1=0,i2=0,r=-1;

    if (!nargin || nargout>1) usage(FL,
    "ERR Invalid number of I/O arguments.");

    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
    "input not valid QSpace object");

    if (nargin>1) mxGetNumber(argin[1], i1);
    if (nargin>2) mxGetNumber(argin[2], i2);

 
    if (mxIsQSpace(argin[0])) {
       const QSpace<double,double> A(argin[0],'r');
       argout[0]=traceQS(A,i1,i2,nargin);
    }
    else {
       const QSpace<double,wbcomplex> A(argin[0]);
       argout[0]=traceQS(A,i1,i2,nargin);
    }
}

template<class TQ, class TD>
mxArray* traceQS(
    const QSpace<TQ,TD> &A,
    unsigned i1,
    unsigned i2,
    unsigned narg
){
   unsigned r;
   QSpace<TQ,TD> A2;

   if (A.isEmpty()) {
      A2.init2scalar(0);
      return A2.toMx();
   }

   r=A.rank(); 

   if (narg==1 && r==2) { i1=1; i2=2;  } else
   if (narg!=3) wblog(FL,"ERR Invalid number of I/O arguments.");

   A.contract(r-i1+1,r-i2+1,A2);

   return A2.toMx();
}

