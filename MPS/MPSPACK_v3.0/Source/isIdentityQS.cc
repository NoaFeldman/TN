char USAGE[] =
/* ==================================================================
 * mpsIsIdentityQS.cc */

"   Usage: i = mpsIsIdentityQS(A [,eps])                         \n\
                                                                 \n\
       check whether QSpace A represents a block-diagonal        \n\
       identity matrix. This also checks the corresponding       \n\
       Clebsch-Gordon spaces if any (default: eps=1E-12).        \n\
                                                                 \n\
   (C) AW : May 2006 ; Oct 2010 ; Oct 2014                       \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsIsIdentityQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,r=-1; double eps=1E-12;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (nargin<1) usage(FL,"ERR Invalid number of I/O arguments.");

 
    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
       "first argument requires valid QSpace");

    if (nargin>1)
    if (mxGetNumber(argin[1], eps)) wberror(FL,str);

    if (mxIsQSpace(argin[0])) {
       const QSpace<gTQ,double> A(argin[0],'r');
       i=(A.isIdentityMatrix(eps) ? 1 : 0);
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       i=(A.isIdentityMatrix(eps) ? 1 : 0);
    }

    argout[0]=numtoMx(i);
}

