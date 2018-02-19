char USAGE[] =
/* ==================================================================
 * mpsIsIdentity.cc */

"   Usage: i = mpsIsIdentityCG(A)                                 \n\
                                                                  \n\
       check whether QSpace A has indentity Clebsch-Gordan spaces.\n\
       Note that for full abelian spaces, this is always true.    \n\
                                                                  \n\
   Wb,Jan12,12                                                    \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsIsIdentityCG"
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
       i=A.hasIdentityCGS(eps) ? 1 : 0;
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       i=A.hasIdentityCGS(eps) ? 1 : 0;
    }

    argout[0]=numtoMx(i);
}

