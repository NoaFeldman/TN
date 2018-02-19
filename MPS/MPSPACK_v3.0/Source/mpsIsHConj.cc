char USAGE[]=
/* ================================================================== */
"   Usage: [i,istr]=mpsIsHConj(A [,eps]);                         \n\
                                                                  \n\
   check whether QSpace is hermitian within threshold eps (1E-12).\n\
   A must be even-rank object as hyperindex is allowed.           \n\
                                                                  \n\
Wb,Aug17,06                                                       \n\
";
// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsGetDim"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){

    QSpace<gTQ,gTD> A;
    char vflag=0; int i; double eps=1E-12;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (nargin<1) usage(FL,"ERR Invalid number of I/O arguments.");

    if (!mxIsQSpace(FL,argin[0],'c'))
    wblog(FL,"ERR First argument requires valid QSpace.");

    A.init(FL,argin[0],'r');

    if (nargin>1) {
       wbstring mark(nargin);
       for (i=1; i<nargin; ++i) {
          if (mxIsChar(argin[i])) {
             wbstring s(argin[i]);
             if (s=="-v") vflag='v';
             else wberror(FL,str);
          }
          else if (mxIsNumber(0,0,argin[i])) {
             if ((++mark[i])>1 || mxGetNumber(argin[1], eps))
             wberror(FL,str);
          }
       }
    }

    i=A.isHConj(0,0,eps,vflag);
    argout[0]=numtoMx(i);

    if (nargout>1) argout[1]=mxCreateString(i ? "" : str);
};

