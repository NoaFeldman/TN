char USAGE[] =
/* ================================================================ */
"   Usage: nrm = mpsNormQS(A)                                     \n\
                                                                  \n\
       calculates the (Frobenius) norm of given QSpace, i.e.      \n\
       |A|^2 = norm(A)^2 = trace(A'*A), and appropriately         \n\
       generalized to arbitray rank. In the presence of CGCs, this\n\
       implies |A|^2 = sum_i (|data{i}|^2 * prod_j(|cgs(i,j)|^2)).\n\
                                                                  \n\
    (C) AW : Jun 2010 ; Oct 2014                                  \n\
";
// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsNormQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    double nrm2; unsigned r=-1;

    if (nargin==0 || isHelpIndicator(argin[0])) { usage(); return; }

    if (nargin!=1) usage(FL,"ERR invalid number of I/O arguments.");
    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,"valid QSpace required");

    if (mxIsQSpace(argin[0])) {
       const QSpace<gTQ,double> A(argin[0],'r');
       nrm2=A.norm2();
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       wbcomplex z2=A.norm2();
       if (fabs(z2.i)>1E-14) wblog(FL,
          "ERR %s() got imaginary part (%s)",FCT,z2.toStr().data);
       nrm2=z2.r;
    }

    argout[0]=numtoMx(sqrt(nrm2));
};

