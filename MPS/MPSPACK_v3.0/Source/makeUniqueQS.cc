char USAGE[] =
/* ==================================================================
 * mpsGetDim.cc */

"   Usage: size = mpsMakeUnique(A [,opts])                        \n\
                                                                  \n\
       MEX wrapper routine for QSpace::makeUnique()               \n\
                                                                  \n\
   Options                                                        \n\
                                                                  \n\
       -[qQ]    quite mode                                        \n\
       -[vV]    verbose mode (wrt. outer multiplicity)            \n\
                                                                  \n\
   AWb (C) Nov 2007 ; 2012                                        \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsMakeUnique"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    wbvector<unsigned> D;
    unsigned r=-1;

    wbMatrix<gTQ> Q;
    wbvector<unsigned> d;
    OPTS opts; char qflag=-1;

    str[0]=0;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (!nargin || nargout>1)
    usage(FL,"ERR Invalid number of I/O arguments.");

   opts.init(argin+1,nargin-1);

   if (opts.getOpt("-q")) qflag='q'; else {
   if (opts.getOpt("-Q")) qflag='Q'; else {
   if (opts.getOpt("-v")) qflag='v'; else {
   if (opts.getOpt("-V")) qflag='V'; }}}

   opts.checkAnyLeft();

 
    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
    "input not valid QSpace object");

    if (mxIsQSpace(argin[0])) {
       QSpace<gTQ,double> A(argin[0]);
       A.makeUnique(qflag);
       argout[0]=A.toMx();
    }
    else {
       QSpace<gTQ,wbcomplex> A(argin[0]);
       A.makeUnique(qflag);
       argout[0]=A.toMx();
    }
};

