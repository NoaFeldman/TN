char USAGE[] =
/* ==================================================================
 * mpsGetDim.cc */

"   Usage: [q,dd,dc] = mpsGetQDim(A,dim)                          \n\
                                                                  \n\
       get dimension for every Q in QIDX for given QSpace object  \n\
       dim is dimension to check (dim=='op' can be used to        \n\
       consider both dimensions of an operator simultaneously)    \n\
                                                                  \n\
   Data returned                                                  \n\
                                                                  \n\
       q    symmetry labels present (unique; m x nsymlabels)      \n\
       dd   corresponding dimension of multiplet data (m)         \n\
       dc   corresponding dimension of Clebsch-Gordan (m x nsym)  \n\
                                                                  \n\
   AW (C) Jun 2007 ; Wb,Oct12,10                                  \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsGetQDim"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    wbvector<INDEX_T> D;
    unsigned k=0, r=-1; char isop=0;

    wbMatrix<gTQ> Q;
    wbvector<INDEX_T> dd;
    wbMatrix<INDEX_T> dc;

    str[0]=0;

    if (nargin!=2 || nargout>3) usage(FL,
       "ERR invalid number of I/O arguments.");

 
    try {
       if (!mxIsQSpace(FL,argin[0],r,'c'))
       wblog(FL,"ERR Input not valid QSpace object.");
    }
    catch (...) {
       wblog(FL,"ERR Input not valid QSpace object.");
    }

    if (mxGetNumber(argin[1],k,'q')==0) {
       if (k==0 || (int)k>(int)r) wblog(FL,
          "ERR 2nd argument out of bounds (%d/%d)",k,r);
       k--;
    }
    else {
       if (!mxIsChar(argin[1])) wberror(FL,str);
       if (mxGetString(argin[1],str,12))
       wblog(FL,"ERR failed to read string (arg #2) ???");

       if (strcmp(str,"op")) wberror(FL,"invalid 2nd argument");
       isop=1;
    }

    if (mxIsQSpace(argin[0])) {
       const QSpace<gTQ,double> A(argin[0],'r');
       if (isop)
            A.getQDim(  Q,dd,&dc);
       else A.getQDim(k,Q,dd,&dc);
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0],'r');
       if (isop)
            A.getQDim(  Q,dd,&dc);
       else A.getQDim(k,Q,dd,&dc);
    }

    argout[0]=Q.toMx();

    if (nargout>1) argout[1]=wbvector<double>(dd).toMx();
    if (nargout>2) argout[2]=wbMatrix<double>(dc).toMx();
};

