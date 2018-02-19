char USAGE[] =
/* ==================================================================
 * mpsPlusQS.cc */

"   Usage: C = mpsPlusQS(A, B [, bfac, opts])                    \n\
                                                                 \n\
       adds one QSpace to another as in C = A+bfac*B ([bfac=1])  \n\
       block-diagonal identity matrix.                           \n\
                                                                 \n\
   Options                                                       \n\
                                                                 \n\
      '-v'  verbose mode (e.g. on remaining outer multiplicity)  \n\
                                                                 \n\
   (C) AW : Jun 2006 ; Oct 2014                                  \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsPlusQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

template<class TQ, class TD>
inline QSpace<TQ,TD>& MPS_PLUS_QS(
   QSpace<TQ,TD> &C, const mxArray *argin[], int nargin, char ref=0
){
   const QSpace<TQ,TD> A(FL,argin[0],ref);
   const QSpace<TQ,TD> B(FL,argin[1],ref);
   


   double bfac=1; str[0]=0;
   unsigned l=nargin-1; char vflag=0;
   if (l>=2 && mxIsChar(argin[l])) { str[0]=0;
      mxGetString(argin[l],str,128);
      if (!strcmp(str,"-v")) { vflag='v'; --l; --nargin; }
      else wblog(FL,
         "ERR %s() invalid input option '%s' [%d]",str,FCT,l+1
      );
   }
   if (nargin>2) if (mxGetNumber(argin[2], bfac)) wblog(FL,str);

   C.Cat(FL,A,B,TD(1),TD(bfac),
      vflag ? vflag : 1
   );

   C.SkipZeroData();
   C.NormCG();

   return C;
};


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    char isra, isrb;

    if (nargin==0 || isHelpIndicator(argin[0])) { usage(); return; }

    if (nargin<2 || nargin>4)  wblog(FL,
       "ERR invalid number of arguments (%d)",nargin);


    isra=mxIsQSpace(argin[0]);
    isrb=mxIsQSpace(argin[1]);

    if (isra && isrb) {
       QSpace<gTQ,double> C;
       MPS_PLUS_QS(C,argin,nargin,'r');
       argout[0]=C.toMx();
    }
    else {
       QSpace<gTQ,wbcomplex> C;
       MPS_PLUS_QS(C,argin,nargin);
       argout[0]=C.toMx();
    }
};


