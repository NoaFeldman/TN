char USAGE[] =
/* ================================================================== */
"   Usage: A = permuteQS(A,P [,'conj'])                            \n\
                                                                   \n\
       permute input QSpace using given permutation P.             \n\
                                                                   \n\
       Optional trailing 'conj' also applies (complex) conjugation \n\
       (note that this also affects real QSpaces in that qdir and  \n\
       itags are altered!).                                        \n\
                                                                   \n\
       For convenience, P [,'conj'] may also be represented as     \n\
       single string, e.g. [2 1],'conj' is equivalent to '2,1;*'   \n\
       or '2 1;*' where the convention on string notation          \n\
       follows that of contraction indices [ctrIdx].               \n\
                                                                   \n\
   AW (C) Aug 2006 ; May 2010 ; Oct 2014                           \n";

// This is a MEX wrapper routine for MATLAB.
/* ================================================================== */
  
// Wb,Dec09,14 :: added optional conj as trailing arument or string.


#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "permuteQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"


template<class TD>
mxArray* PERMUTE_QS(
   const char *F, int L,
   const QSpace<gTQ,TD> &A, const wbperm &P, char conj=0
){
   QSpace<gTQ,TD> B; unsigned r=A.rank(F_L);

   if (P.len!=r) wblog(F,L,
      "ERR invalid permutation (len=%d/%d)",P.len,r);
   A.permute(P,B);
   if (conj) { B.Conj(); B.SortDegQ(); }

   return B.toMx();
};


void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   wbperm P; unsigned r=-1;
   char isr=1, conj=0;

   if (nargin && isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin<2 || nargin>3) usage(FL,"ERR invalid number of I/O arguments.");

   isr=mxIsQSpace(argin[0]);
   if (!isr && !mxIsQSpace(FL,argin[0],r,'c')) wblog(FL,
      "ERR invalid QSpace as first argument");

   if (mxIsChar(argin[1])) {
      ctrIdx q(FL,argin[1]);
      P.init(FL,q); conj=q.conj;

      if (nargin>2) usage(FL,"ERR invalid usage "
        "(got too many input arguments for char permutation)"
      );
   }
   else {
      P.init(FL,argin[1], 1);
      if (!P.isValidPerm()) { P.print("P");
         wblog(FL,"ERR invalid permutation (2nd arugment)");
      }
      if (nargin>2) {
         if (mxIsEqual(argin[2],"conj")) { conj=1; }
         else usage(FL,"ERR invalid last argument (expecting 'conj')");
      }
   }

   if (isr) {
      const QSpace<gTQ,double> A(argin[0],'r');
      argout[0]=PERMUTE_QS(FL,A,P,conj);
   }
   else {
      const QSpace<gTQ,wbcomplex> A(argin[0]);
      argout[0]=PERMUTE_QS(FL,A,P,conj);
   }
}

