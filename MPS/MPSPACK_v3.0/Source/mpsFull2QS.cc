char USAGE[]=
/* ================================================================== */
"[M,I,QF] = mpsFull2QS(A [,R])                             \n\
                                                           \n\
   convert QSpace to full 2D array                         \n\
   with the indices split into two groups of equal size.   \n\
                                                           \n\
   Q can be enlarged using the optional QSpace R whos      \n\
   data will be used only to determine the size for newly  \n\
   acquired entries (extra space is initialized to zero).  \n\
   Optional 3rd output argument QF is I data expanded.     \n\
                                                           \n\
Wb,Sep05,06                                                \n\
";
// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

#include "wblib.h"

#define QTYPE double
#define _TQ QTYPE

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

   wbMatrix<_TQ> QB;
   wbMatrix<unsigned> SB;
   mxArray *a;
   unsigned r=-1; int i; str[0]=0;

   if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin<1 || nargin>2) wblog(FL,
   "ERR Invalid number of I/O arguments.");

   for (i=0; i<nargin; i++) {
      try { mxIsQSpace(FL,argin[i],r,'c'); }
      catch (...) {
         wblog(FL,"ERR input arg #%d must be of type QSpace\n%s",
         i+1, i ? "and of SAME rank as arg #1" : "");
      }
   }

   for (i=0; i<nargin; i++) if (!mxIsQSpace(argin[0])) break;
   if (i==nargin) {

      QSpace<_TQ,double> A(argin[0],'r');
      wbarray<double> M;

      if (nargin>1) {
         QSpace<_TQ,double> Aref(argin[1],'r');
         A.ExtendQ(Aref);
      }
      A.toFull2(M,QB,SB); a=M.toMx();
   }
   else {

      QSpace<_TQ,wbcomplex> A(argin[0]);
      wbarray<wbcomplex> M;

      if (nargin>1) {
         QSpace<_TQ,wbcomplex> Aref(argin[1]);
         A.ExtendQ(Aref);
      }
      A.toFull2(M,QB,SB); a=M.toMx();
   }

   argout[0]=a;

   if (nargout>1) {
      mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

      mxAddField2Scalar(FL,S, "QB", QB.toMx());
      mxAddField2Scalar(FL,S, "SB", SB.toMx());

      argout[1]=S;
   }

   if (nargout>2) {
      unsigned q,l,n; i=0;
      wbvector<unsigned> S;
      wbMatrix<_TQ> QF;

      SB.recProd(S); n=S.sum();
      QF.init(n,QB.dim2);

      for (q=0; q<QB.dim1; q++) { n=S[q];
      for (l=0; l<n; l++, i++) QF.recSet(i,QB,q); }

      argout[2]=QF.toMx();
   }
}

