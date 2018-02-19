char USAGE[] = 
/* =============================================================
 * mat2cellc.cc */

"   Usage: [A [,I,D]] = uniquerows(Q [,opts])                  \n\
                                                               \n\
      MEX file that makes rows of array unique. Moreover, it   \n\
      returns a cell array of index sets of rows with the same \n\
      data. D specifies the multiplicity of the unique rows.   \n\
                                                               \n\
   Options                                                     \n\
                                                               \n\
     '-1'  return only first index in each group in I          \n\
                                                               \n\
   Wb (C) Sep 2006                                             \n";

/* This is a MEX-file for MATLAB.
 * ============================================================= */

#include "wblib.h"

/* gateway routine  */

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   wbvector<unsigned> D;
   wbperm P; mxArray *a=0;
   char uflag=0; // Wb,Jun01,11

   if (nargin && isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin>1) {
      wbstring s(argin[1]);
      if (s=="-1") { uflag=1; --nargin; }
   }

   if (nargin!=1) usage(FL,"Invalid number of input arguments");
   if (nargout>3) wberror(FL,"Invalid number of output arguments");

   if (mxIsDouble(argin[0])) {
      wbMatrix<double> Q; Q.init(argin[0]);
   // sort to original index order on degenerate space
      Q.groupRecs(P,D);
      a=Q.toMx();
   }
   else if (mxIsChar(argin[0])) {
      wbMatrix<char> Q; Q.init(argin[0]);
   // sort to original index order on degenerate space
      Q.groupRecs(P,D);
      a=Q.toMx();
   }
   else {
      wblog(FL,"ERR invalid type of input argument #1 (%s)",
      mxGetClassName(argin[0]));
   }

   argout[0]=a; // MX_USING_ANS

   if (nargout>1) {
      if (uflag) {
         unsigned i,l,d; wbindex I(D.len);
         for (l=i=0; i<D.len; i++, l+=d) { d=D[i]; I[i]=P[l]; }
         argout[1]=I.toMx();
      }
      else {
         unsigned i,l,d;
         wbvector<unsigned> I;
         mxArray *S = mxCreateCellMatrix(1,D.len);

         for (l=i=0; i<D.len; i++, l+=d) {
            d=D[i]; I.init2ref(d, P.data+l); I+=1; // make index 1-based
            mxSetCell(S,i,I.toMx());
         }
         argout[1]=S;
      }
   }

   if (nargout>2) argout[2]=D.toMx();

   return;
}

