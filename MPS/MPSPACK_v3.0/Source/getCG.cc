
char USAGE[] =
/* =============================================================
 * clebsch.cc */

"   Usage:                                                    \n\
                                                              \n\
   S=getCG([opts], 'qtype',q1,q2 [,q]);                       \n\
                                                              \n\
      get CData for given symmetries. If q is not specified,  \n\
      all possible q's resulting from q1 and q2 are returned. \n\
                                                              \n\
   Options                                                    \n\
                                                              \n\
     '-h'  display this usage and exit                        \n\
     '-v'  verbose                                            \n\
                                                              \n\
   (C) Wb,Oct12,09 ; Wb,Apr26,10                              \n";


#define gTQ  double
#define QTYPE gTQ

#include "wblib.h"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
   unsigned k=0; mxArray *a=NULL;
   char vflag=0;

   if (nargin && isHelpIndicator(argin[k])) { usage(); return; }
   if (nargin && mxIsEqual(argin[k],"-v")) { vflag=1; nargin--; k++; }

   if (nargin<3 || nargin>4 || !mxIsChar(argin[k]) ||
      !mxIsVector(argin[k+1]) || !mxIsVector(argin[k+2]) ||
     (nargin>3 && !mxIsVector(argin[k+3]))
   ){
       if (nargin || k || nargout) wblog(FL,"ERR invalid usage");
       else { usage(); return; }
   }

   QVec qvec(FL,argin[k]);
   qset<gTQ> J1(FL,argin[k+1],&qvec), J2(FL,argin[k+2],&qvec), J;

   if (nargin>3) J.init(FL,argin[k+3],&qvec);

   if (J1.len<=1) {
      if (J1.len==0 || J2.len!=J1.len || qvec.len!=J1.len)
      wblog(FL,"ERR invalid usage (%d/%d/%d)",J1.len,J2.len,qvec.len);

      if (J.isEmpty()) {
         wbvector < CData<gTQ,unsigned,double>* > S;
         gCG.getCData(FL,qvec[0],J1,J2,S);
         a=S.toMxP();
      }
      else {
         const CData <gTQ,unsigned,double> &S=gCG.getCData(FL,qvec[0],J1,J2,J);
         a=S.toMx();
      }
   }
   else wblog(FL,
   "ERR multiple Q-labels not implemented yet (%d)",J1.len);

   if (vflag) gCG.printStatus(FL);

   argout[0]=a;
}

