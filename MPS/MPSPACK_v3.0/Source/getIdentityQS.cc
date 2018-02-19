char USAGE[] =
/* ================================================================ */
"   Usage:                                                      \n\
                                                                \n\
   [A,CG]=getIdentityQS([opts], Op1 [,i1] [,Op2[,i2]] [,tag]);  \n\
                                                                \n\
     if Op2 [,i2] is not specified                              \n\
        get plain identity operator from given QSpace.          \n\
     if Op2 [,i2] is specified:                                 \n\
        get tensor product space of input spaces defined by     \n\
        Op1 and Op2. If i[12] is not specified, rank-2 objects  \n\
        are assumed taking (the first) two indizes.             \n\
                                                                \n\
     An itag name for the combined output space (dim3) may be   \n\
     specified as last argument (with conj flag being ignored). \n\
                                                                \n\
   Options                                                      \n\
                                                                \n\
     '-h'  display this usage and exit                          \n\
     '-v'  verbose                                              \n\
                                                                \n\
   (C) AW : Apr 2010 ; Oct 2014                                 \n";

// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getIdentityQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
   char e=0, vflag=1, zflag=0, isra, isrb=-1;
   int k=-1, narg=nargin;
   wbindex Ia,Ib; iTag t;
   mxArray *a=NULL;

   if (nargin && isHelpIndicator(argin[0])) { usage(); return; }

   if (nargin==1 && mxIsEqual(argin[0],"--getCG")) {
      argout[0]=gCG.toMx();
      if (nargout>1) { argout[1]=gRG.toMx(); }
      return;
   }

   if (nargin && mxIsEqual(argin[0],"-v")) {
      vflag='V'; ++argin; --narg;
   }
   if (narg && mxIsChar(argin[narg-1])) { --narg;
      if (mxIsEqual(argin[narg],"-0")) zflag=1;
      else t.init(FL,argin[narg]).deConj();
   }

   if (narg<1 || nargout>2) e=__LINE__; else
   if (narg>1) {
      if (mxIsIndex(argin[1])) { 
         Ia.init(argin[1],1);
         narg-=2; if (narg) k=2; 
      } else { narg-=1; if (narg) k=1; }
   }

   if (narg) {
      if (narg>2) e=__LINE__; else
      if (!e && narg>1) {
         if (!mxIsIndex(argin[k+1])) e=__LINE__; else
         Ib.init(argin[k+1],1);
      }
   }

   if (e) {
      if (nargin || nargout) wblog(__FILE__,e,"ERR invalid usage");
      else { usage(); return; }
   }


   isra=mxIsQSpaceVec(FL,argin[0]);
      if (!isra && !mxIsQSpaceVec(FL,argin[0],-1,0,0,'C')) wblog(FL,
      "ERR invalid input operator Op1%N%N%s",str);

   if (k>0) {
      isrb=mxIsQSpaceVec(FL,argin[k]);
      if (!isrb && !mxIsQSpaceVec(FL,argin[k],-1,0,0,'C'))
      wblog(FL,"ERR invalid input operator Op2%N%N%s",str);
   }

   if (isra) {
      wbvector< QSpace<gTQ,double> > A;
      mxInitQSpaceVec(FL,argin[0],A,'r');

      if (isrb>0) {
         wbvector< QSpace<gTQ,double> > B;
         QSpace<gTQ,double> C;

         mxInitQSpaceVec(FL,argin[k],B,'r');
         C.initIdentityCG(A,Ia,B,Ib,vflag); if (t!=0) C.SetTag(3,t,'l');
         a=C.save2Mx(vflag);
      }
      else if (isrb==0) {
         wbvector< QSpace<gTQ,wbcomplex> > B;
         QSpace<gTQ,wbcomplex> C;

         mxInitQSpaceVec(FL,argin[k],B);
         C.initIdentityCG(A,Ia,B,Ib,vflag); if (t!=0) C.SetTag(3,t,'l');
         a=C.save2Mx(vflag);
      }
      else {
         QSpace<gTQ,double> C;
         if (!Ib.isEmpty()) wblog(FL,"ERR invalid usage (%d)",Ib.len);
         if (t!=0) wblog(FL,
            "ERR %s() got tag specified (%s) !?",FCT,t.toStr().data);
         a=C.initIdentityCG(A,Ia,zflag).save2Mx();
      }
   }
   else {
      wbvector< QSpace<gTQ,wbcomplex> > A;
      mxInitQSpaceVec(FL,argin[0],A,'r');

      if (isrb>0) {
         wbvector< QSpace<gTQ,double> > B;
         QSpace<gTQ,double> C;

         mxInitQSpaceVec(FL,argin[k],B,'r');
         C.initIdentityCG(A,Ia,B,Ib,vflag); if (t!=0) C.SetTag(3,t,'l');
         a=C.save2Mx(vflag);
      }
      else if (isrb==0) {
         wbvector< QSpace<gTQ,wbcomplex> > B;
         QSpace<gTQ,wbcomplex> C;

         mxInitQSpaceVec(FL,argin[k],B);
         C.initIdentityCG(A,Ia,B,Ib,vflag); if (t!=0) C.SetTag(3,t,'l');
         a=C.save2Mx(vflag);
      }
      else {
         QSpace<gTQ,wbcomplex> C;
         if (!Ib.isEmpty()) wblog(FL,"ERR invalid usage (%d)",Ib.len);
         if (t!=0) wblog(FL,
            "ERR %s() got tag specified (%s) !?",FCT,t.toStr().data);
         a=C.initIdentityCG(A,Ia,zflag).save2Mx();
      }
   }

   argout[0]=a;
   if (nargout>1) { argout[1]=gCG.toMx(); }
};


