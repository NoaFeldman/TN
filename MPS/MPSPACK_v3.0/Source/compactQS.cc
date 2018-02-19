char USAGE[] =
/* ================================================================ */
"   Usage:                                                      \n\
                                                                \n\
   S=compactQS([opts],'qtype',Q1,Q2,Q,A [,perm]);               \n\
                                                                \n\
      reduce matrix elements of the rank-3 IROP.                \n\
      Here A is still specified in full (dense) tensor format   \n\
      with matrix elements given by                             \n\
                                                                \n\
         A(i1,i2,io) := <i1|A(io)|i2>                           \n\
                                                                \n\
      NB! this assumes the index (1,2,op) specified on the lhs! \n\
      The symmetry types are specified by qtype. Moreover,      \n\
      all states i1 and i2 must be already cast into symmetry   \n\
      eigenstate grouped into multiplets (see getSymStates.cc)  \n\
      with q- and z-labels both present and interleaved in      \n\
      each input state label sets Q*.                           \n\
                                                                \n\
      The state spaces are combined in the order (Q2+Q =>Q1),   \n\
      thus using the Clebsch Gordan coefficients (Q1,Q|Q2).     \n\
      This is a natural order in that, Q=-1 annihilates a       \n\
      particle on a |ket> state (the annihilation operator      \n\
      on a ket state increases charge, i.e. Q would be +1,      \n\
      in case Q1 and Q had been combined into Q!).              \n\
      This therefore represents a natural index order.          \n\
                                                                \n\
   Options                                                      \n\
                                                                \n\
      perm  permutation to be applied to {Q1,Q2,Q} and A        \n\
            prior to combining first two indizes into third.    \n\
                                                                \n\
      '-h'  display this usage and exit                         \n\
      '-v'  verbose                                             \n\
                                                                \n\
   (C) AW : Apr 2010 ; Jul 2011 ; Oct 2014                      \n";

// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "compactQS"
#endif


#define LOAD_CGC_QSPACE


#include "wblib.h"



void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
   unsigned i,k=0;
   char vflag=0;


   if (nargin && isHelpIndicator(argin[k])) { usage(); return; }

   if (nargin && mxIsEqual(argin[k],"-v")) { vflag=1; } else
   if (nargin && mxIsEqual(argin[k],"-V")) { vflag=2; }
   if (vflag) { nargin--; k++; if (vflag>1) { CG_VERBOSE=1; } }

   if (nargout>2 || nargin<5 || nargin>6 || !mxIsChar(argin[k]) ||
      !mxIsDblMat(0,0,argin[k+1]) || !mxIsDblMat(0,0,argin[k+2]) ||
      !mxIsDblMat(0,0,argin[k+3]) || !mxIsDblArr(0,0,argin[k+4])
   ){
       if (nargin || k || nargout) wblog(FL,"ERR invalid usage");
       else { usage(); return; }
   }

   QVec qvec(FL,argin[k++]);
   wbMatrix<gTQ> Q1(FL,argin[k++]), Q2(FL,argin[k++]), Q(FL,argin[k++]);
   wbarray<double> D3(FL,argin[k++]); D3.addSingletons(3);

   QSpace<gTQ,double> A;
   wbperm P; if (nargin>5) P.init(FL,argin[k++],1);
   wbvector< wbMatrix<gTQ>* > qq(3);

   qq[0]=&Q1; qq[1]=&Q2; qq[2]=&Q;
   if (!P.isEmpty()) {
      if (P.len!=3) wblog(FL,"ERR invalid permutation length (%d/3)",P.len);
      qq.Permute(P); D3.Permute(P);
   }
   
   for (i=0; i<3; ++i)
   Wb::FixRational(FL,qq[i]->data,qq[i]->numel(),'r',4,1024);


   A.initOpZ_WET(FL,qvec,*(qq.data[0]),*(qq.data[1]),*(qq.data[2]),D3);

   if (Q.allEqual(0)) {

      wbMatrix<gTQ> QQ(A.QIDX); QQ.getBlocks(0,1,A.QDIM,A.QIDX);
      wbarray<double> **dd=A.DATA.data;
      CRef<gTQ> *cg=A.CGR.data;
      unsigned n=A.CGR.numel(), ok=1;

      for (i=0; i<A.DATA.len; ++i) {
         if (dd[i]->SIZE.len!=3 || dd[i]->SIZE.data[2]!=1) {
            wblog(FL,
              "WRN failed to reduce to scalar op (%d)", dd[i]->SIZE.len,
               dd[i]->SIZE.len==3 ? dd[i]->SIZE.data[2] : -1);
            ok=0; break;
         }
      }

      if (ok) {
      for (i=0; i<n; ++i) {
         if (!cg[i].cgr) { wblog(FL,
            "WRN %s() got cgr==NULL (%d)",FCT,cg[i].cgw.len); 
            continue;
         }

         const CRef<gTQ> &ci=cg[i];
         if (ci.cgw.len!=1 || ci.cgr->qdir.len!=3) wblog(FL,
            "ERR %s() invalid CRef data\n%s",FCT,ci.toStr().data);
         if (ci.Size(2)!=1) wblog(FL,
            "ERR failed to reduce to scalar op\n%s having P=[%s]",
            ci.toStr().data, ci.cgp.toStr().data
         );
      }

      for (i=0; i<A.DATA.len; ++i) dd[i]->SIZE.len=2;

      for (i=0; i<n; ++i) {
         cg[i].Reduce2Identity(
            A.qtype.allAbelian() ? 0 : '!'
         );
      }

      A.otype=QS_NONE;
      A.itags.len=2;

   }}

   A.NormCG();

   argout[0]=A.toMx();
   if (nargout>1) argout[1]=gCG.toMx();

   if (vflag>1 && nargout<2) {
      mxPutAndDestroy(0,0,gCG.toMx(),"gCG","base");
   }

#ifdef __WBDEBUG__
   Wb::MemCheck(FL);

   qvec.init(); P.init(); qq.init();
   Q1.init(); Q2.init(); Q.init(); D3.init();
   gCG.init(); A.init();

   Wb::MemCheck(FL);
#endif
}

