char USAGE[] =
/* ======================================================================== */
"   Usage: [[NRG_DATA,] info] = ...                                        \n\
            NRGWilsonQS(H0, A0, Lambda, ff, FC, Z, [, gg, FL, OPTS])       \n\
                                                                           \n\
   NRG run using (non-)Abelian symmetries for a Wilson chain of length L   \n\
   based on the standard iterative diagonalization prescription.           \n\
                                                                           \n\
      H0      impurity / starting Hamiltonian (rank-2 QSpace)              \n\
              H0 will be diagonalized upon initialization hence does not   \n\
              have to be diagonal on input. The resulting basis            \n\
              transformation is multiplied onto A0.                        \n\
                                                                           \n\
      A0      basis used for H0 in LRs order (rank-3 QSpace)               \n\
              where the local space s(igma) must already have the          \n\
              (equivalent) space of a Wilson site (note that A0 is used    \n\
              to obtain the matrix elements of FC at the first iteration.  \n\
                                                                           \n\
      Lambda  used to rescale H_k by sqrt(Lambda); therefore this should be\n\
              consistent with the Lambda used to obtain ff and gg below    \n\
                                                                           \n\
      ff      strength of couplings along Wilson chain as they appear in   \n\
              the physical Hamiltonian, i.e. they must not be rescaled yet!\n\
              (double, (L-1) x Nf)                                         \n\
                                                                           \n\
      FC      set of Nf annihilation operators used for coupling (QSpace)  \n\
              e.g. in case of all-abelian symmetries (default) this becomes\n\
              H = F1'*ZF2 + F1*(ZF2)' + h.conj. (zflag=1, see below)       \n\
              with F1=F2=FC such that operator appears twice together      \n\
              with its hermitian conjugate.                                \n\
                                                                           \n\
      Z       (diagonal) operator required for the calculation of the      \n\
              fermionic hopping matrix elements (QSpace). It is defined by \n\
              ferminonic signs (-1)^ns required in the context as described\n\
              for FC above, where ns is the number of fermionic particles  \n\
              for the local states s(igma).                                \n\
                                                                           \n\
      zflag   variation of how to apply the fermionic signs Z on the       \n\
              hopping terms in the Hamilonian determined by the operators  \n\
              FC (default: 1 e.g. as described for FC above in the context \n\
              of all-abelian quantum numbers). Valid settings:             \n\
                                                                           \n\
              0: ignore Z alltogether                                      \n\
              1: H = F1'*ZF2 + F1*(ZF2)' + hermitian conjugate (default)   \n\
              2: intended for particle-hole-symmetric cases only:          \n\
                 H = (ZF1)'*(ZF2) at every second iteration only, WITHOUT  \n\
                 adding a further hermitian conjugate (already included!)  \n\
                 starting with NO Z-operators applied to A0.local (i.e.    \n\
                 same as for zflag==1), therefore starting with F1'*F1     \n\
                 at first iteration (i.e. using F1 to update operators,    \n\
                 and again F1 when building enlarged Hamiltonian.          \n\
                 For every other iteration, H = (Z*F1)'*(Z*F2) is applied  \n\
                 using Z with BOTH operators, consistent with particle-    \n\
                 hole symmetry!                                            \n\
              3: same as zflag=2, but Z is applied at every other iteration\n\
                 already starting with the very first iteration.           \n\
                                                                           \n\
      gg      strength of local operators FL along the chain similar to ff \n\
              i.e. not rescaled! (double, (L-1) x Nl, not including H0)    \n\
      FL      set of Nl local operators along the chain (QSpaces [Nl])     \n\
                                                                           \n\
   NB! FC and FL are considered the same throughout the chain, yet there   \n\
       strength is modified by ff and gg, respectively.                    \n\
                                                                           \n\
   NB! if the couplings ff/gg are the same for all FC and FL               \n\
       they may be represented by a single vector only.                    \n\
                                                                           \n\
   NB! The length L of the Wilson chain is deterimend by the length of ff, \n\
       i.e. L=length(ff)+1                                                 \n\
                                                                           \n\
   Options / Flags                                                         \n\
                                                                           \n\
    'Nkeep',. max number of states to keep at every NRG iteration (512)    \n\
    'NKEEP',. vector that specifies Nkeep individually for first few sites \n\
                                                                           \n\
    'Etrunc', use given energy threshold rather than Nkeep to trunctate    \n\
              unless Nkeep is reached (<=0 uses Nkeep only; default: 0).   \n\
    'ETRUNC', vector that specifies Etrunc individually for first few sites\n\
    'ET1'     Etrunc for 1st with trunctation (default: 1.2*Etrunc)        \n\
              if ET1<0, then ET1=Etrunc is used.                           \n\
                                                                           \n\
    'Estop',. stop dynamically once stable fixed point is reached, i.e.    \n\
              energies EK(1:Nkept/4) for iterations k and k+2 are the same \n\
              to within given value for Estop (default: 0, i.e. not used). \n\
    '-Estop'  same as 'Estop',1E-4                                         \n\
                                                                           \n\
    'NEE',..  number of states to store for energy flow diagram (Nkeep*Ns) \n\
    'fout',.  store data for every site into output file <fout>_##.mat     \n\
              if specified, the first RETURN argument is skipped (./NRG/NRG).\n\
    'deps',.. restore degeneracy within eps (1E-12)                        \n\
    'db',..   require minimum distance towards discarded energies (-1)     \n\
              (db<0 searches for maximum distance in vicinity of Nkeep)    \n\
    'dmax'..  maximum number of extra states to include for db (Nkeep/10)  \n\
              (mostly a safeguard in case db has been set too large).      \n\
    '-IdF'    build local Id from F-ops, fully ignoring A0 (NB! by default,\n\
              sigma-space of A0 is considered a proper Wilson bath site!)  \n\
    'ionly'   info only (incl. Eflow) i.e. do not store iterative NRG data \n\
    '-v'      verbose mode                                                 \n\
                                                                           \n\
   Output                                                                  \n\
                                                                           \n\
    NRG_DATA  MPS data (only if 'fout' was not specified; otherwise this   \n\
              data is stored into specified file structure):               \n\
                                                                           \n\
     .A[KT]   states K(ept) / T(hrown, discarded) at current NRG iteration \n\
              according to NRG truncation scheme (rank-3 QSpaces in LRs order)\n\
     .H[KT]   diagonal Hamiltonians in the basis of AK/AT. Since these are \n\
              (approximate) eigenstates, their diagonal is stored only).   \n\
                                                                           \n\
    info      info structure containing the following fields:              \n\
                                                                           \n\
     .EE      (rescaled) energy data for energy flow diagram               \n\
     .E0      substracted ground state energy at energy scale of iteration \n\
     .QS      local state space determined from FC/FL operators            \n\
     .FCC     sum of tensor product of all coupling operator sets          \n\
     .Lambda  Lambda used in the calculation (same as input)               \n\
     .ff      coupling along Wilson chain (rescaled, as used)              \n\
     .gg      local energies along Wilson chain (rescaled, as used)        \n\
     .L       length of Wilson chain (length(ff)+1)                        \n\
     .TN      (temperature) scale associated with chain length N           \n\
     .EScale  standard NRG energy scale used along the Wilson chain        \n\
     .stamp   time stamp upon finishing the calculation                    \n\
                                                                           \n\
   Alternative usage: NRGWilsonQS('MATfilename')                           \n\
   NRGWilsonQS('fpara.mat' [, 'fout_tag'])                                 \n\
                                                                           \n\
      This allows to read all input parameters from given fpara.mat file   \n\
      (MatLab binary) which allows pure standalone run of NRG without      \n\
      requiring MatLab (except for the API for binary *.mat I/O).          \n\
                                                                           \n\
      This way all variables and options MUST be defined by their names    \n\
      shown above; moreover, this usage must include a specification       \n\
      for fout_tag as all output data will be written there.               \n\
      fout_tag is either specified as 2nd argument above, else it should   \n\
      be defined in terms of the variable 'fout' in the parameter file     \n\
      fpara.mat (default for fout_tag: './NRG/NRG').                       \n\
                                                                           \n\
   AW (C) Apr 2006 ; May 2010 ; Nov 2014                                   \n";

// ChangeLog
// Wb,Apr22,10: added/updated deps and db to markSet

// Wb,Apr13,13: keep A0 and H0 for iter=0 exactly as provided
// by input to avoid confusion when reading NRG_00.mat data.
// tags: KEEP_AH0


#define PROG_TAG "NRG"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "NRGWilsonQS"
#endif

   #define WB_SPARSE_CLOCK


#define LOAD_CGC_QSPACE
#include "wblib.h"

template <class TQ, class TD>
void nrgBuildH4_abelian(
   QSpace<TQ,TD> &H4,
   const QSpace<TQ,TD> &HK,
   const QSpace<TQ,TD> Es,
   const double *f,
   const wbvector< QSpace<TQ,TD> > &F1K,
   const wbvector< QSpace<TQ,TD> > &F2,
   const double *g,
   const wbvector< QSpace<TQ,TD> > &FZ
);

template <class TQ, class TD>
void nrgBuildH4_cg(
   QSpace<TQ,TD> &H4,
   QSpace<TQ,TD> &A4,
   const QSpace<TQ,TD> &HK,
   const double *f,
   const wbvector< QSpace<TQ,TD> > &F1K,
   const wbvector< QSpace<TQ,TD> > &F2,
   char addHC,
   const double *g,
   const wbvector< QSpace<TQ,TD> > &FZ
);

template <class TQ, class TD>
void updateFOps(
   wbvector< QSpace<TQ,TD> > &F12,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   const wbvector< QSpace<TQ,TD> > &FC,
   char lflag=0
);

template <class TQ, class TD>
void nrgDispIter(
   const unsigned iter,
   const unsigned itermax,
   const wbvector<TD> &E4, const TD &EK,
   const QSpace<TQ,TD> &HK,
   const char disp
);

void checkGSDeg(const wbvector<double> &E4);

double getPhysE0(const wbvector<double> &DE);

   static void CleanUp(void) { Wb::CleanUp(FL); };

   double Lambda=1;

   double FNfac=-1;


inline double nrgScale(int iter) {


   if (FNfac<=0) wblog(FL,"ERR %s() invalid FNfac=%g",FCT,FNfac);

   double x = pow(Lambda, -double(iter)/2.) * 0.5*(Lambda+1) * FNfac;

   return (iter ? x : 1.0);
}


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,j,k,l,r,m=0,n,N, iter=0, Nkeep=256;
    double dbl, Etrunc1, Etrunc=0, Estop=0, deps=1E-12, db=-1;
    char disp=1, nostore=0, toFile=0, wf=-1, checkIdA=1; 
    char akflag=0, zflag=1, cgflag=0;
    unsigned gotNK=0, gotTR=0;
    int dmax=-1, NEE=-1;

    wbvector<double> E4,E0,ET,dd;
    wbMatrix<unsigned> D4;
    wbMatrix<double> ff,gg,EE,EK;
    wbMatrix<gTQ> QS;
    wbMatrix<int> NK;


    int idAK,idAT,idHK,idHT, idE0, idES;

    wbvector< QSpace<gTQ,gTD> > FC, FZ, F1, F2, F1K, HKALL;
    QSpace<gTQ,gTD> AK, AT, H4,A4, FCC, FX, Z;
    QSpace<gTQ,double> HK, HT, H0, A0, ID;

    const wbperm PA;

    OPTS opts;
    wbstring fout;
    mxArray *S, *a;

    const char* vpass[] = { "param", "Gamma" };
    unsigned npass=sizeof(vpass)/sizeof(const char*);
    mxArray* apass[npass];

    time_t tstart=time(NULL);

    wbSigHandler SIG(FL);
    WbClock nrgTime_1("build H4");
    WbClock nrgTime_2("eig(H4)");
    WbClock nrgTime_3("get AK/AT");
    WbClock nrgTime_4("update ops");
    WbClock nrgTime_5("mat I/O");

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    
    FNfac=-1;

    #ifdef __WB_MEM_CHECK__
       Wb::MemCheck(FL,"start");
    #endif

#ifdef MATLAB_MEX_FILE
    mexAtExit(CleanUp);
#endif
    Wb::CleanUp(NULL,1);
    
    if (nargin<=2) {

       if (nargin<1)
       usage(FL,"Invalid number of I/O arguments.");

       if (!mxIsChar(argin[0])) wblog(FL,
       "ERR Expecting name of parameter file for call with one argument.");

       if (mxGetString(argin[0],str,128)) wblog(FL,
       "ERR Error reading string of arg #1");

       opts.init(str);"options" by reference to file

       opts.getOpt("H0",a,'!'); H0.init(FL,a);
       opts.getOpt("A0",a,'!'); A0.init(FL,a);
       opts.getOpt("Lambda", Lambda,'!');

       opts.getOpt("ff",a,'!'); ff.init(FL,a);
       opts.getOpt("FC",a,'!'); mxInitQSpaceVec(FL,a,FC);
       opts.getOpt("Z", a,'!'); if (!mxIsEmpty(a)) Z.init(FL,a);

       opts.getOpt("gg",a);  if (a) gg.init(a);
       opts.getOpt("FL",a);  if (a) mxInitQSpaceVec(FL,a,FZ);

       for (i=0; i<npass; ++i) {
          opts.getOpt(vpass[i],a);
          apass[i] = a ? mxDuplicateArray(a) : NULL;
       }
 
       if (nargin>1) {
          if (!mxIsChar(argin[1]) || mxGetString(argin[1],str,128))
          wblog(FL,"ERR reading string of arg #2 (fout)");
          fout=str;
       }
    }
    else {

       m=6; if (nargin<int(m) || nargout>2)
       usage(FL,"Invalid number of I/O arguments.");

       wbindex iq("1 5 6", 1);
       if (mxIsEmpty(argin[5])) iq.len=2;

       for (i=0; i<iq.len; ++i) {
          try { r=-1; mxIsQSpace(FL,argin[iq[i]],r,0,-2); }
          catch (...) {
             wblog(FL,"ERR invalid QSpace at arg #%d",iq[i]+1);
          }
       }

       H0.init(FL,argin[0]);

       i=1; if (!mxIsEmpty(argin[i])) { r=-1;
          if (!mxIsQSpace(FL,argin[i],r))
          wblog(FL,"ERR Invalid QSpace at arg #%d", i+1);
          A0.init(FL,argin[i]); 
       }

       if (mxGetNumber(argin[2], Lambda)) wberror(FL,str);
       ff.init(FL,argin[3]);

       mxInitQSpaceVec(FL,argin[4], FC);
       if (!mxIsEmpty(argin[5])) Z.init(FL,argin[5]);

       if (nargin>=int(m+2))
       if ((mxIsEmpty(argin[m]) && mxIsEmpty(argin[m+1])) ||
           mxIsQSpaceVec(0,0,argin[m+1])
       ){
          gg.init(argin[m]);
          mxInitQSpaceVec(FL,argin[m+1], FZ); m+=2;
       }

#ifdef MATLAB_MEX_FILE
       for (i=0; i<npass; ++i)
       apass[i]=mexGetVariable("caller",vpass[i]);
#else
       npass=0;
#endif

       opts.init(argin+m, nargin-m);
    }

    wblog(FL,"=== %s ====================================",PROG);

    if (opts.getOpt("-q")) disp= 0; else
    if (opts.getOpt(FL,"-v")) disp=10; else
    if (opts.getOpt(FL,"-V")) disp=20;
    opts.getOpt(FL,"deps",deps);
    opts.getOpt(FL,"db",  db  );
    opts.getOpt(FL,"dmax",dmax);

    opts.getOpt(FL,"Nkeep",Nkeep);
    opts.getOpt(FL,"Etrunc", Etrunc ); Etrunc1=1.2*Etrunc;
    opts.getOpt(FL,"ET1",Etrunc1);     if (Etrunc1<0) Etrunc1=Etrunc;

    opts.getOpt(FL,"zflag",zflag);
    opts.getOpt(FL,"NEE", NEE);

    if (opts.getOpt(FL,"-Estop","using Estop=1E-4")) Estop=1E-4;
    else opts.getOpt(FL,"Estop",Estop);

    if (fout.isEmpty()) {
       opts.getOpt(FL,"fout", fout);
    }

    opts.getOpt(FL,"NKEEP",NK);
    opts.getOpt(FL,"ETRUNC",ET);

    FNfac=-1;
    opts.getOpt(FL,"FNfac", FNfac);


    nostore=opts.getOpt(FL,"ionly");
    if (opts.getOpt(0,0,"-IdF")) checkIdA=0;

    opts.checkAnyLeft();
 
    #ifdef __WB_MEM_CHECK__
       Wb::MemCheck(FL,"info");
    #endif


    if (!H0.isHConj() || H0.rank(FL)!=2) wblog(FL,
       "ERR H0 must be Hermitian rank-2 object!\n%s",str);

    if (H0.itags.isEmpty()) {
       H0.init_itags(FL,"nrg:HK",0);
    }

    i=H0.skipZeroOffDiag(1E-14,'b'); if (i) {
       wblog(FL," *  skipped %d off-diag. blocks from H0",i);
#ifdef MATLAB_MEX_FILE
       H0.put("H0_");
#endif
    }


    AK=A0;
    if (AK.isEmpty()) { akflag=1;
       wblog(FL,"WRN building A0 (empty input) from H0");
       AK.initIdentity(H0);
       AK.PrependSingletons(3);
       AK.Permute("1 3 2");
    }
    else {
       if (!AK.isConsistentR(3)) wblog(FL,
          "ERR invalid %s (%s)", nargin>1 ? "arg #2" : "QSpace A0", str);

       AK.contract(FL,"1,3;*",AK,"1,3",FX);
       if (!FX.isIdentityMatrix()) {
          MXPut(FL,"a").add(AK,"AK").add(FX,"FX");
          wblog(FL,"ERR A0 does not describe orthonormal basis");
       }

       HT.init2DiffOp(AK,1,H0,0);
       if (!HT.isEmpty()) { m=HT.DATA.len;
          HT.Append2AndDestroy(FL,H0); sprintf(str," *  "
            "increased H0 by %d diagonal zero-block%s (%ld->%ld)",
             m, m!=1 ? "s":"", H0.DATA.len-m, H0.DATA.len);
          wblogBuf.push(FL,str);
       }

    }

    if (A0.itags.isEmpty()) {
       A0.init_itags(FL,"nrg:AK",0);
    }

    if (H0.qtype.noneAbelian() && zflag<=1) wblog(FL,
       "WRN %s() having zflag=%d with sym=%s\n"
       "(hint: particle-hole symmetry requires zflag>=2)",
       PROG, zflag, H0.qtype.toStr().data
    ); 


    if (!(gotNK=NK.numel())) { NK.init(1,1); NK[0]=Nkeep; }


  { QSpace<gTQ,gTD> UK,UT,HK_,HT_;

    H0.EigenSymmetric(
       UK, UT, HK_,HT_,E4, D4, NK[0]
    );

    if (!HT_.isEmpty()) wblog(FL,"WRN %s() "
       "already truncate local Hamiltonian H0 (%g) !??",FCT,NK[0]);

    HK=H0; AK=A0;
  }




    for (k=0; k<FC.len; ++k) {
       cgflag=FC[k].gotCGS(FL); r=FC[k].rank(FL);
       if ((cgflag<=0 && r>2) || (cgflag>0 && r>3)) wblog(FL,
          "ERR %s() got rank-%d operator (%d; '%s')",PROG,r,cgflag,
          FC[k].qtype.toStr().data
       );
    }

    if (zflag) {
       if (Z.isEmpty()) wblog(FL,
          "ERR %s() got empty Z-operator for fermionic signs",PROG);
       if (!Z.isDiagMatrix()) wblog(FL,
          "ERR %s() got not diagonal operator Z",PROG);
    }


    if (zflag<0 || zflag>3) wblog(FL,
       "ERR %s() invalid zflag=%d",PROG,zflag);
    m=FC.len; F1.init(m); F2.init(m);

    for (i=0; i<m; ++i) {
       if (zflag<=1)
            FC[i].hconj(F1[i]);
       else F1[i]=FC[i];

       if (zflag)
          { Z.contract(FL,2,FC[i],1, F2[i]);
            F2[i].otype=FC[i].otype; }
       else F2[i]=FC[i];
    }

    if (1) { QSpace<gTQ,gTD> FY;

       if (FC.len && FC[0].rank(FL)==2) {
          for (i=0; i<F1.len; ++i) {
             F1[i].TensorProd(F2[i],FX).hconj(FY);
             FX.Append2AndDestroy(FL,FCC);
             FY.Append2AndDestroy(FL,FCC);
          }
       }
       else if (FC.len && FC[0].otype==QS_OPERATOR) {
          wbperm P("1,3,2,4");

          for (i=0; i<F1.len; ++i) {
             if (zflag<=1) {
                F1[i].contract(FL,"3", F2[i],"3",FX,P);
                FX.hconj(FY);
             }
             else {
                F1[i].contract(FL,"3*",F1[i],"3",FX,P);
                F2[i].contract(FL,"3*",F2[i],"3",FY,P);
             }
             FX.Append2AndDestroy(FL,FCC);
             FY.Append2AndDestroy(FL,FCC);
          }
       }
       else wblog(FL,"ERR invalid FC operator set");

       FCC.makeUnique(); FCC.SkipZeroData(1E-14);

       if (!FCC.isHConj()) {
          MXPut(FL,"I4").add(F1,"F1").add(F2,"F2").add(FCC,"FCC")
             .add(zflag,"zflag");
          wblog(FL,"ERR %s() coupling "
             "F'*F does not yield hermitian H\n'%s'",PROG,str
          );
       }
    }

    if (ff.isEmpty() || FC.isEmpty()) wblog(FL,
       "ERR %s() got empty couplings (%dx%d,%d)",
        PROG,ff.dim1,ff.dim2,FC.len);
    if (!ff.isNormal()) wblog(FL,
       "ERR %s() got NaN in couplings ff !?",PROG);


    if (F1.len>1 && ff.isVector()) {
       wbvector<gTD> x(ff.numel(), ff.data);
       ff.init(x.len,F1.len);
       for (i=0; i<ff.dim1; ++i)
       for (j=0; j<ff.dim2; ++j) ff(i,j)=x[i];
    }
    else if (ff.isVector() && ff.dim1==1) {
       SWAP(ff.dim1,ff.dim2);
    }

    if (!gg.isEmpty() && gg.isVector()) {
       if (FZ.len>1) {
          wbvector<gTD> x(gg.numel(), gg.data);
          gg.init(x.len,FZ.len);
          for (i=0; i<gg.dim1; ++i)
          for (j=0; j<gg.dim2; ++j) gg(i,j)=x[i];
       }
       else if (gg.isVector() && gg.dim1==1)
       SWAP(gg.dim1,gg.dim2);
    }

    N=1+ff.dim1;

    if (!gg.isEmpty() && ff.dim1!=gg.dim1) wblog(FL,
       "ERR chain length mismatch between ff and gg (%d,%d)",
        ff.dim1, gg.dim1);
    if (ff.dim2!=F1.len) wblog(FL,
       "ERR dimension of ff inconsistent with number of F1 ops (%d,%d)",
        ff.dim2, F1.len);
    if (gg.dim2!=FZ.len) wblog(FL,
       "ERR dimension of gg inconsistent with number of FL ops (%d,%d)",
        gg.dim2, FZ.len);
    if (!gg.isNormal()) wblog(FL,
       "ERR %s() got NaN in local energies gg !??",PROG);

    if (F1.len!=F2.len) wblog(FL,"ERR %d/%d",F1.len,F2.len);
    if (!Z.isEmpty() && (Z.QDIM!=HK.QDIM || Z.rank(FL)!=2)) wblog(FL,
       "ERR severe rank mismatch of Z (%dx%d/%d)",
        Z.QIDX.dim1, Z.QIDX.dim2, Z.QDIM
    );

    for (k=0; k<F1.len; ++k) {
       if (F1[k].SkipZeroData()) { wblog(FL,
          "WRN coupling operator set #%d is zero (skipped)", k+1);
           continue;
       }
       if (F1[k].isEmpty()) { wblog(FL,
          "WRN no contribution by coupling operator set #%d !??", k+1);
           for (i=0; i<ff.dim1; ++i) ff(i,k)=0;
           continue;
       }
       if (F1[k].QDIM!=HK.QDIM || !F1[k].isOperator()) wblog(FL,
          "ERR severe (rank) mismatch of F1 (q=%d,%d; r=%d/%s/%d)",
           F1[k].QDIM, HK.QDIM, F1[k].rank(FL), QS_STR[F1[k].otype],
           F1[k].isOperator()
       );
       if (F2[k].QDIM!=HK.QDIM || !F2[k].isOperator()) wblog(FL,
          "ERR severe (rank) mismatch of F2 (q=%d,%d; r=%d/%s/%d)",
           F2[k].QDIM, HK.QDIM, F2[k].rank(FL), QS_STR[F2[k].otype],
           F2[k].isOperator()
       );
       F1[k].checkQ(FL,F2[k]);
       if (F1[k].QIDX.dim2!=F2[k].QIDX.dim2 || 
           F1[k].otype!=F2[k].otype) wblog(FL,
          "ERR severe F[12] operator rank mismatch (%d: %d/%d; %d/%d; %s/%s)",
           k+1, F1[k].QIDX.dim2, F2[k].QIDX.dim2, F1[k].QDIM, F2[k].QDIM,
           QS_STR[F1[k].otype], QS_STR[F2[k].otype]
       );

       if (k) {
          char ic=F1[k].gotCGS(FL); if (cgflag!=ic)
          wblog(FL,"ERR %s() CG inconsistency (%d/%d)",PROG,ic,cgflag);
       }
       else cgflag=F1[k].gotCGS(FL);
    }

    for (k=0; k<FZ.len; ++k) {
       if (FZ[k].SkipZeroData()) { wblog(FL,
          "WRN local operator #%d contains zero data - skip.", k+1);
           continue;
       }
       if (FZ[k].isEmpty()) { wblog(FL,
          "WRN no contribution by local operator #%d !??", k+1);
           for (i=0; i<gg.dim1; ++i) gg(i,k)=0;
           continue;
       }
       if (FZ[k].QDIM!=HK.QDIM || FZ[k].rank(FL)!=2) { wblog(FL,
          "ERR severe local operator inconsistency (%d,%d)",
           FZ[k].QDIM, FZ[k].rank(FL));
       }
       if (!FZ[k].isHConj()) wblog(FL,
       "ERR FL[%d] not a Hermitian rank-2 object", k+1);

       if (k) {
          char ic=FZ[k].gotCGS(FL); if (cgflag!=ic)
          wblog(FL,"ERR %s() CG inconsistency (%d/%d)",PROG,ic,cgflag);
       }
       else cgflag=FZ[k].gotCGS();

       if (FZ[k].rank(FL)!=2) wblog(FL,
          "ERR %s() got non-scalar rank-%d local operator",PROG,
          FZ[k].rank(FL)
       );
    }

    if (Nkeep>1024 && fout.isEmpty()) {
       fout="./NRG/NRG"; wblog(FL,
       "WRN saving data to file for Nkeep=%d (>1024)", Nkeep);
    }
    toFile = !fout.isEmpty();

    if (nargout==0 && !toFile) { printf("\n>> %s:%d\n"
       ">> nobody wants anything - done :)\n\n",
       FL); return;
    }

    if (!Nkeep) {
       wblog(FL,"<i> %NNkeep=%d life is simple - exit B)", Nkeep);
       return;
    }

    if (Estop<0) wblog(FL,"ERR invalid Estop=%g",Estop);
    if (Estop>0.1) wblog(FL,"ERR invalid Estop=%g (expected <0.1)",Estop);


    if (FNfac<0) {
       FNfac=1;
       FNfac=ff(ff.dim1-1,0)/nrgScale(ff.dim1);
       if (FNfac<0.1 || FNfac>10) {
       wblog(FL,"WRN %s() got FNfac=%.3g!??",PROG,FNfac); }
    }

    wblog(FL," *  %48R","-");
    l=sprintf(str,"Lambda=%g, L=%d",Lambda,N);
    if (Etrunc>0)
         { l+=sprintf(str+l,", Etrunc=%.3g (@%d)",Etrunc,Nkeep); }
    else { l+=sprintf(str+l,", Nkeep=%d",Nkeep); }
    if (!toFile)
         { l+=sprintf(str+l,", internal"); }
    wblog(FL," *  %s\n *  sym=%s",str,A0.qtype.toStr('V').data);

    str[0]=0; l=0;
    if (disp!=1) l+=sprintf(str+l," disp=%d",disp);
    if (nostore) l+=sprintf(str+l," %s",toFile? "NOSTORE":"noStore");
    if (l) wblog(FL," *  flags:%s",str);

    wblog(FL,"=== ================================================");

    if (!gg.isEmpty())
       wblog(FL," *  got %d local operator%s",FZ.len,FZ.len!=1 ? "s":"");
    wblogBuf.flush();

    if (toFile) {
       wbstring cwd(128); if (getcwd(cwd.data, 127)==NULL) wblog(FL,
       "ERR cwd length exceeds maximum length %d\n%s",cwd.len,cwd.data);

       if (nostore) wblog(FL,"NB! saving info only");
       else {
          for (iter=0; iter<1024; ++iter) {
              sprintf(str, "%s_%02d.mat", fout.data, iter);
              if (remove(str)) break;
          }   if (iter>1000) wblog(FL,"WRN iter=%d ???", iter);
       }

       char *s1=str, *s2=str+32; s2[0]=0;

       if (nostore)
            strcpy (s1,"info");
       else sprintf(s1,"{%02d}",N);

       if (iter) sprintf(s2," (%d files removed)",iter);

       if (disp) printf("\n");
       printf("   pwd: %s\n"
              "   out: %s_%s.mat%s\n",
       repHome(cwd).data, repHome(fout).data, s1, s2);
       if (disp) printf("\n");
    }
    else {
       wblog(FL,"TST pwd: out: toFile=%d (0x%lx: %s)",
       toFile, fout.data, fout.data ? fout.data:"");
    }


    { wbvector< wbMatrix<gTQ> > qq(F1.len+FZ.len);
      wbvector<INDEX_T> sz;

      for (k=i=0; i<F1.len; ++i) F1[i].getQDim(qq[k++],sz);
      for (  i=0; i<FZ.len; ++i) FZ[i].getQDim(qq[k++],sz);

      QS.CAT(1,qq).makeUnique();

      if (checkIdA) { wbMatrix<gTQ> Qs;
         AK.getQsub(2,Qs);
         Qs.makeUnique(); if (Qs!=QS) {
#ifdef MATLAB_MEX_FILE
         Qs.Print("A0->QS"); QS.Print("FX->QS");
#endif
         if (akflag) wblog(FL,
            "ERR local QIDX inconsistency%N%N    Hint: "
            "A0 not specified and failed to construct default from H0");
         else wblog(FL,
            "ERR local QIDX inconsistency in specified A0%N%N    Hint: "
            "ensure that H0 contains complete basis even if data is zero"
         );
      }}
    }

    ID.QIDX.Cat(2,QS,QS); ID.QDIM=QS.dim2; ID.setupDATA();

    for (k=0; k<QS.dim1; ++k) {
       gTQ *qk=ID.QIDX.rec(k);

       for (i=0; i<F1.len; ++i) { if (F1[i].findDimQ(qk,n)) break; }
       if (i==F1.len) {
          for (i=0; i<FZ.len; ++i) { if (FZ[i].findDimQ(qk,n)) break; }
          if (i==FZ.len) {
             MXPut(FL).add(ID,"ID").add(F1,"F1").add(FZ,"FZ");
             wblog(FL,"ERR failed to find Id.Q(%d,:) in F[1L]",k+1);
          }
       }
       ID.DATA[k]->initIdentity(n);
    }

    ID.qtype=AK.qtype;

    if (cgflag>0 || !H0.CGR.isEmpty())
         ID.initIdentityCGS();
    else ID.itags.init_qdir("+-");

    if (checkIdA) {
       AK.contract(FL,"1,2;*",AK,"1,2",FX);
       if (!FX.DATA.len) wblog(FL,"ERR AK/AK contracted to empty!??");
       FX*=(1/FX.DATA[0]->data[0]);

       try { dbl=FX.normDiff2(ID); }
       catch (...) {
          MXPut(FL).add(FX,"FX").add(ID,"ID");
          throw;
       }
       if (dbl>1E-16) { 
          MXPut(FL,"i").add(ID,"ID").add(FX,"id");
          wblog(FL,"ERR Id inconsistency (%.4g)",dbl);
       }
    }


    dbl=1/nrgScale(0);

    if (dbl!=1) wblog(FL,"ERR %s() got nrgScale(0)=%g !??",PROG,dbl);

    dbl=pow(Lambda,-N/2.);
    if (ff.dim1 && (fabs(ff.recMax(ff.dim1-1) / ff.recMax(0)) > 10*dbl))
    wblog(FL,
       "WRN ff must fall off exponentially\n[%.3g .. %.3g; %.3g] !??",
        ff.recMax(N-1), ff.recMax(0), dbl);

    for (i=0; i<ff.dim1; ++i) { dbl=1/nrgScale(i+1);
       for (j=0; j<ff.dim2; ++j) ff(i,j)*=dbl;
       for (j=0; j<gg.dim2; ++j) gg(i,j)*=dbl;
    }

    S=mxCreateStructMatrix(1, toFile || nostore ? 1 : N, 0, NULL);

    idAK = mxAddField(FL, S,"AK");
    idAT = mxAddField(FL, S,"AT");
    idHK = mxAddField(FL, S,"HK");
    idHT = mxAddField(FL, S,"HT");

    idE0 = mxAddField(FL, S,"E0");
    idES = mxAddField(FL, S,"ES");

    if (idAK<0 || idAT<0 || idHK<0 || idHT<0 || idE0<0 || idES<0)
    wblog(FL,
      "ERR failed to add fields to structure (%d,%d,%d,%d; %d,%d)",
       idAK,idAT,idHK,idHT, idE0,idES
    );

    if (!NK.isVector()) wblog(FL,"ERR invalid NKEEP");
    else {
       wbvector<int> nk(NK.numel(),NK.data);
       NK.init(N, cgflag>0 ? 4:2);
       for (n=MIN(unsigned(nk.len),N), i=0; i<n; ++i) {
          NK.setRec(i,nk[i]);
       }
       for (; i<N; ++i) NK.setRec(i,Nkeep);
       NK(N-1,0)=0;
    }

    if (NEE<0) NEE=2*Nkeep;

    EE.init(N, (unsigned)NEE); EE.set(NAN);
    E0.init(N); EK.init(N,3); EK.set(NAN);
    HKALL.init(N);


    #ifdef __WB_MEM_CHECK__
       Wb::MemCheck(FL,"info");
    #endif

    for (iter=0; iter<N; ++iter) { SIG.call99();

       if (iter) { nrgTime_1.resume();
          HK *= (nrgScale(iter-1) / nrgScale(iter));

          if (cgflag<=0) { try {
             nrgBuildH4_abelian(H4, HK, ID,
                ff.rec(iter-1), F1K, F2,
                gg.dim1 ? gg.rec(iter-1) : NULL, FZ);
             }
             catch (...) {
                wblog(FL,"ERR %s() iter=%d",PROG,iter); 
             }

          }
          else { try {
             nrgBuildH4_cg(H4,A4, HK,
                ff.rec(iter-1), F1K,
                (zflag<=1 || !wf) ? F2 : F1,
                zflag<=1 ? 1 : 0,
                gg.dim1 ? gg.rec(iter-1) : NULL, FZ);
             }
             catch (...) {
                MXPut(FL,"q").add(F1K,"F1K").add(F1,"F1").add(F2,"F2")
               .add(FZ,"FZ"); wblog(FL,"ERR iter=%d",iter); 
             }

             if (iter==1) { if (!H4.isHConj(0,0,1E-12,'v')) {
                MXPut(FL,"a").add(HK,"HK").add(ff.getRec(iter-1),"ff")
                  .add(gg.dim1? gg.getRec(iter-1) : wbvector<double>(),"gg")
                  .add(F1K,"F1").add((zflag<=1 || !wf) ? F2 : F1,"F2")
                  .add(FZ,"FZ");
                wblog(FL,
                  "ERR H4 got non-hermitian operator setting !??\n"
                  "hint: got correct zflag=%d ?\n%s",zflag,str);
                }
             }
          }

          nrgTime_1.stop();
          nrgTime_2.resume();

          #ifdef __WB_MEM_CHECK__
             Wb::MemCheck(FL,"info");
          #endif

          H4.EigenSymmetric(
             AK, AT, HK, HT, E4, D4, NK(iter,0),
             iter<ET.len && ET[iter]>0 ?  ET[iter]
             : ((Etrunc1>0 && !gotTR) ? Etrunc1 : Etrunc),
             E0.data+iter, PA, deps, db, dmax
          );

          if (!gotTR) {
             if (NK(iter,0)<(int)E4.len) gotTR=iter+1; else {
                double x=1.1*(iter+1>=ET.len ? Etrunc : ET[iter+1]);
                if (E4.end()>x && (ET.len || Etrunc1>x)) { gotTR=iter+2; }
             }
          }

          #ifdef __WB_MEM_CHECK__
             Wb::MemCheck(FL,"info");
          #endif

          nrgTime_2.stop();
          nrgTime_3.resume();

          if (cgflag>0) { INDEX_T nx;
             n=AK.getDim(AK.rank(FL)-1, &nx); NK(iter,1)=nx;
             NK(iter,2)=D4.colSum(0);
             NK(iter,3)=D4.colSum(1);

             A4.contract(2,AK,1,AK,"1,3,2"); if (!AT.isEmpty())
             A4.contract(2,AT,1,AT,"1,3,2");
          }
          else {
             AK.Permute("2,3,1"); if (!AT.isEmpty())
             AT.Permute("2,3,1");
             NK(iter,1)=D4.colSum(0);
          }

          nrgTime_3.stop();
       }
       
       AK.init_itags(FL,"nrg:AK",iter);
       AT.init_itags(FL,"nrg:AT",iter);
       HK.init_itags(FL,"nrg:HK",iter);
       HT.init_itags(FL,"nrg:HT",iter);

       HKALL[iter]=HK;

       wf=(char(iter%2)==(zflag-2));

       nrgTime_4.resume();


       updateFOps(F1K, AK,AK, 
          (zflag<=1 || wf) ? F1 : F2,
          iter+1>=N
       );

       nrgTime_4.stop();

       memcpy(
          EE.rec(iter), E4.data,
          MIN(size_t(NK(iter,0)),EE.dim2)*sizeof(double)
       );

       if (iter>18 && Estop>0) {
          unsigned nk=MIN(MIN(int(EE.dim2),NK(iter-2,0))/4,NK(iter,0)/4);
          double x=sqrt(rangeNormDiff2(EE.rec(iter-2), EE.rec(iter), nk));
          if (x<Estop && iter+2<N) {
             wblog(FL,"%N==> "
              "NB! energy flow converged at %.3g (N=%d->%d)",x,N,iter+2);
             N=iter+2;
             Estop=-Estop;

             EE.dim1=N; EK.dim1=N; NK.dim1=N;
             E0.len=N; HKALL.len=N;
             ff.dim1=N; if (gg.dim1>N) { gg.dim1=N; }
          }
       }

       n=NK(iter,0);
       EK(iter,0)=(n && n<E4.len ? E4[n-1] : E4.end());
       EK(iter,1)=(     n<E4.len ? E4[n  ] : NAN);
       EK(iter,2)=E4.end();

       dbl = !ISNAN(EK(iter,1)) ? EK(iter,1) :  EK(iter,0);
       nrgDispIter(iter,N,E4, dbl, HK, disp ?
         (disp>1 ? disp : (Etrunc<=0 ? 0 : (dbl>=0.8*Etrunc ? -1 : +1)))
          : -1
       );


       if (!nostore) {
          k = toFile ? 0 : iter;

          mxReplaceField(FL, S, k, idAK, AK.toMx());
          mxReplaceField(FL, S, k, idAT, AT.toMx());
          mxReplaceField(FL, S, k, idHK, HK.toMx());
          mxReplaceField(FL, S, k, idHT, HT.toMx());

          mxReplaceField(FL, S, k, idE0, numtoMx(E0[iter]));
          mxReplaceField(FL, S, k, idES, numtoMx(nrgScale(iter)));


          if (toFile) {
             Wb::UseClock NT5(&nrgTime_5);

             sprintf(str, "%s_%02d.mat", fout.data, iter);

             if (disp>10) wblog(FL,
                "I/O saving data to file `%s'", basename(str));

             Wb::matFile F(FL,str,"w");
             for (n=mxGetNumberOfFields(S), i=0; i<n; ++i) {
                F.put(FL,
                   mxGetFieldNameByNumber(S,i),
                   mxGetFieldByNumber(S,0,i)
                );



             }
             F.close();
          }
       }
       doflush();

       #ifdef __WB_MEM_CHECK__
          Wb::MemCheck(FL,"info");
          wblog(FL,"--- %50R","-");
       #endif
    }

    mxArray *aC=gCG.toMx(), *aR=gRG.toMx(), *aM=gCX.toMx();
    mexPutVariable("caller","gCG",aC);
    mexPutVariable("caller","gRG",aR);
    mexPutVariable("caller","gCX_nrg",aM);

    if (disp) printf("\n");


    for (dbl=1E99, i=0; i<EK.dim1; ++i) {
       if (!ISNAN(EK(i,1)) && EK(i,1)>EK(i,0)) {
          if (dbl>EK(i,1)) dbl=EK(i,1); 
       }
    }
    if (Etrunc) str[0]=0;
    else sprintf(str," (Etr=%.6g)",dbl);

    sprintf(str+128, (Etrunc>0 && 0.9*Etrunc<dbl) ? " * ":"WRN");
    m=NK.colMax(0,i);

    if (NK.dim2==4) {
       n=NK.colMax(2,j); wblog(FL,
       "%s NK=%d/%d (%d/%d)%s",str+128,m,n,NK(i,1),NK(j,3),str);
    }
    else {
       n=NK.colMax(1,j); wblog(FL,
       "%s NK=%d/%d%s",str+128,m,n,str);
    }

    if (Etrunc) wblog(FL,"%s Etr=%.4g [%g]",str+128,dbl,Etrunc);


    if ((!toFile && nargout>0) || (toFile && nargout>1))
    argout[0] = S;

    sprintf(str,"NRG data obtained using %s",PROG);
    wbstring ver(str);


    MXPut I(0,0,"Inrg");
    I.add(ver,"istr")
     .add(wbTimeStamp(),"stamp");
    I.addP(MXPut(0,0)
         .add (wbstring().time_sys(tstart),"started")
         .add (wbstring().time_sys(),"finished")
         .addP(nrgTime_1.toMx(),"build_H4")
         .addP(nrgTime_2.toMx(),"eig_H4")
         .addP(nrgTime_3.toMx(),"get_AKT")
         .addP(nrgTime_4.toMx(),"update_ops")
         .addP(nrgTime_5.toMx(),"io")
     .toMx(),"usage");

#ifdef __WB_MEM_CHECK__
    I.addP(Wb::gML.totStr('l'),"MEM");
#endif

    MXPut IS(0,0);
      IS.add(ID,"ID").add(H0,"H0").add(A0,"A0").add(FC,"FC").add(Z,"Z")
        .add(ff,"ff").add(F1,"F1").add(F2,"F2").add(gg,"gg");
      if (!FZ.isEmpty()) IS.add(FZ,"FL"); else IS.add(dd.init(),"FL");
      IS.add(FCC,"FCC").add(zflag,"zflag");

    dd.init(N); for (i=0; i<N; ++i) dd[i]=nrgScale(i);

    I.addP(EE.toMx('r'),"EE").add(HKALL,"HK").add(E0,"E0").add(EK,"EK")
     .add(dd,"EScale").add(getPhysE0(E0),"phE0")
     .addP(IS.toMx(),"ops");

    dd.init(4);
      dd[0]=CG_EPS1; dd[2]=CG_SKIP_EPS1;
      dd[1]=CG_EPS2; dd[3]=CG_SKIP_EPS2;
    I.add(Lambda,"Lambda")
     .addP(MXPut(0,0)
        .add(dd,"cg_eps").add(deps,"deps").add(db,"db")
        .add(dmax,"dmax").add(-Estop,"Estop").add(disp,"disp")
        .add(toFile,"toFile").addP(A0.qtype.toStr().data,"sym")
     .toMx(),"paras");

    for (i=0; i<npass; ++i) if (apass[i]) {
       I.addP(apass[i],vpass[i]); apass[i]=NULL;
    }

    I.addP( MXPut(0,0)
       .add(Etrunc,"Etrunc").add(Etrunc1,"Etrunc1")
       .add(gotTR,"itrunc").add(ET,"ETRUNC").add(FNfac,"FNfac")
     .toMx(),"Itr")
     .add(Nkeep,"Nkeep").add(NK,"NK").add(N,"N");

    I.save2(S);


    if (toFile) {
       sprintf(str, "%s_info.mat", fout.data);

       if (disp>1) wblog(FL,
          "I/O saving info data to `%s'", basename(str));

       Wb::matFile F(FL,str,"w");

       for (n=mxGetNumberOfFields(S), i=0; i<n; ++i)
       matPutVariable(F.mfp,
          mxGetFieldNameByNumber(S,i),
          mxGetFieldByNumber(S,0,i)
       );

       matPutVariable(F.mfp,"gCG",aC);
       matPutVariable(F.mfp,"gRG",aR);
       matPutVariable(F.mfp,"gCX",aM);
    }

    mxDestroyArray(aC);
    mxDestroyArray(aR);
    mxDestroyArray(aM);

    if (nargout>=2)        argout[1]=S; else
    if (toFile && nargout) argout[0]=S;
    else mxDestroyArray(S);


    if (!disp) {
       nrgTime_1.reset();
       nrgTime_2.reset();
       nrgTime_3.reset();
       nrgTime_4.reset();
    }

    if (disp) { wblog(FL,"FIN %s() I/O time usage: %s",PROG,
       nrgTime_5.gettime("asm").data);
       CleanUp();
    }
    nrgTime_5.reset();




#ifdef __WB_MEM_CHECK__

    wblog(FL,"%NMTR %s() check memory management (before cleanup)",PROG);

    Wb::MemCheck(FL,"info");

    opts.init();

    ff.init(); gg.init(); dd.init();
    E4.init(); EE.init(); E0.init(); EK.init();
    D4.init(); QS.init(); NK.init(); fout.init();

    FC.init(); FZ.init(); FCC.init(); FX.init();
    F1.init(); F2.init(); F1K.init(); Z.init();

    AK.init(); AT.init(); A0.init(); A4.init(); ID.init();
    HK.init(); HT.init(); H0.init(); H4.init(); HKALL.init();


    gCG.BUF.clear(); gCG.buf3.clear();
    gRG.buf.clear(); gCX.clear();

    Wb::MemCheck(FL,"info");

    wblog(FL,"MTR %s() check memory management (after cleanup)%N",PROG);





#else
#endif


};


template <class TQ, class TD>
void nrgBuildH4_abelian(
    QSpace<TQ,TD> &H4,
    const QSpace<TQ,TD> &HK,
    const QSpace<TQ,TD> Es,
    const double *f,
    const wbvector< QSpace<TQ,TD> > &F1K,
    const wbvector< QSpace<TQ,TD> > &F2,
    const double *g,
    const wbvector< QSpace<TQ,TD> > &FZ
){
    unsigned i; char xflag=0;

    QSpace<TQ,TD> HX,HY, EL, FX;
    wbMatrix<TQ> QQ;

    EL.initIdentity(HK,'r');


if (!H4.isHConj()) wblog(FL,"XXX %s() H4 not h.conj. !??",FCT); 

    Es.TensorProd(HK,HX).save2(H4);

    if (HK.QIDX.dim1==1 || FZ.len) {
       H4.ExpandDiagonal(2,4); xflag=1;
    }

if (!xflag) { H4.ExpandDiagonal(2,4); xflag=1; }
if (!H4.isHConj()) {
 MXPut(FL,"q0").add(H4,"H4").add(EL,"EL").add(HK,"HK").add(Es,"Es");
 wblog(FL,"ERR %s() H4 not h.conj. !??",FCT); 
}

    for (i=0; i<FZ.len; ++i) { if (g[i]==0) continue;
       FX=FZ[i]; FX*=g[i];
       FX.TensorProd(EL,HX).Append2AndDestroy(FL,H4);
if (!H4.isHConj()) wblog(FL,"ERR %s() H4 not h.conj. (%d) !??",FCT,i);
    }


if (!H4.isHConj()) wblog(FL,"XXX %s() H4 not h.conj. !??",FCT);

    for (i=0; i<F1K.len; ++i) { if (f[i]!=0) {
       FX=F2[i]; FX*=f[i];
       FX.TensorProd(F1K[i],HX); HX.hconj(HY);
       HX.Append2AndDestroy(FL,H4);
       HY.Append2AndDestroy(FL,H4);
if (!H4.isHConj()) wblog(FL,"XXX %s() H4 not hconj (%d) !?",FCT,i);
    }}

    H4.makeUnique();
if (!H4.isHConj()) wblog(FL,"XXX %s() H4 not h.conj. !??",FCT);

    if (!xflag)
    H4.ExpandDiagonal(2,4);

if (!H4.isHConj()) {
 MXPut(FL,"q").add(H4,"H4").add(EL,"EL").add(HX,"HX").add(Es,"Es");
 wblog(FL,"ERR %s() H4 not h.conj. !??",FCT); 
}

};


template <class TQ, class TD>
void nrgBuildH4_cg(
   QSpace<TQ,TD> &H4,
   QSpace<TQ,TD> &A4,
   const QSpace<TQ,TD> &HK,
   const double *f,
   const wbvector< QSpace<TQ,TD> > &F1K,
   const wbvector< QSpace<TQ,TD> > &F2,
   char addHC,
   const double *g,
   const wbvector< QSpace<TQ,TD> > &FZ
){
   unsigned i;
   QSpace<TQ,TD> Q,FX,HX,AX;

   for (i=0; i<F1K.len; ++i) if (F1K[i].isEmpty()) wblog(FL,
       "WRN got empty space F1K[%d]",i+1);
   for (i=0; i<F2.len; ++i) if (F2[i].isEmpty()) wblog(FL,
       "WRN got empty space F2[%d]",i+1);
   for (i=0; i<FZ.len; ++i) if (FZ[i].isEmpty()) wblog(FL,
       "WRN got empty space FZ[%d]",i+1);

   A4.initIdentityCG(F1K,F2,"1 3 2");


   Q=HK; Q.ExpandDiagonal();
   Q.contract(2,A4,1,AX);

   #ifdef __WB_MEM_CHECK__
      Wb::gML.printSize(FL,'l');
   #endif

   A4.contract("1,3;*",AX,"1,3",H4);

   #ifdef __WB_MEM_CHECK__
      Wb::gML.printSize(FL,'l');
   #endif

   for (i=0; i<F1K.len; ++i) { if (f[i]==0) continue;
      F2[i].times(f[i],Q); A4.contract(FL,3,Q,2,FX);

      if (F2[i].rank(FL)==2) {
         if (!addHC) wblog(FL,"WRN got addHC=%d !?",addHC);
         F1K[i].contract(FL,2,FX,1,AX);
      }
      else {
         if (addHC)
              F1K[i].contract(FL,"2,3",  FX,"1,4",AX);
         else F1K[i].contract(FL,"1,3;*",FX,"1,4",AX);
      }

      #ifdef __WB_MEM_CHECK__
         Wb::gML.printSize(FL,'l');
      #endif

      try {
         A4.contract(FL,"1,3;*",AX,"1,3",HX); }
      catch (...) {
         MXPut(FL,"i4").add(A4,"A4").add(AX,"AX").add(HX,"HX")
          .add(FX,"FX").add(F1K[i],"F1K").add(i+1,"i"); 
         wblog(FL,"ERR %s()",FCT);
      }

      #ifdef __WB_MEM_CHECK__
         Wb::gML.printSize(FL,'l');
      #endif

      if (addHC) { QSpace<TQ,TD> HY;
      HX.hconj(HY).Append2AndDestroy(FL,H4); }
      HX.Append2AndDestroy(FL,H4);
   }

   AX.init();

   for (i=0; i<FZ.len; ++i) { if (g[i]==0) continue;
      FZ[i].times(g[i],FX);
      A4.contract(FL,"1,3",A4.contract(FL,3,FX,2,Q),"1,3",HX);
      HX.Append2AndDestroy(FL,H4);
   } 

   H4.makeUnique();
   
   HX.init2DiffOp(A4,1,H4,0);
   HX.Append2AndDestroy(FL,H4);

};




template <class TQ, class TD>
void updateFOps(
    wbvector< QSpace<TQ,TD> > &F12,
    const QSpace<TQ,TD> &A1,
    const QSpace<TQ,TD> &A2,
    const wbvector< QSpace<TQ,TD> > &FC,
    char lflag
){
    unsigned i;
    QSpace<TQ,TD> Xk;

    if (A1.QDIM!=A2.QDIM) wblog(FL,
    "ERR %s() severe data inconsistency (%d,%d,%d)",FCT,A1.QDIM,A2.QDIM);

    F12.init(FC.len);

    for (i=0; i<FC.len; ++i) {

       A2.contract(FL,3,FC[i],2,Xk);
       A1.contract(FL,"1,3;*",Xk,"1,3",F12[i]);
       F12[i].otype=FC[i].otype;

       if (F12[i].isEmpty() && !lflag) wblog(FL,
       "WRN got empty space F12[%d] (iter %d/%d)",i+1);
    }

    if (WbUtil<TD>::isComplex()) wblog(FL,"WRN %s() "
       "expecting AK to be real! (%s)",FCT,getName(typeid(TD)).data
    );
};


template <class TQ, class TD>
void nrgDispIter(
    const unsigned iter,
    const unsigned N,
    const wbvector<TD> &E4, const TD &EK,
    const QSpace<TQ,TD> &HK,
    const char disp
){
    static unsigned D4last=0;
    unsigned i,j, m=HK.QIDX.dim1, n=HK.QDIM;
    wbvector<TQ> q, qmin, qmax;

    wbvector<unsigned> S, I;
    wbarray<unsigned> nn;

    char fmt[8]="%3d ";

    if (!HK.isConsistent()) wblog(FL,"ERR %s",str);
    if (HK.isEmpty()) { m=n=0; }

    qmin.init(n); qmax.init(n); S.init(n);

    for (i=0; i<n; ++i) {
       HK.QIDX.getCol(i,q);
       qmin[i]=q.min(); qmax[i]=q.max();
       S[i]=unsigned(qmax[i]-qmin[i]+1);
    }

    nn.init(S); I.init(n);

    for (i=0; i<m; ++i) {
       for (j=0; j<n; ++j) I[j]=unsigned(HK.QIDX(i,j)-qmin[j]);
       nn(I)+=HK.DATA[i]->SIZE.max();
    }

#ifdef __WBDEBUG__
    i=1;
#else
    if (disp>=0 && (
        (E4.last()==EK) ||
        (iter<3 || iter+3>N) ||
        fabs((double(E4.len)-D4last)/MAX(unsigned(E4.len),D4last)) > 0.20
    )) i=1; else i=0;
#endif

    n=nn.sum();
    if (n!=E4.len)
         sprintf(str,"NK=%d/%ld, EK=%.2f",n, E4.len, EK);
    else sprintf(str,"NK=%d (EK=%.2f)",n, EK);

    wblog(FL,"NRG %02d: Q=[%s; %s]; %s %s",
       iter, qmin.toStr().data, qmax.toStr().data, str, i ? "":"\r\\");
    D4last=MAX(1U,unsigned(E4.len));

    if (disp<=10 || HK.isEmpty()) return;
    if (iter>0) checkGSDeg(E4);
    if (nn.isEmpty() || nn.SIZE.len!=2) return;

    m=nn.SIZE[0];
    n=nn.SIZE[1];

    for (i=0; i<m; ++i) { printf("\n   ");
    for (j=0; j<n; ++j) {
        if (nn(i,j)) printf(fmt,nn(i,j));
        else printf("    ");
    }}

    printf("\n\n"); fflush(0);
};


void checkGSDeg(const wbvector<double> &E4) {
   unsigned i,n;
   double dE=E4.aMin(1);

   for (n=i=0; i<E4.len; ++i) if (ABS(E4[i])<1E-10) ++n;

   if (n==1)
            wblog(FL, " *  ground state is unique (%g).", dE); else
   if (n>1) wblog(FL, " *  ground state is not unique (%d; %g).", n,dE);
   else     wblog(FL, "ERR no ground state found (%d) ???", n);
}


double getPhysE0(const wbvector<double> &DE){

   double E=0.;
   for (unsigned i=0; i<DE.len; ++i) E+=(nrgScale(i)*DE[i]);

   return E;
}


template <class TQ, class TD>
void setupPOP(
   QSpace<TQ,TD> &P,
   const QSpace<TQ,TD> &AK,
   const QSpace<TQ,TD> &UK
){
   if (P.isEmpty()) return;
   wblog(FL,"ERR");
}


