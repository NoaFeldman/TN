char USAGE[] =
/* ======================================================================== */
"   Usage: [omega, A0 [,Info]] = ...                                    \n\
          fdmNRG(NRG_DATA [,Inrg], B, C, Z0, [,opts])                   \n\
                                                                        \n\
     This routine calculates arbitrary correlation function (CF)        \n\
     in the full density matrix Numerical Renormalization Group         \n\
     (fdm-NRG) framework at arbitrary temperature T for arbitrary       \n\
     abelian or non-abelian symmetries.                                 \n\
                                                                        \n\
        A(w) = <B||C'>_w = FFT(<{B(t),C'}>).  (' equiv ^dagger)         \n\
                                                                        \n\
     If B is specified as empty [], it will be automatically assumed    \n\
     equal to C, i.e. B=C as is the case for spectral functions;        \n\
     anti-commutator (default) or commutator will be chosen according   \n\
     to cflags below (see also zflags for Fermionic signs).             \n\
                                                                        \n\
     This routine uses data from a previous standard NRG run with       \n\
     the result written as a matrix product state in QSpace format.     \n\
     It must contain the state space kept (AK) as well as the one       \n\
     thrown / discarded (AT) at every NRG iteration.                    \n\
                                                                        \n\
   See also: NRGWilsonQS                                                \n\
                                                                        \n\
   Input:                                                               \n\
                                                                        \n\
     NRG_DATA   structure containing data of previous NRG run; if this  \n\
                contains path/file, then the data will be read from the \n\
                files path/file_##.mat for every site (see NRGWilsonQS) \n\
                                                                        \n\
      .AK       MPS state of k(ept) space (regular NRG, QSpace)         \n\
      .AT       MPS state of t(runcated) NRG space (QSpace)             \n\
      .HK       effective NRG eigenspectra (kept, QSpace)               \n\
      .HT       effective NRG eigenspectra (truncated, QSpace)          \n\
                                                                        \n\
     Inrg       NRG info structure (required if NRG_DATA contains       \n\
                structure; otherwise <NRG_DATA>_info.mat will be accessed.\n\
                required elements within Inrg structure / data:         \n\
                                                                        \n\
       .Lambda  NRG discretization parameter for the conduction band.   \n\
                Lambda is needed to undo the energy scaling permformed  \n\
                in obtaining the NRG data at the first place.           \n\
                                                                        \n\
       .EScale  energy scale fore each Wilson shell                     \n\
                                                                        \n\
       .E0      (relative) ground state energy offsets in the same      \n\
                energy scales as H[KT]; this is needed for evaluating   \n\
                the correct weights of the density matrix               \n\
                                                                        \n\
     B, C       operators acting on the impurity used to calculate the  \n\
                correlation between them as in A(w)=FFT(<{B(t),C'}>).   \n\
                                                                        \n\
     Z0         Z0 operator for first site (impurity) required for      \n\
                fermionic setup where Z0=(-1)^(number of particles)     \n\
                                                                        \n\
     NB! If the operators B or C are specified by name, they will be    \n\
         read from workspace and stored in NRG_DATA together with       \n\
         the other NRG MPS data (disabled by flag `nostore')            \n\
                                                                        \n\
   Options / Flags (default method: fDM)                                \n\
                                                                        \n\
    'fDM'       new and improved DM-NRG built from full DM (default)    \n\
    'bul'       original way of calculating CF in NRG (R.Bulla)         \n\
    'hof'       original way of calculating CF in DM-NRG (W.Hofstetter) \n\
    'NRho', k   builds density matrix from iteration k only             \n\
    'keven'     include even (odd) shells only for correlation function \n\
    'kodd'      with 1st site (k=0) counting as odd (for bul or hof only)\n\
                                                                        \n\
    'T',...     temperature at which correlation function is calculated;\n\
                default: about an order of magnitude larger than the    \n\
                energy scale at the end of the Wilson chain             \n\
                                                                        \n\
    'partial'   store G1 and G2 separatly (e.g. for detailed balance    \n\
                or the calculation of expectation values)               \n\
    'disp',..   more detailed output (1)                                \n\
    'nostore'   do not store operator elements of B or C with NRG_DATA  \n\
    'locRho'    if Rho needs to be calculated it is stored in RAM only  \n\
                and will be deleted after the program is finished       \n\
    'noRHO'     do not calculate Rho, eg. calculate partition function  \n\
                only (relevant with empty ops B and C only)             \n\
    'rhoNorm',. allows to specify rhoNorm explicitely as input          \n\
                for testing purposes                                    \n\
                                                                        \n\
    'zflags'    on which operators in <B||C> to apply fermionic signs Z0\n\
                ex. 'zflags', [ 1 0 ... ]  would apply Z0 to the first  \n\
                but not to the second operator (default: all 0 if Z0    \n\
                is empty, 1 otherwise).                                 \n\
    'cflags'    where to use regular commutator instead of default      \n\
                anti-commutator for correlation functions (specifically \n\
                relevant if there is a delta peak at omega=0; by default,\n\
                inverse of zflags).                                     \n\
                                                                        \n\
    'nlog',...  number of bins per decade for log. discretization (256) \n\
    'emin',...  minimum energy for omega binning (TN/1000)  ^1)         \n\
    'emax',...  maximum energy for omega binning (10.)      ^1)         \n\
    'binOffset',log. binning starting from |omega|>=binOffset(0.)       \n\
                                                                        \n\
    'symDMA'    symmetrize contributions to spectral function (0)       \n\
    'mspec,'..  number of spectral moments to calculate (0)             \n\
                                                                        \n\
   ^1) NB! first (last) bin contains all data below (above) emin (emax) \n\
                                                                        \n\
   Output                                                               \n\
                                                                        \n\
     om         omega corresponding to the energy binning [-emax,+emax] \n\
     A0         raw data for spectral function (i.e. discrete in om)    \n\
                                                                        \n\
     Info       info structure: specific fields are                     \n\
       .ver     version / mode of code                                  \n\
       .T       actual temperature used                                 \n\
       .rho     FDM weight distribution for given temperature           \n\
       .RHO     reduced thermal density matrix for A0 (R-basis).        \n\
       .Om      Omega(T) = -T*log(Z), (grand)canonical potential.       \n\
       .A4      individual spectral data for each of the two            \n\
                contributions to the (anti)commutator (hence A4 has     \n\
                twice as many  columns as A0 above).                    \n\
       .a4      (1st row) integrated spectral integrals over the two    \n\
                individual contributions to the (anti)commutator of     \n\
                given correlator (2nd row may be ignored; it refers     \n\
                to check w.r.t. detailed balance).                      \n\
       .reA0    contains Re(G(omega->0)) for the calculated correlation \n\
                functions; note that correlators with cflags!=0, e.g.   \n\
                susceptibilities, acquire in addition to the Kramers    \n\
                Kronig transformation (principal value!) also a         \n\
                Matsubara correction for Ea==Eb.                        \n\
       .symfac  IROP factors for given spectral functions.              \n\
                                                                        \n\
   ALTERNATIVE USAGE                                                    \n\
                                                                        \n\
   fdmNRG('--mat','nrgdata', 'setup.mat' [,'B'], 'C')                   \n\
                                                                        \n\
     B and C are the same as above (B=C if B is not specified).         \n\
     This usage allows to read all input parameters from given MAT files\n\
     (MatLab binaries); all variables and options MUST be defined in    \n\
     <setup.mat> by their names shown above; e.g. the operators B and   \n\
     C are specified by their names as defined in <setup.mat>.          \n\
                                                                        \n\
     The output will be written / updated to the file set               \n\
     `<nrgdata>_##.mat' and `<nrgdata>_info.mat'.                       \n\
                                                                        \n\
     Output will be appended to setup.mat.                              \n\
                                                                        \n\
   See also NRGWilsonQS.                                                \n\
                                                                        \n\
   AW (C) May 2006 ; May 2010 ; Nov 2014                                \n\
";

char USAGE_2[] =
"Alternative usage: fdmNRG('nrgdata', 'setup.mat', 'B', 'C')";

/* change log

   changed (C1,C2) or (F,C) to (B,C) consistently
   Wb,Nov24,14 

*/




#define PROG_TAG "FDM"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "fdmNRG"
#endif

#include <float.h>

#define TST_TCHI_0


#define  EPS 1E-12
#define DEPS DBL_EPSILON

#define LOAD_CGC_QSPACE
#include "wblib.h"

#include "nrgdata.hh""nrgdata_v0.hh"
#include "dmrho.cc"
#include "spectral.hh"

template <class TQ, class TD>
void dmNRGIter(
    const char *F_, int L,
    const unsigned iter,
    Spectral &ASpec,
    const QSpace<TQ,TD> &Rho,
    const NRGData<TQ,TD> &H,
    const NRGData<TQ,TD> &B,
    const NRGData<TQ,TD> &C,
    const wbvector<char> &cc,
    const char *XX,
    const char *isbuf=0
);

template <class TQ, class TD>
inline void dmNRGIter(
    const char *F, int L,
    const unsigned iter,
    Spectral &ASpec,
    const QSpace<TQ,TD> &Rho,
    const NRGData<TQ,TD> &H,
    const NRGData<TQ,TD> &C,
    const wbvector<char> &cc,
    const char *XX,
    const char *isbuf=0
){
    dmNRGIter(F,L,iter, ASpec, Rho, H, C, C, cc, XX, isbuf);
}

template <class TQ, class TD>
void dmNRGIter(
    const char *F, int L,
    const unsigned iter,
    Spectral &ASpec,
    const QSpace<TQ,TD> &Rho,
    const QSpace<TQ,TD> &HK,
    const QSpace<TQ,TD> &HT,
    const wbvector< QSpace<TQ,TD> > &BKT,
    const wbvector< QSpace<TQ,TD> > &BTK,
    const wbvector< QSpace<TQ,TD> > &CKT,
    const wbvector< QSpace<TQ,TD> > &CTK,
    const wbvector<char> &cc,
    const char *isbuf=0,
    char dblock=0
);

template <class TQ, class TD>
void updateSpec(
    const char *F, int L,
    Spectral &ASpec,
    const unsigned is,
    QSpace<TQ,TD> &RX,
    const QSpace<TQ,TD> &H1,
    const QSpace<TQ,TD> &H2,
    const double Escale,
    const char *isbuf=0,
    char dblock=0,
    char mflag=0
);

mxArray* addTOA(const char *F, int L, const mxArray *S,
   double T,
   const wbvector<double> &om,
   const wbMatrix<double> &A0,
   const char *Fname="",
   const char *Cname=""
);

   static void CleanUp(void) { Wb::CleanUp(FL); };

   double Lambda=1;

void nrgGetE0(const wbvector<double> &dE, wbvector<double> &E0) {

   if (E0.len!=dE.len) E0.init(dE.len);
   if (!E0.len) { wblog(FL,"WRN %s() got empty E0 !??",FCT); return; }
     


   unsigned k=dE.len-1;
   for (E0[k--]=0; k<dE.len; --k) {
      E0[k]=E0[k+1]-dE[k+1]*nrgScale(k+1);
   }
};


unsigned symDMA=0;


void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   unsigned i,e,k,m,iter, N=0, lflag, NRho=0, dloc=0, Dk, nlog=256;
   char disp=1, store=1, gotops=1, mspec=0;
   char noRHO=0, calcRho=0, locRho=0, calcOps=0, version=4;
   char rzero, keven, kodd, Tflag=0, matflag=0, partial=0, PARTIAL=0;
   double TN, Om=0, dbl=0, emin=-1, emax=-1, Delta=0., T=-1;
   const unsigned r1=2, r2=4;

   wbMatrix<double> A0;
   wbvector<double> w, a1, a2;
   wbvector<unsigned> D;

   Spectral ASpec;
   mxArray *S, *a; 
   const mxArray *ac, *Inrg=NULL;

   wbstring vstr(32), fname(128);
   char vtag[4], vbuf[4], itag[16];

   NRGData<gTQ,gTD> H("H"), A("A"), RHO("RHO"), B("B"), C("C");
   wbvector< QSpace<gTQ,gTD> > B0, C0;
   QSpace<gTQ,gTD> Rho, Z0;

   wbMatrix<gTD> rhoNorm, rhoNormI, dd;
   wbvector<double> om,DE,E0;
   wbvector<char> zz,cc;
   OPTS opts;

   wbSigHandler SIG(FL);

#ifdef WB_CLOCK
   WbClock wbc_dmNRG("fdmNRG");
   wbc_dmNRG.start();
#endif

   NRG_N=0; NRG_ITER=0; gES.init(); 

   if (nargin && isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin<2) wblog(FL,"ERR invalid usage");

   if (mxIsChar(argin[0])) {
      mxGetString(argin[0],str,8); str[7]=0;
      if (!strcmp(str,"--mat")) { matflag=1; --nargin; ++argin; }
   }

#ifdef MATLAB_MEX_FILE
   mexAtExit(CleanUp);
#endif
   Wb::CleanUp(NULL,1);


   if (matflag) {

      if (nargout || (nargin!=3 && nargin!=4)) wblog(FL,
      "ERR invalid number of I/0 arguments (%d)",nargin);

      for (i=0; i<(unsigned)nargin; ++i) if (!mxIsChar(argin[i]))
      wblog(FL,"ERR invalid arg #%d %N%N%s%N", i+1, USAGE_2);

      i=mxGetString(argin[0],fname.data,fname.len-1);
      if (i || fname.isEmpty()) wblog(FL,
      "ERR reading string of arg #1 (%d)", fname.len);

      if (mxGetString(argin[1],str,MAX_STRLEN))
      wblog(FL,"ERR error reading string of arg #2 (%d)", MAX_STRLEN);

      opts.init(str);"options" by reference to file

      if (nargin==4) {
         k=2;
         if (mxGetString(argin[k], str, MAX_STRLEN)) wblog(FL,
            "ERR error reading string of arg #%d",k+1);
         B.initName(str);
         opts.getOpt(str,a,'!'); if (a) mxInitQSpaceVec(FL,a,B0,0,&r1,&r2);
      }

      k=nargin-1;
      if (mxGetString(argin[k],str,MAX_STRLEN)) wblog(FL,
         "ERR error reading string of arg #%d", k+1);
      C.initName(str);
      opts.getOpt(str,a,'!'); if (a) mxInitQSpaceVec(FL,a,C0,0,&r1,&r2);

      opts.getOpt("Z0",a); if (!a) opts.getOpt("Z",a,'!');
      if (a) { unsigned r=-1;
         if (!mxIsQSpace(FL,a,r)) wblog(FL,"ERR invalid Z0 (%s)", str);
         Z0.init(FL,a);
      }
   }

   else {

      if (nargout>3 || nargin<4) wblog(FL,
         "ERR invalid number of I/O arguments (%d/%d)",nargout,nargin);

      if (mxIsChar(argin[0])) {
         i=mxGetString(argin[0], fname.data, fname.len-1);
         if (i || fname.isEmpty()) wblog(FL,
         "ERR error reading file tag from arg #1");
      }

      ac=argin[1];
      if (!ac || mxGetNumberOfElements(ac)!=1
         || mxGetFieldNumber(ac,"Lambda")<0
         || mxGetFieldNumber(ac,"EScale")<0
         || mxGetFieldNumber(ac,"E0"    )<0
      ){
         k=1; if (fname.isEmpty())
         wblog(FL,"ERR valid Inrg structure required (2nd argument)");
      }
      else { Inrg=ac; k=2; }

      for (i=1; i<=2; ++i, ++k) {
         if (mxIsEmpty(argin[k])) {
            continue;
         }

         if (mxIsChar(argin[k])) { try {
            if (mxGetString(argin[k],str,MAX_STRLEN))
            wblog(FL,"ERR error reading string of arg #%d",k+1);

            ac=mexGetVariablePtr("caller",str);

            if (!ac) wblog(FL,
              "ERR reading QSpace `%s' (arg #%d) from workspace",str,k+1);

            if (i==1)
                 { mxInitQSpaceVec(FL,ac,B0,0,&r1,&r2); B.initName(str); }
            else { mxInitQSpaceVec(FL,ac,C0,0,&r1,&r2); C.initName(str); }
            }
            catch (...) {
               wblog(FL,"ERR reading %s operator (arg #%d)",
               i==1 ? "1st":"2nd", k+1);
            }
         }
         else { try {
            wbvector< QSpace<gTQ,gTD> > &X=(i==1 ? B0 : C0);
            mxInitQSpaceVec(FL,argin[k],X,0,&r1,&r2);
          } catch (...) {
            wblog(FL,"ERR reading 1st operator (arg #%d)",k+1);
          }}
      }
      Z0.init(FL,argin[k++]);

      opts.init(argin+k, nargin-k);
   }

   opts.getOpt("disp", disp);
   if (!disp) {
        WbClock XC("X"); XC.set_TSC_TO_SEC(); }
   else printf("\n");

   if (!fname.isEmpty()) {
   }
   else {
      wblog(FL,"<i> internal mode (no data I/O)");
      if (nargout<2) wblog(FL,
      "ERR fdmNRG requires at least two output arguments (%d)",nargout);
   }

   A  .setupIO(FL, fname.data, argin[0], disp ? 'i' : 0);
   H  .setupIO(FL, fname.data, argin[0]);
   RHO.setupIO(FL, fname.data, argin[0]);
   B  .setupIO(FL, fname.data, argin[0]);
   C  .setupIO(FL, fname.data, argin[0]);

   if (fname.isEmpty()) {
      if (!Inrg) { opts.getOpt("Inrg", a); Inrg=(const mxArray*)a; }
      if (!Inrg) wblog(FL,"ERR Inrg structure required (2nd argument)\n"
         "since no NRG file structure was specified");
      if (mxGetNumberOfElements(Inrg)!=1) wblog(FL,
         "ERR invalid input Inrg (scalar structure expected)");
      wblog(FL," *  reading Lambda, E0, EScale from Inrg structure");

      ac=mxGetField(Inrg,0,"Lambda"); if (ac) mxGetNumber(ac,Lambda);
      else wblog(FL,"ERR invalid input Inrg (missing field 'Lambda')");

      ac=mxGetField(Inrg,0,"E0"); if (ac) DE.init(FL,ac);
      else wblog(FL,"ERR invalid input Inrg (missing field 'E0')");

      ac=mxGetField(Inrg,0,"EScale"); if (ac) gES.init(FL,ac);
      else wblog(FL,"ERR invalid input Inrg (missing field 'EScale')");
   }
   else {
      if (disp) {
         for (i=fname.len; i>0; --i) { if (fname[i]=='/') { ++i; break; }}
         wblog(FL," *  reading Lambda, E0, EScale from %s_info.mat",
         fname.data+i);
      }

      ac=RHO.getMxInfo("Lambda"); if (ac) mxGetNumber(ac,Lambda);
      else wblog(FL,"ERR failed to read Lambda from <NRG>_info.mat");

      ac=RHO.getMxInfo("E0"); if (ac) DE.init(FL,ac);
      else wblog(FL,"ERR failed to read E0 from <NRG>_info.mat");

      ac=RHO.getMxInfo("EScale"); if (ac) gES.init(FL,ac);
      else wblog(FL,"ERR failed to read EScale from <NRG>_info.mat");
   }


   if (Lambda<=1.) wblog(FL,"ERR invalid Lambda=%g.",Lambda);

   H.checkVec(FL,"K",2, "NRG data (HK)");

   N=NRG_N;
   
   if (N==0) wblog(FL,"ERR empty NRG data set !??");
   if (N!=DE.len || N!=gES.len) wblog(FL,
      "ERR mismatch in chain length (E0[%d], EScale[%d], N=%d)",
       DE.len,gES.len, N);

   SIG.call99();
   dloc=getDLoc(A,&SIG,&Dk, 0);

   A.init(FL,"K",0);

   TN=nrgScale(N-1);


   store =!opts.getOpt("nostore");
   locRho =opts.getOpt("locRho" );
   noRHO  =opts.getOpt("noRHO"  );
   calcRho=opts.getOpt("calcRho"); if (locRho) calcRho=1;
   calcOps=opts.getOpt("calcOps");

   opts.getOpt("nlog", nlog);
   opts.getOpt("emin", emin);
   opts.getOpt("emax", emax);
   opts.getOpt("binOffset",Delta);
   opts.getOpt("NRho", NRho);

   keven=opts.getOpt("keven");
   kodd =opts.getOpt("kodd" ); m=0;

   if (opts.getOpt("bul")) { ++m; version=0; }
   if (opts.getOpt("hof")) { ++m; version=1; }
   if (opts.getOpt("fDM")) { ++m; version=4; }

   if (m>1) wblog(FL,
   "ERR multiple definitions of methods (%d)",m);

   symDMA=opts.getOpt("symDMA");
   opts.getOpt(FL,"mspec",mspec);

   opts.getOpt("cflags",a); if (a) { cc.init(FL,a,"cflags"); }
   opts.getOpt("zflags",a); if (a) { zz.init(FL,a,"zflags"); }

   if (NRho) {
      if (NRho>N) wblog(FL,"ERR NRho out of bounds (%d/%d)", NRho, N);
      if (version>3) version=3;
   }
   else if (version<4) NRho=N;

#ifndef MATLAB_MEX_FILE
   if (Wb::GetEnv(FL,"DMA_T",T)>=0)
   Tflag=1;
#endif

   opts.getOpt("T", T);
   partial=opts.getOpt("partial");
   PARTIAL=opts.getOpt("PARTIAL");

   opts.getOpt("rhoNorm",a);
   if (a) rhoNormI.init(FL,a);

   opts.checkAnyLeft();

   if (emin<=0) {
      emin=1E-3*pow(Lambda,-double(N)/2);
   }

   if (1) {
      double e1=0., e2=0.;

      H.init(FL,"K",0); e1=getErange(H.K)*nrgScale(0);
      H.init(FL,"K",1); e2=getErange(H.K)*nrgScale(1); dbl=MAX(e1,e2);

      if (emax<=0) {
         if (dbl<10) { emax=10; }
         else { emax=10*ceil(0.12*dbl); }
      }
      else { e1=emax; lflag=1;
         if (emax<10) { emax=10; } else
         if (emax<1.2*dbl) {
            emax=10*ceil(0.12*dbl);
         }
         else { lflag=0; }

         if (lflag) wblog(FL,
         "WRN adjusting emax=%g => %g (having E<=%.3g)",e1,emax,dbl); 
      }
   }

   if (disp) {
      strcpy(str,PROG);
      wblog(FL,"=== %s",strpad(str,'=',48)); str[0]=0;
      if (disp>1 ) sprintf(str,", disp=%d", disp);
      if (!store ) strcat(str, ", nostore");
      if (locRho ) strcat(str, ", locRho");
      if (calcOps) strcat(str, ", calcOps");
      if (symDMA ) strcat(str, ", symDMA");
      if (version==0) strcat(str, ", bul");
      if (version==1) strcat(str, ", hof");
      if (keven  ) strcat(str, ", keven");
      if (kodd   ) strcat(str, ", kodd");
      if (partial) strcat(str,", partial");
      if (PARTIAL) strcat(str,", PARTIAL");

      wblog(FL," *  Lambda=%g, L=%d, D=%d, dloc=%d\n *  sym=%s",
         Lambda,N,Dk, dloc, A.K.qtype.toStr('V').data);
      if (str[0]) wblog(FL," *  %s",str+2); fflush(0);
   }


   switch (version) {
     case 0:

       vstr.cpy(FL,"bare NRG (R. Bulla)");
       strcpy(vtag,"NRG"); break;

     case 1:

       if (keven || kodd)
       wblog(FL,"WRN keven/kodd with conventional fdmNRG !??");

       if (NRho==N)
            vstr.cpy(FL,"original DM-NRG (W. Hofstetter)");
       else vstr.printf(FL,"DM-NRG (W. Hofstetter; N=%d/%d)",NRho,N);
       strcpy(vtag,"HOF"); break;

     case 3:

       vstr.printf(FL,"fdm-NRG (single shell NRho=%d)", NRho);
       strcpy(vtag,"DMK"); break;

     case 4:

       vstr.cpy(FL,"fdm-NRG (full density matrix)");
       strcpy(vtag,"FDM"); break;

     default: wblog(FL,"ERR invalid version %d.", version);
   }

   if (version<=1) strcpy(vbuf,"buf"); else vbuf[0]=0;

   if ((keven || kodd) && version>1) {
      wblog(FL,"WRN %s ignored unless using conventional NRG\n"
        "hof or bul; using %s (%d)", keven ? "keven":"kodd",
         vstr.data, version
      ); keven=kodd=0;
   }
   else if (keven && kodd) {
      wblog(FL,"WRN Both keven and kodd are set (ignoring kodd)");
      kodd=0;
   }

   if (version>=4) {
      dbl=sqrt(2/(log(Lambda)*log(double(dloc))));
      m=(unsigned)ceil(5*dbl);

      k=( m+1<N ? N-m-1 : 0 );
      if (int(m)<0) { m=0; }; dbl=nrgScale(k);
      if (T<0) { T=dbl;
         if (disp) wblog(FL," *  using default T=%.4g (N-%d)", T, N-k);
      }
      else if (T>0 && T<dbl)
      wblog(FL,"WRN T=%.3g < T_{N-%d}=%.3g", T, N-k, dbl);
   }
   else if (T<0) T=0;

   if (disp) {
      str[0]=0; if (Delta!=0) sprintf(str,", binOffset=%g",Delta);
      wblog(FL,
      " *  %48R\n"
      " *  method      : %s\n"
      " *  temperature : %.4g (= %.4g TN)\n"
      " *  omega range : %.4g .. %.4g (%d/dec%s)\n"
      "=== %48R", "-",
      vstr.data, T, T/TN, emin, emax, nlog, str, "="); fflush(0);
   }

   SIG.call99();
   if (version>1) {
      A.checkVec(FL,"T",3,"NRG data (AT)"); SIG.call99();
      H.checkVec(FL,"T",2,"NRG data (HT)"); SIG.call99();
   }


   if (C0.isEmpty()) gotops=0;

   if (B0.len && B0.len!=C0.len) wblog(FL,
      "ERR Length mismatch of 1st with 2nd operator (%d,%d).",
       B0.len, C0.len);
   if (B0.len && B0[0].QDIM!=C0[0].QDIM) wblog(FL,
      "ERR QDIM mismatch of 1st with 2nd operator (%d,%d).",
       B0[0].QDIM, C0[0].QDIM);

   C.initOp(FL, C0, version>1, NULL, "op2=C");
   B.initOp(FL, B0, version>1, &C,   "op1=B");

   if (calcOps){
      C.forceCalc();
      B.forceCalc();
   }
   if (!store) { C.store=0; if (B.nrgIdx.len) B.store=0; }

   if (disp) {
      if (!gotops) 
         wblog(FL," *  got empty ops: calculate FDM/Z only");
      else {
         if (B.nrgIdx.len)
         B.info(FL ); else {
            wblog(FL," *  B = C (calculating <C(t)||C'>)"); }
         C.info(FL );
      }
   }


   if (C.store) C.updatePara(FL,"KK", 0, "KK_");
   if (B.store) B.updatePara(FL,"KK", 0, "KK_");

   if (disp) {
      if (!C.calc && !B.calc)
      wblog(FL,"<i> using stored operator set.");
   }

   if (B.nrgIdx.len && B.KK.isEqual(C.KK)) {
      wblog(FL,"<i> got B==C (-> compute C only)");
      B.nrgIdx.init();
   }


   if (!zz.len || !cc.len) {
      if (cc.len) {
         zz=cc; for (i=0; i<zz.len; ++i) zz[i]=(cc[i] ? 0:1);
      }
      else if (zz.len) {
         cc=zz; for (i=0; i<cc.len; ++i) cc[i]=(zz[i] ? 0:1);
      }
      else {
         zz.init(C0.len); cc.init(C0.len);
         if (Z0.isEmpty()) cc.set(1); else zz.set(1);
      }
   }

   if (zz.len!=C0.len) {
      if (zz.len==1) { zz.Resize(C0.len); zz.set(zz[0]); } else wblog(FL,
      "ERR length mismatch of zflags with [BC] (%d/%d)",zz.len,C0.len);
   }
   if (cc.len!=C0.len) {
      if (cc.len==1) { cc.Resize(C0.len); cc.set(cc[0]); } else wblog(FL,
      "ERR length mismatch of cflags with [BC] (%d/%d)",cc.len,C0.len);
   }

   if (disp) {
      if (zz^cc) {
         if (zz.len==1)
              wblog(FL," *  zflags=%d", zz[0]);
         else wblog(FL," *  zflags=[%s]", zz.toStrf("","").data);
      }
      else
         wblog(FL," *  zflags=[%s], cflags=[%s]",
         zz.toStrf("","").data, cc.toStrf("","").data
      );
   }

   if (!Z0.isEmpty()) { e=0;
      try {
         if (C.calc>0) {
            for (i=0; i<C.KK.len && !e; ++i)
            if (zz[i]) e+=applyZ0(Z0,C.KK[i]);
         }
         if (B.calc>0) {
            for (i=0; i<B.KK.len && !e; ++i)
            if (zz[i]) e+=applyZ0(Z0,B.KK[i]);
         }
         if (e) wblog(FL,"ERR %s() e=%d",FCT,e);
      }
      catch(...) {
         wblog(FL,
           "ERR QIDX of Z0 does not match Q-spaces of %s\n"
           "Hint: unset Z0 to empty if not relevant.",
            e==1 ? C.name.data : B.name.data
         );
      }
   }
   else
   if (zz.anyUnequal(0)) wblog(FL,"ERR got zflags for empty Z0!");

   ASpec.init(
     "A", 2*C.KK.len, emin, emax, nlog, Delta, T,
      PARTIAL ? N:0
   );
   if (mspec) {
      ASpec.calcSpectralMoments(FL,mspec);
   }

   if (version<=1) ASpec.initBuf();


   nrgGetE0(DE,E0);

   if (version>=1) {
      double Eref=E0[NRho>0 ? NRho-1 : N-1];
      E0-=Eref;
      Om=initRHO(rhoNorm, T, H, E0, dloc, NRho, disp) + Eref;
   }
   else {
      H.init(FL,"K",0); getEData(H.K,a1,D);
      E0.set(0); E0[0]=-a1.min()*nrgScale(0);

      initRHO(rhoNorm, T, H, E0, 0, NRho, disp);
   }
   SIG.call99();

   if (!rhoNormI.isEmpty()) {
      if (rhoNorm.hasSameSize(rhoNormI)) {
         wblog(FL,"RHO init rhoNorm from input args");
         rhoNorm=rhoNormI;
      }
      else wblog(FL,
      "ERR size mismatch of input rhoNorm\n%dx%d vs. %dx%d",
       rhoNorm.dim1,rhoNorm.dim2,rhoNormI.dim1,rhoNormI.dim2);
   }

   if (version && (gotops || !noRHO)) { wbtop().runningLarge(FL);

      dd.init(RHO.getMxInfo("rhoNorm"));

      if (calcRho || dd!=rhoNorm || RHO.checkVec(FL,"",2) ||
         (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T)
      ){
         int q=-1;
            if (calcRho) q=1; else
            if (dd!=rhoNorm) q=2; else
            if (RHO.checkVec(FL,"",2)) q=3; else
            if (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T) q=4;

         wblog(FL,"%s %scalculate density matrices RHO (%d) \\",
            vtag, !dd.isEmpty() ? "(re)":"",q);

         if (locRho && !RHO.MX) {
            printf(" internally (locRho=%d)",locRho);
            RHO.MX=mxCreateStructMatrix(1,NRG_N,0,NULL);
            RHO.MXloc=1; locRho=99;
         }

         printf("\n"); SIG.call99();
         nrgUpdateRHO(rhoNorm, A, H, RHO, NRho, dloc, &SIG);
         printf("\n"); SIG.call99();
      }
      else wblog(FL,"%s using stored density matrices",vtag);
   }
   
   wbtop().runningLarge(FL);

   if (locRho!=99) locRho=0;


   if (disp) {
      wblog(FL,"TST CG_EPS     =[%g, %g] \r\\",CG_EPS1, CG_EPS2);
      wblog(FL,"TST CG_SKIP_EPS=[%g, %g] \r\\",CG_SKIP_EPS1, CG_SKIP_EPS2);
   }

   if (gotops)
   for (iter=0; iter<N; ++iter) { SIG.call99();

      NRG_ITER=iter;
      rhoNorm.getRec(iter,w); rzero=(w<DEPS);
      
      H.init(FL,"K",iter);
      H.init(FL,"T",iter);

      if (C.calc>0 || B.calc>0) {
         A.init(FL,"K",iter);
         A.init(FL,"T",iter);
      }

      m=4;

      C.updateOp(FL,A,m,iter);
      B.updateOp(FL,A,m,iter);

      C.storeOp(FL,iter,m);
      B.storeOp(FL,iter,m);

      snprintf(itag,16,"%s %02d",vtag,iter+1); lflag=0;

      if (version<=1 && H.T.isEmpty()) {
         wblog(FL,"%s: all kept -> skip to next iteration%s",
         itag, disp<2 ? "\r\\" : "");
         continue;
      }

      if (version==0) {
         wblog(FL,"%s: (plain NRG)%-30s\r\\", itag,
         C.calc>0 || B.calc>0 ? " + operator update" : ""); ++lflag;
      }

      if (w[1]>=DEPS) {
         if (version) {
            if (disp>2)
                 wblog(FL,"%s: (RhoT)", itag);
            else wblog(FL,"%s: (RhoT)%30s\r\\", itag,"");
            ++lflag;
         }

         initRho(Rho, H.T, iter, w[1]);
         if (B.nrgIdx.len==0) {
            dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "TT", vbuf);
            dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "TK", vbuf);
         }
         else {
            dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "TT", vbuf);
            dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "TK", vbuf);
         }
      }
      else if (H.T.isEmpty()) {
         if (disp>2)
              wblog(FL,"%s: no discarded space ", itag);
         else wblog(FL,"%s: no discarded space %8s\r\\",itag,"");
         ++lflag;
      }

      if (w[0]>=DEPS) {
         if (version) {
            wblog(FL,"%N%s: include kept space (%.3g)", itag, w[0]);
            ++lflag;
         }

         initRho(Rho, H.K, iter, w[0]);

         if (B.nrgIdx.len==0) {
            dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "KK", vbuf);
            dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "KT", vbuf);
         }
         else {
            dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "KK", vbuf);
            dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "KT", vbuf);
         }
      }

      if (rzero && !NRho && version!=1 && !H.T.isEmpty()) {
         if (disp>2)
              wblog(FL,"%s: (dRho=%.3g)",itag,w.sum());
         else wblog(FL,"%s: (dRho=%.3g)%8s\r\\",itag,w.sum(),"");
         ++lflag;
      }

      if (version>=1) {
         RHO.init(FL,"A",iter);

         if (!lflag) {
            if (RHO.A.isEmpty()) {
               sprintf(str,"%s: (EMPTY RHO)",itag);
            }
            else {
               sprintf(str,"%s: (RHO)",itag);
               ++lflag;
            }

            if (disp>2) wblog(FL,"%s",str);
            else wblog(FL,"%-50s\r\\",str);
         }

         if (!RHO.A.isEmpty()) {
         if (!H.T.isEmpty()) {
            if (B.nrgIdx.len==0)
                 dmNRGIter(FL,iter, ASpec, RHO.A, H,    C, cc, "KT", vbuf);
            else dmNRGIter(FL,iter, ASpec, RHO.A, H, B, C, cc, "KT", vbuf);
         }
         if (version==1) {
            if (B.nrgIdx.len==0)
                 dmNRGIter(FL,iter, ASpec, RHO.A, H,    C, cc, "KK", vbuf);
            else dmNRGIter(FL,iter, ASpec, RHO.A, H, B, C, cc, "KK", vbuf);
         }}
      }

      if (vbuf[0]) {

         if ((!keven || (iter%2)==1) && (!kodd || (iter%2)==0)) {
            ASpec.crossAddBuf(2, iter<NRho && version);
         }
         else {
            sprintf(str,
            "%s: skipped (%s)", itag, keven ? "keep even" : "keep odd");
            if (disp>2)
                 wblog(FL,"%s",str);
            else wblog(FL,"%s \r\\",str);
         }
         ASpec.initBuf();
      }

      if (PARTIAL) ASpec.saveIter(FL,iter);

      if (!lflag) {
         if (!B.calc && !C.calc) break;
         if (disp>2)
              wblog(FL,"%s: (update operators)",itag);
         else wblog(FL,"%s: (update operators)%20s\r\\",itag,"");
         ++lflag;
      }
      doflush();
   }

   if (disp>1) printf("\n");


   ASpec.applyIROPfac(
   B.TT.len ? B.TT : (B.KK.len ? B.KK : (C.TT.len ? C.TT : C.KK)), &a1);

   S=mxCreateStructMatrix(1,1,0,NULL);
   mxAddField2Scalar(FL,S,"ver",    vstr.toMx());
   mxAddField2Scalar(FL,S,"stamp",  wbTimeStamp().toMx());

   mxAddField2Scalar(FL,S,"paras", MXPut(0,0)
      .add(B0,   "B"        )
      .add(C0,   "C"        )
      .add(zz,   "zflags"   )
      .add(cc,   "cflags"   )
      .add(E0,   "E0"       )
      .add(gES,  "EScale"   )
      .add(emin, "emin"     )
      .add(emax, "emax"     )
      .add(nlog, "nlog"     )
      .add(Delta,"binOffset")
      .add(EPS,"eps").add(DEPS,"deps").add(CG_VERBOSE,"cg_verbose")
      .add(CG_EPS1,"cg_eps1").add(CG_EPS2,"cg_eps2")
      .add(CG_SKIP_EPS1,"cg_skip_eps1").add(CG_EPS2,"cg_skip_eps2")
   .toMx());

   RHO.init(FL,"A",0);
   A.init(FL,"K",  0);

   mxAddField2Scalar(FL,S,"A0",    A.K.toMx());
   mxAddField2Scalar(FL,S,"RHO",   RHO.A.toMx());
   mxAddField2Scalar(FL,S,"rho",   rhoNorm.toMx());
   mxAddField2Scalar(FL,S,"Om",    numtoMx(Om));
   mxAddField2Scalar(FL,S,"T",     numtoMx(T));
   mxAddField2Scalar(FL,S,"symfac",a1.toMx());


   ASpec.getRawData(om,A0);
   if (partial) mxAddField2Scalar(FL,S,"A4",A0.toMx('r'));
   if (PARTIAL) {
      wbarray<double> A1,A2;
      ASpec.getPARTIAL(A1,A2);
      mxAddField2Scalar(FL,S,"A1",A1.toMx());
      mxAddField2Scalar(FL,S,"A2",A2.toMx());
   }
   a1=A0.recSum();

   ASpec.pairRawData();
   ASpec.detailedBalance(T,om,A0);

   if (disp)
   ASpec.dispSumRule();

   if (partial) mxAddField2Scalar(FL,S,"A4f",A0.toMx('r'));
   a2=A0.recSum();

   if (a1.len==a2.len) {
      A0.init(2,a1.len);
      A0.recSetP(0,a1.data);
      A0.recSetP(1,a2.data);

      mxAddField2Scalar(FL,S,"a4",A0.toMx());
   }
   else wblog(FL,
      "WRN severe size inconsistency (%d,%d) !??",a1.len,a2.len);

   mxAddField2Scalar(FL,S,"reA0", ASpec.Ar.toMx());
   mxAddField2Scalar(FL,S,"mspec",ASpec.mspec.toMx());

   a=A0.toMx();
   ASpec.getRawData(om,A0);

   if (disp) {
      doflush();
   }

   if (!fname.isEmpty()) {
      C.updateInfo(FL, "om",      om.toMx('t'));
      C.updateInfo(FL, "a0",      A0.toMx('r'));
      C.updateInfo(FL, "a4",      a);
      C.updateInfo(FL, "ISpec",   S, 'k');

      if (!locRho) {
      C.updateInfo(FL, "rhoT",    numtoMx(T));
      C.updateInfo(FL, "rhoNorm", rhoNorm.toMx());
      }

      if (Tflag) {
         a=addTOA(FL,
           C.getMxInfo("TOA"), T, om, A0, B.name.data, C.name.data
         );
         C.updateInfo(FL,"TOA",a);
      }
   }

   if (nargout>=2) {
      argout[0] = om.toMx('t');
      argout[1] = A0.toMx('r');
   }

   if (nargout>2)
        argout[2]=S;
   else mxDestroyArray(S);

  { mxArray *aM=gCX.toMx();
    mexPutVariable("caller","gCX_fdm",aM);
    mxDestroyArray(aM);
  }

   if (disp) {
      doflush();
   }

   if (disp) {
      wblog(FL,"FIN NRGData I/O time usage: %s",
      tic_nrgIO.gettime("asm").data);
   }

   tic_nrgIO.reset();

#ifdef WB_CLOCK
   wbc_dmNRG.stop();
#endif

   if (disp) { CleanUp(); printf("\n"); }

   if (disp) {
      doflush();
   }

#ifdef __WB_MEM_CHECK__
   A0.init(); w.init(); a1.init(); a2.init();
   D.init(); fname.init();
   H.init(); A.init(); RHO.init(); Rho.init(); Z0.init();
   B.init(); C.init(); B0.init(); C0.init();
   ASpec.init(); opts.init(); gES.init();

   rhoNorm.init(); rhoNormI.init(); dd.init();
   om.init(); DE.init(); E0.init();
   zz.init(); cc.init();

   getBoltzman(E0,D,E0,-1,-1);
   tic_nrgIO.Init();

   Wb::MemCheck(FL,"LIST");
#endif


   if (disp) {
      doflush();
   }

}


mxArray* addTOA(const char *F, int L, const mxArray *S,
   double T,
   const wbvector<double> &om,
   const wbMatrix<double> &A0,
   const char *Fname,
   const char *Cname
){

   mxArray *a, *S2=mxCreateStructMatrix(1,1,0,NULL);
   const char* fn[] = {"T","Om","A0","B","C"};
   unsigned nf=5;

   if (!S) {
      mxAddField2Scalar(F,L, S2, fn[0], numtoMx( T ));
      mxAddField2Scalar(F,L, S2, fn[1], om.toMx('t'));
      mxAddField2Scalar(F,L, S2, fn[2], A0.toMx('r'));
      mxAddField2Scalar(F,L, S2, fn[3], mxCreateString(Fname ? Fname :""));
      mxAddField2Scalar(F,L, S2, fn[4], mxCreateString(Fname ? Cname :""));
   }
   else {
      unsigned i,j, m=mxGetM(S), n=mxGetN(S);
      int id[nf], e=mxGetNumberOfDimensions(S)>2 || (m!=1 && n!=1);

      m*=n;

      for (i=0; i<nf; ++i) { 
         id[i]=mxGetFieldNumber(S,fn[i]);
         if (id[i]<0) ++e;
      }

      if (e) wblog(F,L,"ERR invalid structure TOA !??");

      n=mxGetNumberOfFields(S);
      S2=mxCreateStructMatrix(m+1,1,0,NULL);
      for (i=0; i<n; ++i) mxAddField(S2,mxGetFieldNameByNumber(S,i));

      for (i=0; i<m; ++i)
      for (j=0; j<n; ++j) {
          a=mxGetFieldByNumber(S,i,j); if (!a) continue;
          mxSetFieldByNumber(S2,i,j,mxDuplicateArray(a));
      }

      mxSetFieldByNumber(S2, i, id[0], numtoMx( T ));
      mxSetFieldByNumber(S2, i, id[1], om.toMx('t'));
      mxSetFieldByNumber(S2, i, id[2], A0.toMx('r'));
      mxSetFieldByNumber(S2, i, id[3], mxCreateString(Fname ? Fname :""));
      mxSetFieldByNumber(S2, i, id[4], mxCreateString(Fname ? Cname :""));
   }

   return S2;
}


template <class TQ, class TD>
inline void dmNRGIter(
    const char *F_, int L,
    const unsigned iter,
    Spectral &ASpec,
    const QSpace<TQ,TD> &Rho,
    const NRGData<TQ,TD> &H,
    const NRGData<TQ,TD> &B,
    const NRGData<TQ,TD> &C,
    const wbvector<char> &cc,
    const char *XX,
    const char *isbuf
){
    if (!strcmp(XX,"KK"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.K, H.K,
       B.KK, B.KK,
       C.KK, C.KK, cc, isbuf, 'K'
    ); else
    if (!strcmp(XX,"KT"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.K, H.T,
       B.KT, B.TK,
       C.KT, C.TK, cc, isbuf
    ); else
    if (!strcmp(XX,"TK"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.T, H.K,
       B.TK, B.KT,
       C.TK, C.KT, cc, isbuf
    ); else
    if (!strcmp(XX,"TT"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.T, H.T,
       B.TT, B.TT,
       C.TT, C.TT, cc, isbuf, 'D'
    );
    else wblog(FL,"ERR invalid tag `%s'", XX);
};

template <class TQ, class TD>
void dmNRGIter(
    const char *F, int L,
    const unsigned iter,
    Spectral &ASpec,
    const QSpace<TQ,TD> &Rho,
    const QSpace<TQ,TD> &HK,
    const QSpace<TQ,TD> &HT,
    const wbvector< QSpace<TQ,TD> > &BKT,
    const wbvector< QSpace<TQ,TD> > &BTK,
    const wbvector< QSpace<TQ,TD> > &CKT,
    const wbvector< QSpace<TQ,TD> > &CTK,
    const wbvector<char> &cc,
    const char *isbuf,
    char dblock
){
   unsigned i, d=BKT.len, e=0;
   QSpace<TQ,TD> RX;

   const double Escale = nrgScale(iter);

   if (HK.isEmpty() || HT.isEmpty()) {
      for (i=0; i<BTK.len; ++i) if (!BTK[i].isEmpty()) ++e;
      for (i=0; i<BKT.len; ++i) if (!BKT[i].isEmpty()) ++e;
      for (i=0; i<CTK.len; ++i) if (!CTK[i].isEmpty()) ++e;
      for (i=0; i<CKT.len; ++i) if (!CKT[i].isEmpty()) ++e;
      if (e) { dbstop(FL); wblog(FL,
      "ERR Inconsistency in empty data space (%d) !!", e); }

      return;
   }

   if (BTK.len!=d || CTK.len!=d || CKT.len!=d) wblog(FL,
      "ERR severe length mismatch (B-ops(%d): %d/%d; %d/%d)",
       iter+1, BKT.len, BTK.len, CKT.len, CTK.len);
   if (cc.len!=d) wblog(FL,
      "ERR cc length mismatch (%d/%d) ???", cc.len, d);
   if (typeid(TD)!=typeid(double)) wblog(F,L,
      "ERR ops expected to be real (%s)",getName(typeid(TD)).data);

   if (Rho.isEmpty()) return;

   for (i=0; i<d; ++i) {
       Rho.contract(2, CKT[i],1, RX);
       RX.TimesEl(BKT[i]);

       if (symDMA && BKT[i]!=CKT[i]) {
          QSpace<TQ,TD> RY;

          Rho.contract(1, BKT[i],1, RY);
          RY.TimesEl(CKT[i]);

          RX+=RY; RX*=0.5;
       }

       updateSpec(F_L, ASpec,
       2*i, RX, HK, HT, Escale, isbuf, dblock, cc[i]);
   }

   for (i=0; i<d; ++i) {
       wbperm P; if (CTK[i].isR3Op()) P.initStr("1,3,2");

       CTK[i].contract(2,Rho,1,RX,P);
       RX.TimesEl(BTK[i]);

       if (symDMA && BKT[i]!=CKT[i]) {
          QSpace<TQ,TD> RY;

          BTK[i].contract(2, Rho,2, RY);
          RY.TimesEl(CTK[i]);

          RX+=RY; RX*=0.5;
       }

       if (cc[i]) RX*=-1;

       updateSpec(F_L, ASpec,
       2*i+1, RX, HT, HK, Escale, isbuf, dblock, cc[i]);
   }
};


template <class TQ, class TD>
void updateSpec(
    const char *F, int L,
    Spectral &ASpec,
    const unsigned is,
    QSpace<TQ,TD> &RX,
    const QSpace<TQ,TD> &H1,
    const QSpace<TQ,TD> &H2,
    const double Escale,
    const char *isbuf,
    char dblock,
    char mflag
){
    unsigned i,j,k,m,n,r,s, tobuf=0, e=0; INDEX_T m1,m2;
    wbvector<unsigned> D1, D2, I1, I2;
    wbindex i1h, i2h, i1, i2;
    wbperm P;

    wbMatrix<TQ> Q1h, Q2h, Q1, Q2;
    wbvector<TD> E1, E2;
    wbarray<TD> Omega, CC;
    TD *e1, *e2;

    if (is>=ASpec.Ap.dim1) wblog(F,L,
       "ERR %s() index out of bounds (%d/%d)",FCT,is,ASpec.Ap.dim1);

    if (isbuf && isbuf[0]) {
       if (!strcmp(isbuf,"buf")) tobuf=1;
       else wblog(F,L,"ERR %s() invalid isbuf=>%s<",FCT,isbuf);
    }

    if (RX.isEmpty()) {
       wblog(F_L,"WRN %N%s() empty spectral update !??",FCT);
       return;
    }

    if (RX.rank(FL)==3 || !RX.CGR.isEmpty()) { RX.isConsistent(FL);
       if (RX.CGR.isEmpty()) {
          if (!RX.isAbelian()) { wblog(FL,
             "ERR %s() expecting CG data for rank-3 data (%s)",
             FCT, RX.qtype.toStr().data);
          }
       }
       else {
       for (i=0; i<RX.DATA.len; ++i) {
          wbarray<TD> &D=(*RX.DATA[i]);

          if (D.SIZE.len<2 || D.SIZE.len>3) wblog(FL,
             "ERR invalid rank operator (%s)",D.sizeStr().data);
          if (D.SIZE.len==3) {
             if (D.SIZE[2]!=1) wblog(FL,
                "ERR expecting rank-3 to have dim3=1 (%s)",D.sizeStr().data);
             else D.SIZE.len=2;
          }

          for (j=0; j<RX.CGR.dim2; ++j) {
             if (!RX.CGR(i,j).isScalar()) wblog(FL,
                "ERR %s() got non-scalar cgref (%d,%d): %s",
                FCT,i+1,j+1,RX.CGR(i,j).toStr().data
             );
          }
       }}
    }

    getEData(H1, E1, D1); D1.cumsum0(I1);
    getEData(H2, E2, D2); D2.cumsum0(I2);

    E1 *= Escale;
    E2 *= Escale;

    H1.getQsub(0,Q1h); RX.getQsub(0,Q1);
    H2.getQsub(0,Q2h); RX.getQsub(1,Q2);

    i=matchIndex(Q1h,Q1,i1h,i1, 1,&m1,&m2); if (m1 || i1.len!=Q1.dim1) ++e;
    i=matchIndex(Q2h,Q2,i2h,i2, 1,&m1,&m2); if (m1 || i2.len!=Q2.dim1) ++e;

    if (e) wblog(FL,
       "ERR %s() need unique and exact match (%d)!",FCT,e);
    if (i1.len!=RX.DATA.len || i1.len!=i2.len) wblog(FL,
       "ERR %s() severe QSpace inconsistency (%d/%d; %d)!",
        FCT, i1.len, i2.len, RX.DATA.len
    );

    i1.toPerm(P); P.Invert(); i1h.Select(P);
    i2.toPerm(P); P.Invert(); i2h.Select(P);

    for (k=0; k<i1.len; ++k) {
        i=i1h[k]; m=D1[i];
        j=i2h[k]; n=D2[j]; Omega.init(m,n);

        if (!RX.DATA[k]->hasSameSize(Omega)) {
           MXPut(FL,"i").add(*RX.DATA[k],"RX").add(Omega,"Omega").add(k+1,"k");
           RX.DATA[k]->info("RX[k]"); Omega.info("Omega");
           wblog(FL,"WRN %s() dimension mismatch (%d) !??",FCT,k);
        }

        e1=E1.data+I1[i]; e2=E2.data+I2[j];

        for (r=0; r<m; ++r)
        for (s=0; s<n; ++s) Omega(r,s) =  e2[s]-e1[r];

        if (tobuf) {
           ASpec.Add2buf(Omega, *(RX.DATA[k]), is); }
        else {
           ASpec.Add(
              Omega, *(RX.DATA[k]), is, 
              dblock && RX.hasSameQ(k,0,k,1) ? dblock : 0,
              mflag
           );
        }
    }
};


