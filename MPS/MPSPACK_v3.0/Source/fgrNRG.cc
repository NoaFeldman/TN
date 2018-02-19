char USAGE[] =
/* =====================================================================
 * fgrNRG.cc; (C) Wb, June 2008 -- tags: CCT */

"   Usage: [cc [,Info]] = ...                                           \n\
   fgrNRG(NRG0, NRG2, Lambda, C1 [,C2], tt, [,opts])                    \n\
                                                                        \n\
      full density matrix (FDM) Numerical Renormalization Group (NRG)   \n\
      calculation for Fermi Golden Rule rates.                          \n\
                                                                        \n\
   This routine calculates the time evolution of local operators        \n\
   C1(t)*C2 in fdm-NRG framework at given temperature (default T=TN)    \n\
                                                                        \n\
      cc(t) = < C1(t)C2'> = trace( rho_NRG0 * C1(t)C2')                 \n\
                                                                        \n\
   with the MIXED time evolution C(t) = exp(iH0*t) C exp(-iH*t)         \n\
   since the initial and the final Hamiltonians must be used            \n\
   to account for om=E_final-E_initial (!!)                             \n\
                                                                        \n\
   If C2 is not specified or specified as empty [], it will be          \n\
   automatically assumed equal to C1, as is usually the case            \n\
   for Fermi Golden Rule caluclations, while EOM technique leads        \n\
   to C1!=D2.                                                           \n\
                                                                        \n\
   This routine uses data from previous standard NRG runs with          \n\
   the result written as a matrix product state in QSpace format.       \n\
   It must contain the state space kept (AK) as well as the one         \n\
   thrown away (AT) at every NRG iteration.                             \n\
                                                                        \n\
   The inital state is determined by the density matrix of              \n\
   NRG_DATA1 at given temperature (default T=0K). The time              \n\
   evolution is determined by NRG_DATA2                                 \n\
                                                                        \n\
   See also: NRGWilsonQS and dmNRG                                      \n\
                                                                        \n\
   Input:                                                               \n\
                                                                        \n\
     NRG_DATA*  structure containing data of previous NRG run; if this  \n\
                contains path/file, then the data will be read from the \n\
                files path/file_##.mat for every site (see NRGWilsonQS) \n\
                                                                        \n\
      .AK       MPS state of k(ept) space (regular NRG, QSpace)         \n\
      .AT       MPS state of t(runcated) NRG space (QSpace)             \n\
      .HK       effective NRG eigenspectra (kept, QSpace)               \n\
      .HT       effective NRG eigenspectra (truncated, QSpace)          \n\
                                                                        \n\
     Lambda     NRG discretization parameter for the conduction band.   \n\
                Lambda is needed to undo the energy scaling permformed  \n\
                in obtaining the NRG data at the first place.           \n\
                                                                        \n\
     C[12]      local operators, e.g. acting on the impurity            \n\
                                                                        \n\
     NB! If operators C[12] are specified by name, they will be         \n\
     read from workspace and stored in NRG_DATA together with           \n\
     the other NRG data (can be disabled by flag `nostore').            \n\
                                                                        \n\
   Options / Flags                                                      \n\
                                                                        \n\
    'Eoffset',. extra offset of log. discretization bins (0.)           \n\
                                                                        \n\
    'alpha',..  exponential decay of energies exp(-alpha dE) [0.]       \n\
                this applies alpha in time domain (F. Anders 2005)      \n\
    'sigma',..  applies broadening in frequency domain (log.gauss)      \n\
                NB! either alpha or sigma can be specified but not both.\n\
                default sigma=0.3 (with alpha ignored)                  \n\
                                                                        \n\
    'NRho', k   builds density matrix from iteration k only             \n\
    'T',...     temperature at which correlation function is calculated (0.)\n\
                                                                        \n\
    'nlog',...  number of bins per decade for logarithmic discretization (256)\n\
    'emin',...  minimum energy for omega binning (TN/1000)              \n\
    'emax',...  maximum energy for omega binning (10.)                  \n\
                                                                        \n\
    'nostore'   do not store operator elements of C[12] with NRG_DATA   \n\
    'locRho'    if Rho needs to be calculated it is stored in RAM only  \n\
                and will be deleted after the program is finished       \n\
    'partial'   return partial data from every NRG iteration            \n\
    'raw'       raw data only, i.e. no smoothing of data                \n\
    'RAW'       calculates time dependent data from exact frequency data\n\
                i.e. without binning of energies on log. scale          \n\
    'disp'      more detaild output (0)                                 \n\
                                                                        \n\
   Z0, zflags:                                                          \n\
                                                                        \n\
     Z0         Z0 operator for first site (impurity) required for      \n\
                fermionic setup where Z0=(-1)^(number of particles)     \n\
    'zflags'    on which operators C[12] to apply fermionic signs Z0    \n\
                ex. 'zflags', [ 1 0 ... ]  would apply Z0 to the first  \n\
                but not to the second operator  (default: all).         \n\
                                                                        \n\
   Equivalently instead of Z0, zflags, the operators C[12] can be       \n\
   specified explicitely as cell array where each cell vector is        \n\
   interpreted as being contracted onto { L(1), sigma(1), sigma(2),...} \n\
   NB! for every fermionic SITE in C[12], Z must be already applied     \n\
   upon input! (as in Z*C; this encludes idtentity operators not        \n\
   explicitely specified in C!; this deprecates Z0 and zflags).         \n\
                                                                        \n\
   NB! The first (last) bin contains all data below (above) emin (emax) \n\
                                                                        \n\
   Output                                                               \n\
                                                                        \n\
     cc         expectation value at times specified in tt              \n\
     Info       info structure                                          \n\
                                                                        \n\
   fgrNRG('nrgdata0', 'nrgdata2', 'setup.mat', 'C1', 'C2')              \n\
                                                                        \n\
     This usage allows to read all input parameters from given MAT      \n\
     files (MatLab binaries); all variables and options MUST be         \n\
     defined in setup.mat by their names as indicated above;            \n\
     e.g. the operators C1 and C2 are specified by their name.          \n\
                                                                        \n\
     The output will be written / updated to the file set               \n\
     `nrgdata0_##.mat' and `nrgdata2_info.mat'.                         \n\
                                                                        \n\
     Output will be appended to setup.mat.                              \n\
                                                                        \n\
   See also NRGWilsonQS, fdmNRG_QS, tdmNRG                              \n\
   AWb (C) June 2008                                                    \n\
";

char USAGE_2[] = // Alternative usage
"   fgrNRG('nrgdata0', 'nrgdata2', 'setup.mat', 'C1', 'C2')             \n\
                                                                        \n\
     This usage allows to read all input parameters from given MAT      \n\
     files (MatLab binaries); all variables and options MUST be         \n\
     defined in setup.mat by their names as indicated above;            \n\
     e.g. the operators C1 and C2 are specified by their name.          \n\
                                                                        \n\
     The output will be written / updated to the file set               \n\
     `nrgdata0_##.mat' and `nrgdata2_info.mat'.                         \n\
                                                                        \n\
     Output will be appended to setup.mat.                              \n\
                                                                        \n\
   See also NRGWilsonQS, dmNRG, tdmNRG                                  \n\
                                                                        \n\
   AWb C Wb,Jun28,08                                                    \n\
";





#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "fgrNRG"
#endif

#define gTQ int
#define gTD double

#define  EPS 1E-12
#define DEPS DBL_EPSILON
#define NLEN 128

double nrgScale(int iter);

#include "wblib.h"
#include "tdmdata.hh"
#include "spectral.hh"
#include "nrgdata.hh"
#include "dmrho.cc"

template <class TQ, class TD>
void fgrNRGIter(
   NRGData<TQ,TD> &C1,
   NRGData<TQ,TD> &C2,
   QSpace<TQ,TD> &Rho, char wrho,
   NRGData<TQ,TD> &H0,
   NRGData<TQ,TD> &H,
   const double dEK,
   TDMData &fgrData,
   int iter
);

template <class TQ, class TD>
void updateSpec(
   const char *F, int L,
   TDMData &fgrData,
   unsigned iter, unsigned io,
   const QSpace<TQ,TD> &RX,
   const QSpace<TQ,TD> &H1,
   const QSpace<TQ,TD> &H2,
   const double Escale,
   const double dEK
);

#ifdef MATLAB_MEX_FILE
static void CleanUp(void) {
}
#endif

   double Lambda=1;

void nrgGetE0(const wbvector<double> &dE, wbvector<double> &E0) {
   double rL=1/sqrt(Lambda), dbl=rL;

   if (E0.len!=dE.len)
   E0.init(dE.len);

   E0[0]=dE[0];

   for (unsigned k=1; k<dE.len; k++, dbl*=rL)
   E0[k]=E0[k-1]+dE[k]*dbl;

   E0*=nrgScale(0);
}


void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   unsigned i,k,m,nloc=0,iter, e=0, N=0, QDIM, NRho=0, dloc=0, nlog=256;
   char disp=0, store=1, partial=0, RAW=0, raw=0, gotC2=0;
   char calcRho=0, locRho=0, calcOps=0;
   double TN,dbl, emin, dEK, Eoffset=0, emax=10., wx=0., T=-1;
   double alpha=0, sigma=0;

   char vtag[4], itag[16];
   mxArray *S, *a, *ar=NULL;
   const mxArray *ac, *bc=NULL;

   wbMatrix<double> a0,aa;
   wbMatrix<wbcomplex> at;
   TDMData fgrData;

   wbstring vstr(64), istr(64), fnam0, fname;

   NRGData<gTQ,gTD> A0("A"), H0("H"), A("A"), H("H"), RHO("RHO");
   NRGData<gTQ,gTD> C1("C1"), C2("C2");
   wbvector< wbvector< QSpace<gTQ,gTD> > > CC,DD;
   wbvector< QSpace<gTQ,gTD> > C_,D_,Rho;
   QSpace<gTQ,gTD> Z0;

   wbMatrix<gTD> rhoNorm, dd;
   wbvector<double> om,D0,D2,E0,E2,w,t0;
   wbvector<char> zz;
   OPTS opts;

   wbSigHandler SIG(FL);

   if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
   NRG_N=0; NRG_ITER=0; STRICT_ITER0=0; gES.init(); 

#ifdef WB_CLOCK
   WbClock wbc_fgrNRG(PROG);
   wbc_fgrNRG.start();
#endif
#ifdef MATLAB_MEX_FILE
   mexAtExit(CleanUp);
#endif


#ifndef MAIN

   if (nargin<5 || nargout>2) wblog(FL,
   "ERR Invalid number of I/O arguments (%d,%d).",nargin,nargout);

   for (k=0; k<2; k++) if (mxIsChar(argin[k])) {
      i=mxGetString(argin[k],str,NLEN);
      if (i || !str[0]) wblog(FL,
         "ERR Error reading file tag from arg #%d (%d,%d)",
          k+1, i, strlen(str));
      if (k==0) fnam0=str; else fname=str;
   }

   if (!mxIsDblScalar(argin[k])) wblog(FL,
   "Invalid argument %d (Lambda).",k+1);
   if (mxGetNumber(argin[k++], Lambda)) wberror(FL,str);


   gotC2 = nargin>int(k+2) &&
    ((mxIsStruct(argin[k+1]) || mxIsCell(argin[k+1])));

   ac=argin[k++]; if (gotC2) {
   bc=argin[k  ]; }

   if (ac && mxIsCell(ac)) {
      Wb::initLocalOp(FL,ac,CC,k,  '!'); if (gotC2) {
      Wb::initLocalOp(FL,bc,DD,k+1,'!'); }
   }
   else {
      Wb::initLocalOp(FL,ac,C_,k,  '!'); if (gotC2) {
      Wb::initLocalOp(FL,bc,D_,k+1,'!'); }
   }

   if (gotC2) k++;

   t0.init(FL,argin[k++]);

   opts.init(argin+k, nargin-k);

#else


   if (nargout || nargin<4 || nargin>5) wblog(FL,
   "ERR Invalid number of I/0 arguments (%d,%d).", nargin, nargout);

   for (i=0; i<(unsigned)nargin; i++) if (!mxIsChar(argin[i]))
   wblog(FL,"ERR Invalid arg #%d %N%N%s%N", i+1, USAGE_2);

   for (k=0; k<2; k++) {
      i=mxGetString(argin[k],str,NLEN);
      if (i || !str[0]) wblog(FL,
         "ERR Error reading file tag from arg #%d (%d,%d)",
          k+1, i, strlen(str));
      if (k==0) fnam0=str; else fname=str;
   }

   if ((i=mxGetString(argin[k++],str,NLEN))) wblog(FL,
   "ERR Error reading string of arg #%d (%d,%d)", k+1,i,strlen(str));

   opts.init(str);"options" by reference to file

   opts.getOpt("tt",a,'!'); t0.init(FL,a);
   opts.getOpt("Lambda", Lambda,'!');

   gotC2=nargin>4;

   mxArray *b=NULL;

   if (mxGetString(argin[k++],str,MAX_STRLEN))
   wblog(FL,"ERR failed to read string of arg #%d", k);
   C1.name=str; opts.getOpt(str,a,'!'); if (gotC2) {

   if (mxGetString(argin[k  ],str,MAX_STRLEN))
   wblog(FL,"ERR failed to read string of arg #%d", k+1);
   C2.name=str; opts.getOpt(str,b,'!'); }

   if (a && mxIsCell(a)) {
      Wb::initLocalOp(FL,a,CC,k,  '!'); if (gotC2) {
      Wb::initLocalOp(FL,b,DD,k+1,'!'); }
   }
   else {
      Wb::initLocalOp(FL,a,C_,k,  '!'); if (gotC2) {
      Wb::initLocalOp(FL,b,D_,k+1,'!'); }
   }

   if (gotC2) k++;

#endif

   if ((fnam0.isEmpty() || fname.isEmpty()) && nargout<2) wblog(FL,
   "ERR %s requires at least two output arguments", PROG);


   A0.setupIO(FL, fnam0.data, argin[0], 'i');
   H0.setupIO(FL, fnam0.data, argin[0]);
   RHO.setupIO(FL,fnam0.data, argin[0]);

   A .setupIO(FL, fname.data, argin[1], 'i');
   H .setupIO(FL, fname.data, argin[1]);
   C1.setupIO(FL, fname.data, argin[1]);
   C2.setupIO(FL, fname.data, argin[1]);


   if (Lambda<=1.) wblog(FL,"ERR Invalid Lambda=%g.",Lambda);

   H0.checkVec(FL,"K",2,"NRG data (HK0)");
   H .checkVec(FL,"K",2,"NRG data (HK)" );

   H0.checkVec(FL,"T",2,"NRG data (HT0)");
   H .checkVec(FL,"T",2,"NRG data (HT)");
   A0.checkVec(FL,"T",3,"NRG data (AT0)");
   A .checkVec(FL,"T",3,"NRG data (AT)");

   if (NRG_N==0) wblog(FL, "ERR Input NRG data is empty.");
   N=NRG_N; i=0;

   ac=H0.getMxInfo("E0"); if (ac) D0.init(FL,ac); else i=1;
   ac=H .getMxInfo("E0"); if (ac) D2.init(FL,ac); else i=2;
   if (i) wblog(FL,"ERR failed to retrieve E0 from NRG data (%d)",i);

   if (N!=D0.len || N!=D2.len) wblog(FL,
      "ERR Mismatch in chain length (N=%d, E0: %d %d)",N,D0.len,D2.len);

   ac=RHO.getMxInfo("EScale");
   if (ac==NULL) wblog(FL,"ERR failed to read EScale (info.mat)");
   else {
      gES.init(FL,ac);
      if (N!=gES.len) wblog(FL,
         "ERR mismatch in chain length (EScale: %d/%d)",gES.len,N);
      else wblog(FL," *  using EScale from info.mat");
   }

   dloc=getDLoc(A);

   A.init(FL,"K",0); QDIM=A.K.QDIM;

   TN=nrgScale(N-1); emin=TN/1000.;

#ifndef MATLAB_MEX_FILE
   if (Wb::GetEnv(FL,"DMA_T",T)>=0)
#endif


   opts.getOpt("T", T);
   partial=opts.getOpt("partial");
   raw=opts.getOpt(FL,"raw");
   RAW=opts.getOpt(FL,"RAW");
   opts.getOpt(FL, "alpha",alpha); sigma=(alpha>0 ? 0 : 0.3);
   opts.getOpt(FL, "sigma",sigma);

   opts.getOpt(FL,"Eoffset",Eoffset);

   if (!RAW && !raw && alpha<=0 && sigma>0) {
      alpha=1;
   }

   store =!opts.getOpt("nostore");
   locRho =opts.getOpt("locRho" );
   calcRho=opts.getOpt("calcRho");
   calcOps=opts.getOpt("calcOps");

   opts.getOpt("Z0",a); if (a) {
      if (!mxIsQSpace(a)) wblog(FL,"ERR invalid Z0");
      Z0.init(FL,a);
   }

   opts.getOpt("disp", disp);
   opts.getOpt("nlog", nlog);
   opts.getOpt("emin", emin);
   opts.getOpt("emax", emax);
   opts.getOpt("NRho", NRho);
   
   opts.getOpt("zflags",a); if (a) zz.init(FL,a,"zflags");

   if (NRho>N) wblog(FL,"ERR NRho out of bounds (%d/%d)", NRho, N);

   opts.checkAnyLeft();

   if (NRho==0) {
      dbl=sqrt(2/(log(Lambda)*log(double(dloc))));
      m=(unsigned)ceil(5*dbl+3);

      dbl=nrgScale(N-m-1);
      if (T<0) {
         T=dbl; wblog(FL," *  T=%.4g (=T[N-%d])", T, m);
      }
      else if (T>0 && T<dbl)
      wblog(FL,"WRN T=%.3g < T_{N-%d}=%.3g", T, m, dbl);
   }
   else if (T<0) T=0;

   if (NRho) { vstr.printf(FL,
      "tDM-NRG (single shell NRho=%d for Fermi's golden rule)", NRho);
       strcpy(vtag,"fgK");
   }
   else {
       vstr.cpy(FL,"tDM-NRG (Fermi's golden rule)");
       strcpy(vtag,"FGR");
   }

   str[0]=0;
   if (disp) sprintf(str," disp=%d", disp);
   if (!store ) strcat(str, " nostore");
   if (calcOps) strcat(str, " calcOps");
   if (locRho ) strcat(str, " locRho" );
   if (calcRho) strcat(str, " calcRho");

   wblog(FL,"<i> %47R\n%-12s: %s\n"
       "parameters  : N=%d, QDIM=%d, d=%d%s%s\n"
       "temperature : %.4g (= %.4g TN)\n"
       "fgrData     : %.2g .. %.2g (%d/dec)",
   "=", PROG, vstr.data, N, QDIM, dloc, str[0] ? "," : "", str[0] ? str : "",
   T, T/TN, emin, emax, nlog);

   if (t0.len) wblog(FL,
      "<i> time range  : %.2g .. %.2g (%d vals)",
      t0[0], t0.end(), t0.len
   );

   wblog(FL,"<i> %47R","=");


   if (!CC.isEmpty()) { nloc=0;

      if (zz.len || !Z0.isEmpty()) wblog(FL,
         "WRN %Nsince C is specified as cell array,\n"
         "zflags and Z0 will be ignored!%N"
      );

      if (!DD.len) {
         for (k=i=0; i<CC.len; i++) {
            nloc=MAX(nloc,CC[i].len);
            if (CC[i].len) { if (k<i) CC[k]=CC[i]; k++; }
         }
         if (k<i) { wblog(FL,
            "WRN empty operator sets will be ignored (%d/%d)",i-k,i);
            CC.len=k;
         }
      }
      else {
         if (DD.len!=CC.len) wblog(FL,"ERR C2 since specified "
         "must have the same length as C1 (%d/%d)", DD.len, CC.len);

         for (k=i=0; i<CC.len; i++) {
            nloc=MAX(nloc, MAX(CC[i].len, DD[i].len) );
            if (CC[i].len && DD[i].len) {
               if (k<i) { CC[k]=CC[i]; DD[k]=DD[i]; }
               k++;
            }
         }
         if (k<i) { wblog(FL,
            "WRN empty operator sets will be ignored (%d/%d)",i-k,i);
            CC.len=k; DD.len=k;
         }
      }

      C1.initOp(FL,CC, NULL,"1st");
      C2.initOp(FL,DD, &C1, "2nd");
   }
   else {

      if (D_.len && D_.len!=C_.len) wblog(FL,"ERR C2 since specified "
      "must have the same length as C1 (%d/%d)", D_.len, C_.len);

      C1.initOp(FL,C_); C1.initSym("KK");
      C2.initOp(FL,D_); C2.initSym("KK");

      if (calcOps) { C1.forceCalc(); C2.forceCalc(); }


      if (!zz.data) { zz.init(C_.len); zz.set(1); }
      else {
         if (zz.len!=C_.len) {
            if (zz.len==1) { zz.Resize(C_.len); zz.set(zz[0]); }
            else wblog(FL,
               "ERR Length mismatch of zflags with C (%d/%d)",
                zz.len, C_.len
            );
         }
      }

      if (zz.anyUnequal(0)) { e=0;
         wblog(FL," *  applying Z0 [zflags = %s]", zz.toStr().data);

         if (Z0.isEmpty()) wblog(FL,
            "ERR Z0 not defined (required since zflags is specified!)");
         if ((C1.KK.len && C1.KK[0].rank()!=Z0.rank()) ||
             (C2.KK.len && C2.KK[0].rank()!=Z0.rank()) ) wblog(FL,
            "ERR rank mismatch of Z0 with ops (%d/%d) !??",
             Z0.rank(),C1.KK[0].rank());

         if (C1.calc>0) {
            for (i=0; i<C1.KK.len; i++) if (zz[i])
            e+=applyZ0(Z0,C1.KK[i]);
         }

         if (C2.calc>0) {
            for (i=0; i<C2.KK.len; i++) if (zz[i])
            e+=applyZ0(Z0,C2.KK[i]);
         }

         if (e) wblog(FL,
            "ERR QIDX of Z0 does not match QIDX of %s%s%s\n"
            "Hint: unset Z0 to empty if not relevant.", C1.name.data,
            C2.KK.len ? " or " : "",
            C2.KK.len ? C2.name.data : ""
         );
      }
   }

   if (C1.KK.isEmpty()) wblog(FL,
   "ERR first operator (C1) must be set (yet is empty)");

   if (calcOps) {
      C1.forceCalc();
      C2.forceCalc();
   }
   if (!store) {
      C1.store=0;
      if (C2.nrgIdx.len) C2.store=0;
   }

   if (C2.KK.isEmpty()) { C1.name="C";
      C1.info(FL,"C");
   }
   else {
      C1.info(FL,"C1");
      C2.info(FL,"C2");
   }

   if (1 || C1.store) {
      if (CC.isEmpty())
           C1.updateInfO(FL,"i",C_.toMx());
      else C1.updateInfO(FL,"i",C1.CI.toMx());
   }
   if (1 || (C2.store && C2.CI.len)) {
      if (DD.isEmpty())
           C2.updateInfO(FL,"i",D_.toMx());
      else C2.updateInfO(FL,"i",C2.CI.toMx());
   }

   if (!C1.calc && !C2.calc) wblog(FL,
   "<i> using stored operator sets");

   fgrData.init(FL,"FGR", t0,
      N, C1.KK.len, emin, emax, nlog
   );

   fgrData.checkAvgOm();

   if (RAW) fgrData.raw=1;


   nrgGetE0(D0,E0);
   nrgGetE0(D2,E2);

   dbl=E0[NRho>0 ? NRho-1 : N-1];
   E0-=dbl; initRHO(rhoNorm, T, H0, E0, dloc, NRho);
   E2-=dbl;

   dbl=(E2.end()-E0.end());
   if (dbl!=0 || Eoffset!=0) {
   if (dbl==0 || Eoffset==0)
        wblog(FL,"<i> Eoffset = %.8g",Eoffset+dbl);
   else wblog(FL,"<i> Eoffset = %.8g %+.8g = %.8g",Eoffset,dbl,Eoffset+dbl); }
   Eoffset += dbl;

   dd.init(RHO.getMxInfo("rhoNorm"));

   if (calcRho || dd!=rhoNorm || RHO.checkVec(FL,"",2) ||
       (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T)
   ){
      wblog(FL,"<i> update density matrices\\");

      if (locRho && !RHO.MX) {
         printf(" internally (locRho=%d)",locRho);
         RHO.MX=mxCreateStructMatrix(1,NRG_N,0,NULL);
         RHO.MXloc=1; locRho=99;
      }

      printf("\n\n");
      nrgUpdateRHO(rhoNorm, A0,H0, RHO, NRho, &SIG);
      printf("\n");
   }
   else wblog(FL,"<i> using stored DM set.%N");

   if (locRho!=99) locRho=0;


   Rho.initDef(2);

   for (iter=0; iter<N; iter++) {
      NRG_ITER=iter; SIG.call99();

      rhoNorm.getRec(iter,w);
      wx+=(fabs(w[0])+fabs(w[1]));
      
      A0.init(FL,"K",iter);  H0.init(FL,"K",iter);
      A0.init(FL,"T",iter);  H0.init(FL,"T",iter);
      A .init(FL,"K",iter);  H .init(FL,"K",iter);
      A .init(FL,"T",iter);  H .init(FL,"T",iter);

      snprintf(itag,16,"%s %2d/%d", vtag, iter+1, N);

      C1.updateOp(FL,A0,A,4,iter); C1.storeOp(FL,iter,4);
      C2.updateOp(FL,A0,A,4,iter); C2.storeOp(FL,iter,4);

      dEK=Eoffset+E0[iter]-E2[iter];

      if (iter+2<=nloc && !A.T.isEmpty()) wblog(FL,
         "WRN %Nlocal operator extends to trancated sites (%d/%d)",
          iter+1, nloc-1
      );

      if (wx>DEPS) {
         if (w[0]>DEPS) {
            if (w[1]>DEPS) wblog(FL,
            "%s kept (%.3g) + truncated (%.3g)",itag,w[0],w[1]);
            else wblog(FL, "%s kept (%.3g)",itag,w[0]);
         }
         else if (w[1]>DEPS)
         wblog(FL,"%s truncated (%.3g)%20s\r\\",itag,w[1],"");

         for (i=0; i<2; i++) if (w[i]) {
            const QSpace<gTQ,gTD> &HX = (i==0 ? H0.K : H0.T);

            initRho(Rho[i], HX, iter, w[i]);
            fgrNRGIter(
               C1, C2, Rho[i],(i==0 ? 'K':'T'), H0, H, dEK,
               fgrData, iter
            );
         }
      }
      else {
         if (disp>2) wblog(FL,"%s (RHO)",itag);
         else wblog(FL,"%s (RHO)%30s\r\\",itag,"");
      }

      RHO.init(FL,"A",iter);
      fgrNRGIter(
         C1, C2, RHO.A,'K', H0, H, dEK,
         fgrData, iter
      );
   }

   if (disp<=2) printf("\n\n");


   if (RAW) {
      fgrData.addPartialFourier();
      ar=fgrData.At.toMx('r');
   }

   fgrData.getSpecData(om,a0);

   if (raw) {}
   else if (sigma>0 && t0.len) {
      istr.printf(FL,"frequency broadened TDM data (sigma=%g)", sigma);

      fgrData.getSmoothSpec(sigma);
      fgrData.Fourier(om,aa,t0,at,alpha);
   }
   else if (alpha>0) {
      if (Lambda<=1) wblog(FL,"ERR alpha requires Lambda (%g)",Lambda);
      istr.printf(FL,"time-domain broadened TDM data (alpha=%g)",alpha);

      fgrData.Fourier(om,at,alpha,Lambda);
   }
   else if (t0.len) {
      istr.printf(FL,"plain Fourier transformed data (no broadening)");
      wblog(FL,"<i> %s",istr.data);

      Wb::Fourier(om,a0,t0,at);
   }
   
   S=mxCreateStructMatrix(1,1,0,NULL);
   mxAddField2Scalar(FL,S,"ver",    vstr.toMx());


   mxAddField2Scalar(FL,S,"T",      numtoMx(T));
   mxAddField2Scalar(FL,S,"alpha",  numtoMx(alpha));
   mxAddField2Scalar(FL,S,"sigma",  numtoMx(sigma));
   mxAddField2Scalar(FL,S,"rho",    rhoNorm.toMx());
   mxAddField2Scalar(FL,S,"E0",     E0.toMx());
   mxAddField2Scalar(FL,S,"E2",     E2.toMx());
   mxAddField2Scalar(FL,S,"EScale", gES.toMx());
   mxAddField2Scalar(FL,S,"Eoffset",numtoMx(Eoffset));
   mxAddField2Scalar(FL,S,"avgOm",  fgrData.getAvgOm().toMx());

   mxAddField2Scalar(FL,S,"om",     om.toMx('t'));
   mxAddField2Scalar(FL,S,"a0",     a0.toMx('r'));
   mxAddField2Scalar(FL,S,"aa",     aa.toMx('r'));

   if (ar)
   mxAddField2Scalar(FL,S,"cc_raw", ar);

   if (partial)
   mxAddField2Scalar(FL,S,"RAW",    fgrData.toMx());
   mxAddField2Scalar(FL,S,"emin",   numtoMx(emin));
   mxAddField2Scalar(FL,S,"emax",   numtoMx(emax));
   mxAddField2Scalar(FL,S,"nlog",   numtoMx(nlog));

   { mxArray *F=mxCreateCellMatrix(1,2);
     mxSetCell(F,0,H0.NAME.toMx());
     mxSetCell(F,1,C1.NAME.toMx());
     mxAddField2Scalar(FL,S,"NRG",F);
   }

   mxAddField2Scalar(FL,S,"stamp",  wbTimeStamp().toMx());

   if (!C1.NAME.isEmpty()) {
      C1.updateInfo(FL, "tt",   t0.toMx('t'));
      C1.updateInfo(FL, "cc",   at.toMx('t'));
      C1.updateInfo(FL, "Ifgr", S, 'k');
   }

   if (!RHO.NAME.isEmpty()) {
      RHO.updateInfo(FL, "rhoT",    numtoMx(T));
      RHO.updateInfo(FL, "rhoNorm", rhoNorm.toMx());
   }

   if (nargout)
   argout[0]=at.toMx('t');

   if (nargout>1)
        argout[1]=S;
   else mxDestroyArray(S);

#ifdef WB_CLOCK
   wbc_fgrNRG.stop();
#endif
#ifdef MATLAB_MEX_FILE
   CleanUp();
#endif

   wblog(FL,"I/O NRGWilsonQS time usage: %s%N",
   tic_nrgIO.gettime("asm").data);
   tic_nrgIO.reset();


}



template <class TQ, class TD>
void fgrNRGIter(
   NRGData<TQ,TD> &C1,
   NRGData<TQ,TD> &C2,
   QSpace<TQ,TD> &Rho, char w,
   NRGData<TQ,TD> &H0,
   NRGData<TQ,TD> &H,
   const double dEK,
   TDMData &fgrData,
   int iter
){
   unsigned i,k,n,m=0;
   QSpace<TQ,TD> R2;
   char tags[][3] = { "*K","*T" };
   const double Escale = nrgScale(iter);

   if (iter<0) iter=NRG_ITER;

   n=C1.KK.len;
   if (n!=C1.KT.len || n!=C1.TK.len || n!=C1.TT.len) wblog(FL,
     "ERR severe operator dim mismatch [%d %d %d %d]",
      C1.KK.len, C1.KT.len, C1.TK.len, C1.TT.len
   );
   if (C2.KK.len)
   if (n!=C2.KT.len || n!=C2.TK.len || n!=C2.TT.len) wblog(FL,
     "ERR severe operator dim mismatch [%d %d %d %d]",
      C2.KK.len, C2.KT.len, C2.TK.len, C2.TT.len
   );

   if (w!='K' && w!='T') wblog(FL,"ERR invalid w=%c<%d> !??",w,w);
   for (i=0; i<2; i++) tags[i][0]=w;

   if (w=='K') m=1;

   for (i=m; i<2; i++) {
      const QSpace<TQ,TD>
         &H1 = H0.getQSpace(tags[i][0]),
         &H2 = H .getQSpace(tags[i][1]);

      for (k=0; k<n; k++) {
         const QSpace<TQ,TD>
          &C1k = C1.getQSpace(tags[i],k),
          &C2k = (C2.KK.len ? C2.getQSpace(tags[i],k) : C1k);

         Rho.contract(2,C2k,1,R2);


         R2.TimesEl(C1k);

         if (R2.isEmpty()) continue;

         updateSpec(FL,
            fgrData, iter, k,
            R2, H1, H2, Escale, dEK
         );
      }
   }
}


template <class TQ, class TD>
void updateSpec(
   const char *F, int L,
   TDMData &fgrData,
   unsigned iter, unsigned io,
   const QSpace<TQ,TD> &RX,
   const QSpace<TQ,TD> &H1,
   const QSpace<TQ,TD> &H2,
   const double Escale,
   const double dEK
){
   unsigned i,j,k,m,n,r,s, e=0;
   wbvector<unsigned> D1,D2,I1,I2;
   wbindex i1h,i2h,i1,i2;
   wbperm P;

   wbMatrix<TQ> Q1h,Q2h,Q1,Q2;
   wbvector<TD> E1,E2;
   wbarray<TD> EE;
   TD *e1,*e2;

   getEData(H1, E1, D1); D1.cumsum0(I1);
   getEData(H2, E2, D2); D2.cumsum0(I2);

   E1 *= Escale;
   E2 *= Escale; E2 -= dEK;

   if (RX.isEmpty()) {
      wblog(F,L,"ERR %N empty spectral update ???");
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

   H1.getQsub(0,Q1h); RX.getQsub(0,Q1);
   H2.getQsub(0,Q2h); RX.getQsub(1,Q2);

   i=matchIndex(Q1h,Q1,i1h,i1); if (i || i1.len!=Q1.dim1) {
      wblog(FL,"WRN IDX-1: %d  (H: %d/%d)  %d/%d",
      i,i1h.len,Q1h.dim1, i1.len,Q1.dim1); e++;
   }
   i=matchIndex(Q2h,Q2,i2h,i2); if (i || i2.len!=Q2.dim1) {
      wblog(FL,"WRN IDX-2: %d  (H: %d/%d)  %d/%d",
      i,i2h.len,Q2h.dim1, i2.len,Q2.dim1); e++;
   }

   if (e) wblog(F,L,
      "ERR updateSpec - need unique and exact match (%d)!",e);
   if (i1.len!=RX.DATA.len) wblog(F,L,
      "ERR updateSpec - severe QSpace inconsistency (%d %d)!",
       i1.len, RX.DATA.len);

   i1.toPerm(P); P.Invert(); i1h.Select(P);
   i2.toPerm(P); P.Invert(); i2h.Select(P);

   for (k=0; k<i1.len; k++) {
       i=i1h[k]; m=D1[i];
       j=i2h[k]; n=D2[j]; EE.init(m,n);

       if (!RX.DATA[k]->hasSameSize(EE)) { 
          RX.DATA[k]->info("RX[k]"); EE.info("EE");
          wblog(FL,"WRN dimension mismatch!?? (%d)",k);
       }

       e1=E1.data+I1[i]; e2=E2.data+I2[j];

       for (r=0; r<m; r++)
       for (s=0; s<n; s++) EE(r,s)=e2[s]-e1[r];

       fgrData.Add(EE, *(RX.DATA[k]), iter,io);
   }
}


